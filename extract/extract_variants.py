import csv
import dxpy
import pandas as pd
import pandas.core.series

from typing import List
from pathlib import Path
from os.path import exists

from extract.extract_association_pack import ExtractAssociationPack
from general_utilities.association_resources import process_snp_or_gene_tar, build_transcript_table, get_gene_id, \
    run_cmd
from general_utilities.linear_model import linear_model
from general_utilities.linear_model.linear_model import LinearModelResult
from general_utilities.linear_model.proccess_model_output import merge_glm_staar_runs, process_linear_model_outputs, \
    process_staar_outputs
from general_utilities.linear_model.staar_model import staar_null, staar_genes
from general_utilities.job_management.thread_utility import ThreadUtility


class ExtractVariants:

    def __init__(self, output_prefix: str, association_pack: ExtractAssociationPack):

        self._outputs = []
        self._output_prefix = output_prefix
        self._association_pack = association_pack

        # Define our gene-list and make 'gene_info' objects of them (which are Pandas series classes)
        self._gene_infos = []
        self._chromosomes = set()
        # If we are doing extraction based on individual SNPs or a Gene list, we need to make a 'fake' gene info
        # but find all chromosomes those SNPS/Genes lie on
        if self._association_pack.is_non_standard_tar:
            gene_info, returned_chromosomes = process_snp_or_gene_tar(self._association_pack.is_snp_tar,
                                                                      self._association_pack.is_gene_tar,
                                                                      self._association_pack.tarball_prefixes[0])
            self._gene_infos.append(gene_info)
            self._chromosomes = returned_chromosomes
        else:
            for gene in self._association_pack.gene_ids:
                transcripts_table = build_transcript_table()
                gene_info = get_gene_id(gene, transcripts_table)
                self._gene_infos.append(gene_info)
                self._chromosomes.add(gene_info['chrom'])

    def run_tool(self):

        # 1. Download variant VEP annotations
        print("Loading VEP annotations...")
        thread_utility = ThreadUtility(self._association_pack.threads,
                                       error_message='An extraction thread failed',
                                       incrementor=5,
                                       thread_factor=4)
        for chromosome in self._chromosomes:
            thread_utility.launch_job(self._download_vep,
                                      chromosome=chromosome)
        thread_utility.collect_futures()

        # 2. Filter relevant files to individuals we want to keep
        print("Filtering variant files to appropriate individuals...")
        thread_utility = ThreadUtility(self._association_pack.threads,
                                       error_message='An extraction thread failed',
                                       incrementor=20,
                                       thread_factor=4)

        # if, elif, else simply depends on which type of tarball we are using
        if self._association_pack.is_snp_tar:
            for tarball_prefix in self._association_pack.tarball_prefixes:
                thread_utility.launch_job(self._filter_individuals,
                                          tarball_prefix=tarball_prefix,
                                          chromosome='SNP')
        elif self._association_pack.is_gene_tar:
            for tarball_prefix in self._association_pack.tarball_prefixes:
                thread_utility.launch_job(self._filter_individuals,
                                          tarball_prefix=tarball_prefix,
                                          chromosome='GENE')
        else:
            for chromosome in self._chromosomes:
                for tarball_prefix in self._association_pack.tarball_prefixes:
                    thread_utility.launch_job(self._filter_individuals,
                                              tarball_prefix=tarball_prefix,
                                              chromosome=chromosome)
        thread_utility.collect_futures()

        # 4. Actually collect variant information per-gene
        print("Extracting variant information...")
        thread_utility = ThreadUtility(self._association_pack.threads,
                                       error_message='An extraction thread failed',
                                       incrementor=20,
                                       thread_factor=2)
        for gene_info in self._gene_infos:
            for tarball_prefix in self._association_pack.tarball_prefixes:
                thread_utility.launch_job(self._annotate_variants,
                                          tarball_prefix=tarball_prefix,
                                          gene_info=gene_info,
                                          chromosomes=self._chromosomes if
                                          self._association_pack.is_non_standard_tar else None)

        for result in thread_utility:
            self._outputs.extend(result)

        # 5. And run a linear and STAAR model(s) for all genes
        print("Running linear models...")
        self._run_linear_models()
        self._run_staar_models()

        # This function merges STAAR and GLM results together
        self._outputs.extend(merge_glm_staar_runs(self._output_prefix,
                                                  self._association_pack.is_snp_tar,
                                                  self._association_pack.is_gene_tar))

        # 6. Finally, add the phenotypes/covariates table to the outputs
        out_pheno_path = Path(f'{self._output_prefix}.phenotypes_covariates.formatted.tsv')
        Path('phenotypes_covariates.formatted.txt').rename(out_pheno_path)
        self._outputs.append(out_pheno_path)

    def get_outputs(self) -> List[Path]:
        return self._outputs

    def _download_vep(self, chromosome: str) -> None:

        vep_dx = self._association_pack.bgen_dict[chromosome]['vep']
        dxpy.download_dxfile(vep_dx.get_id(), chromosome + ".filtered.vep.tsv.gz")

    @staticmethod
    def _filter_individuals(tarball_prefix: str, chromosome: str) -> None:
        # And filter the relevant SAIGE file to just the individuals we want so we can get actual MAC
        cmd = f'bcftools view --threads 4 -S /test/SAMPLES_Include.txt -Ob ' \
              f'-o /test/{tarball_prefix}.{chromosome}.saige_input.bcf ' \
              f'/test/{tarball_prefix}.{chromosome}.SAIGE.bcf'
        run_cmd(cmd, is_docker=True, docker_image='egardner413/mrcepid-burdentesting')

    def _annotate_variants(self, tarball_prefix: str, gene_info: pandas.core.series.Series,
                           chromosomes: set) -> List[Path]:

        # This is a bit confusing, so explaining in full.
        # We need to annotate EACH GENE separately EXCEPT when running a SNP/GENE list tarball, SO...
        # 1. If we have a SNP/GENE tar, we may need to load in variant annotations for multiple chromosomes,
        # so we set 'chromosomes'
        # to a set of chromosomes that we extract from the SNP/GENE tar.
        # 2. If just a single gene or gene list (chromosomes = None), only need to load the data for the chromosome
        # that specific Gene is on
        if chromosomes is None:
            chromosome = gene_info['chrom']
            variant_index = pd.read_csv(f'{chromosome}.filtered.vep.tsv.gz',
                                        sep="\t",
                                        dtype={'SIFT': str, 'POLYPHEN': str})
        else:
            variant_index = []
            for chromosome in chromosomes:
                variant_index.append(pd.read_csv(f'{chromosome}.filtered.vep.tsv.gz',
                                                 sep="\t",
                                                 dtype={'SIFT': str, 'POLYPHEN': str}))
            variant_index = pd.concat(variant_index)

        # Need to get the variants from the SAIGE groupfile:
        with Path(f'{tarball_prefix}.{gene_info["chrom"]}.SAIGE.groupFile.txt').open('r') as saige_group_file,\
                Path(f'{tarball_prefix}.{gene_info["SYMBOL"]}.variants.txt').open('w') as var_file:
            var_ids = []
            for line in saige_group_file:
                data = line.rstrip().split("\t")
                if data[0] == gene_info.name:
                    for i in range(1, len(data)):
                        curr_id = data[i].replace('_', ':').replace('/', ':')
                        var_file.write(curr_id + "\n")
                        var_ids.append(curr_id)
                    break

        relevant_vars = variant_index[variant_index['varID'].isin(var_ids)]

        # Filter to the variants for this gene
        cmd = f'bcftools view --threads 2 -i \'ID=@/test/{tarball_prefix}.{gene_info["SYMBOL"]}.variants.txt\' -Ob ' \
              f'-o /test/{tarball_prefix}.{gene_info["SYMBOL"]}.variant_filtered.bcf ' \
              f'/test/{tarball_prefix}.{gene_info["chrom"]}.saige_input.bcf'
        run_cmd(cmd, is_docker=True, docker_image='egardner413/mrcepid-burdentesting')
        cmd = f'bcftools +fill-tags --threads 4 -Ob ' \
              f'-o /test/{tarball_prefix}.{gene_info["SYMBOL"]}.final.bcf ' \
              f'/test/{tarball_prefix}.{gene_info["SYMBOL"]}.variant_filtered.bcf'
        run_cmd(cmd, is_docker=True, docker_image='egardner413/mrcepid-burdentesting')

        # Now get actual annotations back in:
        cmd = f'bcftools query -f \'%ID\\t%MAF\\t%AC\\t%AC_Het\\t%AC_Hom\\n\' ' \
              f'-o /test/{tarball_prefix}.{gene_info["SYMBOL"]}.annotated_vars.txt ' \
              f'/test/{tarball_prefix}.{gene_info["SYMBOL"]}.final.bcf'
        run_cmd(cmd, is_docker=True, docker_image='egardner413/mrcepid-burdentesting')

        # And get a list of individuals with a variant:
        cmd = f'bcftools query -i \"GT=\'alt\'\" -f \'[%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%SAMPLE\\t%GT\n]\' ' \
              f'-o /test/{tarball_prefix}.{gene_info["SYMBOL"]}.carriers.txt ' \
              f'/test/{tarball_prefix}.{gene_info["SYMBOL"]}.final.bcf'
        run_cmd(cmd, is_docker=True, docker_image='egardner413/mrcepid-burdentesting')

        geno_table = pd.read_csv(tarball_prefix + "." + gene_info['SYMBOL'] + ".annotated_vars.txt",
                                 sep="\t",
                                 names=['varID', 'MAF_tested', 'AC_tested', 'AC_tested_Het', 'AC_tested_Hom'])
        geno_table = pd.merge(relevant_vars, geno_table, on='varID', how="left")

        carriers_table = pd.read_csv(tarball_prefix + "." + gene_info['SYMBOL'] + ".carriers.txt",
                                     sep="\t",
                                     names=['CHROM', 'POS', 'varID', 'REF', 'ALT', 'IID', 'GT'])

        variant_file = Path(f'{self._output_prefix}.{tarball_prefix}.{gene_info["SYMBOL"]}.variant_table.tsv')
        carriers_file = Path(f'{self._output_prefix}.{tarball_prefix}.{gene_info["SYMBOL"]}.carriers_formatted.tsv')
        geno_table.to_csv(path_or_buf=variant_file, index=False, sep="\t", na_rep='NA')
        carriers_table.to_csv(path_or_buf=carriers_file, index=False, sep="\t", na_rep='NA')

        return [variant_file, carriers_file]

    def _run_linear_models(self):

        print("Loading data and running null Linear Model")
        null_model = linear_model.linear_model_null(self._association_pack.pheno_names[0],
                                                    self._association_pack.is_binary,
                                                    self._association_pack.found_quantitative_covariates,
                                                    self._association_pack.found_categorical_covariates)

        # 2. Load the tarballs INTO separate genotypes dictionaries
        print("Loading Linear Model genotypes")
        thread_utility = ThreadUtility(self._association_pack.threads,
                                       error_message='A GLM thread failed',
                                       incrementor=10,
                                       thread_factor=2)

        for tarball_prefix in self._association_pack.tarball_prefixes:
            thread_utility.launch_job(linear_model.load_tarball_linear_model,
                                      tarball_prefix=tarball_prefix,
                                      is_snp_tar=self._association_pack.is_snp_tar,
                                      is_gene_tar=self._association_pack.is_gene_tar)

        genotype_packs = {}
        for result in thread_utility:
            tarball_prefix, genotype_dict = result
            genotype_packs[tarball_prefix] = genotype_dict

        # 3. Iterate through every model / gene (in linear_model_pack['genes']) pair and run a GLM
        print("Submitting Linear Models to threads")
        thread_utility = ThreadUtility(self._association_pack.threads,
                                       error_message='A GLM thread failed',
                                       incrementor=500,
                                       thread_factor=1)

        for model in genotype_packs:
            for gene_info in self._gene_infos:
                thread_utility.launch_job(linear_model.run_linear_model,
                                          linear_model_pack=null_model,
                                          genotype_table=genotype_packs[model],
                                          gene=gene_info.name,
                                          mask_name=model,
                                          is_binary=self._association_pack.is_binary,
                                          always_run_corrected=True)

        # As futures finish, write unformatted results:
        fieldnames = ['ENST', 'maskname', 'pheno_name', 'p_val_init', 'n_car', 'cMAC', 'n_model',
                      'p_val_full', 'effect', 'std_err']
        # Binary traits get an additional set of fields to describe the confusion matrix.
        if self._association_pack.is_binary:
            fieldnames.extend(['n_noncar_affected', 'n_noncar_unaffected', 'n_car_affected', 'n_car_unaffected'])

        lm_stats_file = open(self._output_prefix + '.lm_stats.tmp', 'w')
        lm_stats_writer = csv.DictWriter(lm_stats_file,
                                         delimiter="\t",
                                         fieldnames=fieldnames,
                                         extrasaction='ignore')

        lm_stats_writer.writeheader()
        for result in thread_utility:
            finished_gene: LinearModelResult = result
            lm_stats_writer.writerow(finished_gene.todict())
        lm_stats_file.close()

        # 5. Annotate unformatted results and print final outputs
        print("Annotating Linear Model results")
        process_linear_model_outputs(self._output_prefix,
                                     self._association_pack.is_snp_tar, self._association_pack.is_gene_tar,
                                     self._gene_infos)

    def _run_staar_models(self):

        valid_staar_chromosomes = set()

        # Create a file of genes if genes_to_run is !none
        if self._gene_infos is not None:
            with open('staar.gene_list', 'w') as gene_list_file:
                for gene_info in self._gene_infos:
                    gene_list_file.write(gene_info.name + '\n')
                    if gene_info['chrom'] not in valid_staar_chromosomes:
                        valid_staar_chromosomes.add(gene_info['chrom'])
                gene_list_file.close()

        # 1. Run the STAAR NULL model
        print("Running STAAR Null Model(s)...")
        thread_utility = ThreadUtility(self._association_pack.threads,
                                       error_message='A STAAR thread failed',
                                       incrementor=10,
                                       thread_factor=1)
        for phenoname in self._association_pack.pheno_names:
            thread_utility.launch_job(staar_null,
                                      phenoname=phenoname,
                                      is_binary=self._association_pack.is_binary,
                                      found_quantitative_covariates=self._association_pack.found_quantitative_covariates,
                                      found_categorical_covariates=self._association_pack.found_categorical_covariates)
        thread_utility.collect_futures()

        # 2. Run the actual per-gene association tests
        print("Running STAAR masks * chromosomes...")
        thread_utility = ThreadUtility(self._association_pack.threads,
                                       error_message='A STAAR thread failed',
                                       incrementor=10,
                                       thread_factor=1)

        for phenoname in self._association_pack.pheno_names:
            for tarball_prefix in self._association_pack.tarball_prefixes:
                for chromosome in valid_staar_chromosomes:
                    if exists(tarball_prefix + "." + chromosome + ".STAAR.matrix.rds"):
                        thread_utility.launch_job(staar_genes,
                                                  tarball_prefix=tarball_prefix,
                                                  chromosome=chromosome,
                                                  phenoname=phenoname,
                                                  has_gene_info=True)

        # 3. Print a preliminary STAAR output
        print("Finalising STAAR outputs...")
        completed_staar_files = []
        # And gather the resulting futures
        for result in thread_utility:
            tarball_prefix, finished_chromosome, phenoname = result
            completed_staar_files.append(
                f'{tarball_prefix}.{phenoname}.{finished_chromosome}.STAAR_results.tsv')

        # 4. Annotate STAAR output
        process_staar_outputs(completed_staar_files, self._output_prefix,
                              self._association_pack.is_snp_tar, self._association_pack.is_gene_tar, self._gene_infos)
