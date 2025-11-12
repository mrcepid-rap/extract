import csv
import json
from pathlib import Path
from typing import List

import dxpy
import pandas as pd
import pandas.core.series
from general_utilities.association_resources import process_snp_or_gene_tar, build_transcript_table, get_gene_id, \
    process_gene_or_snp_wgs, bgzip_and_tabix
from general_utilities.bgen_utilities.genotype_matrix import generate_csr_matrix_from_bgen
from general_utilities.import_utils.file_handlers.export_file_handler import ExportFileHandler
from general_utilities.import_utils.file_handlers.input_file_handler import InputFileHandler
from general_utilities.job_management.command_executor import build_default_command_executor
from general_utilities.job_management.joblauncher_factory import joblauncher_factory
from general_utilities.job_management.thread_utility import ThreadUtility
from general_utilities.linear_model import linear_model
from general_utilities.linear_model.linear_model import LinearModelResult
from general_utilities.linear_model.proccess_model_output import merge_glm_staar_runs, process_model_outputs
from general_utilities.linear_model.staar_model import staar_null, staar_genes, load_staar_genetic_data
from general_utilities.mrc_logger import MRCLogger
from scipy.io import mmwrite

from extract.extract_association_pack import ExtractAssociationPack


class ExtractVariants:

    def __init__(self, output_prefix: str, association_pack: ExtractAssociationPack):

        self._logger = MRCLogger(__name__).get_logger()
        self._outputs = []
        self._output_prefix = output_prefix
        self._association_pack = association_pack

        # Define our gene-list and make 'gene_info' objects of them (which are Pandas series classes)
        self._gene_infos = []
        self._gene_chunk_map = []
        self._chromosomes = set()

        # build the transcripts table
        self._transcripts_table = build_transcript_table(transcripts_path=self._association_pack.transcript_index)

        # If we are doing extraction based on individual SNPs or a Gene list, we need to make a 'fake' gene info
        # but find all chromosomes those SNPS/Genes lie on
        if self._association_pack.is_non_standard_tar:
            gene_info, returned_chromosomes = process_snp_or_gene_tar(self._association_pack.is_snp_tar,
                                                                      self._association_pack.is_gene_tar,
                                                                      self._association_pack.tarball_prefixes[0])
            self._gene_infos.append(gene_info)
            self._chromosomes = returned_chromosomes
            for chromosome in self._chromosomes:
                self._gene_chunk_map.append((gene_info, chromosome))
        else:
            for gene in self._association_pack.gene_ids:
                gene_info = get_gene_id(gene, self._transcripts_table)
                self._gene_infos.append(gene_info)

                # when working with WGS data, we need to find which chunk the gene lies in (too much data to load all)
                for chunk in self._association_pack.bgen_dict:
                    chromosomes = process_gene_or_snp_wgs(
                        identifier=gene_info.name,  # ENST ID
                        tarball_prefix=self._association_pack.tarball_prefixes[0],
                        chunk=chunk
                    )

                    if chromosomes:
                        self._logger.info(f"{gene_info['SYMBOL']} found in {chunk} ({', '.join(chromosomes)})")
                        # Add all detected chromosome chunk for this gene to self._chromosomes
                        self._chromosomes.add(chunk)
                        self._gene_chunk_map.append((gene_info, chunk))


    def run_tool(self):

        # 1. Download variant VEP annotations
        self._logger.info("Loading VEP annotations...")
        thread_utility = ThreadUtility(self._association_pack.threads,
                                       thread_factor=4)

        thread_utility.launch_job(
            function=self._download_vep,
            inputs={}
        )
        thread_utility.submit_and_monitor()

        # 2. Filter relevant files to individuals we want to keep
        self._logger.info("Filtering variant files to appropriate individuals...")
        thread_utility = ThreadUtility(self._association_pack.threads,
                                       thread_factor=4)

        # if, elif, else simply depends on which type of tarball we are using
        if self._association_pack.is_snp_tar:
            for tarball_prefix in self._association_pack.tarball_prefixes:
                thread_utility.launch_job(
                    function=self._filter_individuals,
                    inputs={
                        'tarball_prefix': tarball_prefix,
                        'chromosome': 'SNP'
                    }
                )
        elif self._association_pack.is_gene_tar:
            for tarball_prefix in self._association_pack.tarball_prefixes:
                thread_utility.launch_job(
                    function=self._filter_individuals,
                    inputs={
                        'tarball_prefix': tarball_prefix,
                        'chromosome': 'GENE'
                    }
                )
        else:
            for chromosome in self._chromosomes:
                for tarball_prefix in self._association_pack.tarball_prefixes:
                    thread_utility.launch_job(
                        function=self._filter_individuals,
                        inputs={
                            'tarball_prefix': tarball_prefix,
                            'chromosome': chromosome
                        }
                    )
        thread_utility.submit_and_monitor()

        # 4. Actually collect variant information per-gene
        self._logger.info("Extracting variant information...")
        thread_utility = ThreadUtility(self._association_pack.threads,
                                       thread_factor=2)
        # Only annotate genes for the chunk(s) they were discovered in
        for gene_info, chromosome in self._gene_chunk_map:
            for tarball_prefix in self._association_pack.tarball_prefixes:
                thread_utility.launch_job(
                    function=self._annotate_variants,
                    inputs={
                        'tarball_prefix': tarball_prefix,
                        'gene_info': gene_info,
                        'chromosomes': chromosome
                    },
                    outputs=[
                        'variant_file', 'carriers_file'
                    ]
                )
        thread_utility.submit_and_monitor()

        for result in thread_utility:
            self._outputs.extend(result)

        # 5. And run a linear and STAAR model(s) for all genes
        self._logger.info("Running linear models...")
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

    def _download_vep(self) -> None:
        """Download the VEP annotations for the relevant chromosome."""

        for chunk_key in self._chromosomes:
            vep_dx = self._association_pack.bgen_dict[chunk_key]['vep']
            output_filename = f"{chunk_key}.filtered.vep.tsv.gz"
            vep_file = vep_dx.get_file_handle()
            # change name to output_filename
            vep_file.rename(output_filename)

    def _filter_individuals(self, tarball_prefix: str, chromosome: str) -> None:
        """Filter the relevant SAIGE file to just the individuals we want so we can get actual MAC.

        :param tarball_prefix: The prefix of the tarball to filter
        :param chromosome: The chromosome to filter on, or 'SNP'/'GENE' if this is a SNP/GENE tarball
        :return: None
        """
        # And filter the relevant SAIGE file to just the individuals we want so we can get actual MAC
        cmd = (
            f"plink2 "
            f"--bgen {tarball_prefix}.{chromosome}.BOLT.bgen ref-first "
            f"--sample {tarball_prefix}.{chromosome}.BOLT.sample "
            f"--keep SAMPLES_Include.txt "
            f"--make-bed "
            f"--out {tarball_prefix}.{chromosome}.saige_input"
        )
        self._association_pack.cmd_executor.run_cmd_on_docker(cmd)

        # Convert PLINK binary to compressed VCF (bgzipped)
        cmd_vcf = (
            f"plink2 "
            f"--bfile {tarball_prefix}.{chromosome}.saige_input "
            f"--recode vcf bgz "
            f"--out {tarball_prefix}.{chromosome}.saige_input"
        )
        self._association_pack.cmd_executor.run_cmd_on_docker(cmd_vcf)

        # Convert compressed VCF to BCF
        cmd_bcf = (
            f"bcftools view "
            f"-O b "
            f"-o {tarball_prefix}.{chromosome}.saige_input.bcf "
            f"{tarball_prefix}.{chromosome}.saige_input.vcf.gz"
        )
        self._association_pack.cmd_executor.run_cmd_on_docker(cmd_bcf)

    def _annotate_variants(self, tarball_prefix: str, gene_info: pandas.core.series.Series,
                           chromosomes: set) -> List[Path]:

        # This is a bit confusing, so explaining in full.
        # We need to annotate EACH GENE separately EXCEPT when running a SNP/GENE list tarball, SO...
        # 1. If we have a SNP/GENE tar, we may need to load in variant annotations for multiple chromosomes,
        # so we set 'chromosomes'
        # to a set of chromosomes that we extract from the SNP/GENE tar.
        # 2. If just a single gene or gene list (chromosomes = None), only need to load the data for the chromosome
        # that specific Gene is on

        variant_index = []
        variant_index.append(pd.read_csv(f'{chromosomes}.filtered.vep.tsv.gz',
                                         sep="\t",
                                         dtype={'SIFT': str, 'POLYPHEN': str}))
        variant_index = pd.concat(variant_index)

        # Need to get the variants from the SAIGE groupfile:
        with Path(f'{tarball_prefix}.{chromosomes}.SAIGE.groupFile.txt').open('r') as saige_group_file, \
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
        # cmd = f'bcftools view --threads 2 -i \'ID=@{tarball_prefix}.{gene_info["SYMBOL"]}.variants.txt\' -Ob ' \
        #       f'-o {tarball_prefix}.{gene_info["SYMBOL"]}.variant_filtered.bcf ' \
        #       f'{tarball_prefix}.{chromosomes}.saige_input.bcf'
        # self._association_pack.cmd_executor.run_cmd_on_docker(cmd)

        # Create a regions file (CHROM, POS, END)
        variants_txt = f"{tarball_prefix}.{gene_info['SYMBOL']}.variants.txt"
        regions_txt = f"{tarball_prefix}.{gene_info['SYMBOL']}.regions.txt"

        with open(variants_txt) as fin, open(regions_txt, "w") as fout:
            for line in fin:
                stripped = line.strip()
                if not stripped:
                    continue
                if ":" not in stripped:
                    continue
                parts = stripped.split(":")
                if len(parts) < 2:
                    continue
                chrom, pos = parts[0], parts[1]
                fout.write(f"{chrom}\t{pos}\t{pos}\n")

        # Ensure BCF is indexed
        index_cmd = f"bcftools index -f {tarball_prefix}.{chromosomes}.saige_input.bcf"
        self._association_pack.cmd_executor.run_cmd_on_docker(index_cmd)

        # Filter to the variants for this gene
        cmd = (
            f"bcftools view --threads 2 "
            f"-R {regions_txt} "
            f"-Ob -o {tarball_prefix}.{gene_info['SYMBOL']}.variant_filtered.bcf "
            f"{tarball_prefix}.{chromosomes}.saige_input.bcf"
        )
        self._association_pack.cmd_executor.run_cmd_on_docker(cmd)

        cmd = f'bcftools +fill-tags --threads 4 -Ob ' \
              f'-o {tarball_prefix}.{gene_info["SYMBOL"]}.final.bcf ' \
              f'{tarball_prefix}.{gene_info["SYMBOL"]}.variant_filtered.bcf'
        self._association_pack.cmd_executor.run_cmd_on_docker(cmd)

        # Now get actual annotations back in:
        cmd = f'bcftools query -f \'%ID\\t%MAF\\t%AC\\t%AC_Het\\t%AC_Hom\\n\' ' \
              f'-o {tarball_prefix}.{gene_info["SYMBOL"]}.annotated_vars.txt ' \
              f'{tarball_prefix}.{gene_info["SYMBOL"]}.final.bcf'
        self._association_pack.cmd_executor.run_cmd_on_docker(cmd)

        # And get a list of individuals with a variant:
        cmd = f'bcftools query -i \"GT=\'alt\'\" -f \'[%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%SAMPLE\\t%GT\n]\' ' \
              f'-o {tarball_prefix}.{gene_info["SYMBOL"]}.carriers.txt ' \
              f'{tarball_prefix}.{gene_info["SYMBOL"]}.final.bcf'
        self._association_pack.cmd_executor.run_cmd_on_docker(cmd)

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

        self._logger.info("Loading data and running null Linear Model")
        null_model = linear_model.linear_model_null(phenotype=self._association_pack.pheno_names[0],
                                                    phenofile=self._association_pack.final_covariates,
                                                    is_binary=self._association_pack.is_binary,
                                                    ignore_base=self._association_pack.ignore_base_covariates,
                                                    found_quantitative_covariates=self._association_pack.found_quantitative_covariates,
                                                    found_categorical_covariates=self._association_pack.found_categorical_covariates)

        # 2. Load the tarballs INTO separate genotypes dictionaries
        self._logger.info("Loading Linear Model genotypes")
        thread_utility = ThreadUtility(self._association_pack.threads,
                                       thread_factor=2)

        for chromosome in self._chromosomes:
            for tarball_prefix in self._association_pack.tarball_prefixes:
                thread_utility.launch_job(
                    function=linear_model.load_linear_model_genetic_data,
                    inputs={
                        'tarball_prefix': tarball_prefix,
                        'tarball_type': self._association_pack.tarball_type,
                        'bgen_prefix': chromosome,
                    },
                    outputs=['tarball_prefix', 'genetic_data']
                )
        thread_utility.submit_and_monitor()

        genotype_packs = {}
        for result in thread_utility:
            tarball_prefix = result['tarball_prefix']
            genotype_dict = result['genetic_data']
            genotype_packs[tarball_prefix] = genotype_dict

        # 3. Iterate through every model / gene (in linear_model_pack['genes']) pair and run a GLM
        self._logger.info("Submitting Linear Models to threads")
        thread_utility = ThreadUtility(self._association_pack.threads,
                                       thread_factor=1)

        for model in genotype_packs:
            for gene_info in self._gene_infos:
                thread_utility.launch_job(function=linear_model.run_linear_model,
                                          inputs={
                                              'linear_model_pack': null_model,
                                              'genotype_table': genotype_packs[model],
                                              'gene': gene_info.name,
                                              'mask_name': model,
                                              'is_binary': self._association_pack.is_binary,
                                              'always_run_corrected': True
                                          },
                                          outputs=['gene_dict']
                                          )
        thread_utility.submit_and_monitor()

        # As futures finish, write unformatted results:
        fieldnames = ['ENST', 'mask_name', 'pheno_name', 'p_val_init', 'n_car', 'cMAC', 'n_model',
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
        finished_genes = []
        for result in thread_utility:
            finished_gene: LinearModelResult = result['gene_dict']
            lm_stats_writer.writerow(finished_gene.todict())
            finished_genes.append(finished_gene)
        lm_stats_file.close()

        # 5. Annotate unformatted results and print final outputs
        self._logger.info("Annotating Linear Model results")
        process_model_outputs(input_models=finished_genes,
                              output_path=Path(f'{self._output_prefix}.lm_results.tsv'),
                              tarball_type=self._association_pack.tarball_type,
                              transcripts_table=self._transcripts_table)

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
        self._logger.info("Running STAAR Null Model(s)...")
        thread_utility = ThreadUtility(self._association_pack.threads,
                                       thread_factor=1)
        for phenoname in self._association_pack.pheno_names:
            thread_utility.launch_job(function=staar_null,
                                      inputs={
                                          'phenofile': self._association_pack.final_covariates,
                                          'phenotype': phenoname,
                                          'is_binary': self._association_pack.is_binary,
                                          'ignore_base': self._association_pack.ignore_base_covariates,
                                          'found_quantitative_covariates': self._association_pack.found_quantitative_covariates,
                                          'found_categorical_covariates': self._association_pack.found_categorical_covariates,
                                          'sex': self._association_pack.sex,
                                          'sparse_kinship_file': self._association_pack.sparse_grm,
                                          'sparse_kinship_samples': self._association_pack.sparse_grm_sample
                                      },
                                      outputs=['staar_null_model']
                                      )
        thread_utility.submit_and_monitor()

        # 2. Run the actual per-gene association tests
        self._logger.info("Running STAAR masks * chromosomes...")

        # set the job launcher
        launcher = joblauncher_factory(download_on_complete=True)

        for phenoname in self._association_pack.pheno_names:
            for tarball_prefix in self._association_pack.tarball_prefixes:
                for chromosome in self._chromosomes:

                    staar_data = load_staar_genetic_data(
                        tarball_prefix=tarball_prefix,
                        bgen_prefix=chromosome
                    )

                    # Build mapping
                    valid_gene_ids = {gene.name for gene in self._gene_infos}
                    genes_per_chunk = {
                        chunk: [gene for gene in genes.keys() if gene in valid_gene_ids]
                        for chunk, genes in staar_data.items()
                    }

                    # we only want to run for chunks that have genes to run
                    for chunk, gene_list in genes_per_chunk.items():
                        if not gene_list:
                            continue

                        # send the genetic information to each subjob
                        subset_staar_data = {chunk: staar_data[chunk]}
                        chunk_json_path = Path(f"{tarball_prefix}.{chunk}.staar_chunk.json")
                        with chunk_json_path.open("w") as f:
                            json.dump(subset_staar_data, f, default=lambda o: list(o) if isinstance(o, set) else o)

                        # set the chunk that we are working with
                        working_chunk = self._association_pack.bgen_dict[chromosome]

                        # export the files we will need
                        exporter = ExportFileHandler(delete_on_upload=False)
                        null_model = exporter.export_files(f'{phenoname}.STAAR_null.rds')
                        staar_samples = exporter.export_files(f'{tarball_prefix}.{chromosome}.STAAR.samples_table.tsv')
                        variants_table = exporter.export_files(
                            f'{tarball_prefix}.{chromosome}.STAAR.variants_table.tsv')
                        chunk_file = exporter.export_files(chunk_json_path)
                        transcripts_table = Path("transcripts_table.tsv")
                        self._transcripts_table.to_csv(transcripts_table, sep='\t', index=True)
                        transcripts_table = exporter.export_files(transcripts_table)

                        launcher.launch_job(
                            function=multithread_gene_model,
                            inputs={
                                'null_model': null_model,
                                'pheno_name': phenoname,
                                'tarball_prefix': tarball_prefix,
                                'chromosome': chromosome,
                                'genes': gene_list,
                                'chunk_file': chunk_file,
                                'bgen': working_chunk['bgen'],
                                'index': working_chunk['index'],
                                'sample': working_chunk['sample'],
                                'staar_samples': staar_samples,
                                'staar_variants': variants_table,
                                'tarball_type': self._association_pack.tarball_type,
                                'transcripts_table': transcripts_table
                            },
                            outputs=['output_model']
                        )
        launcher.submit_and_monitor()

        # Gather all the results
        completed_staar_chunks = []
        for result in launcher:
            df = pd.read_csv(result['output_model'], sep='\t', index_col=0)
            completed_staar_chunks.append(df)
        # Combine them all
        combined_staar = pd.concat(completed_staar_chunks, axis=0)
        combined_staar = combined_staar.sort_values(by='start')
        combined_staar.to_csv(f'{self._output_prefix}.staar_results.tsv', sep='\t', index=True)
        output_tsv = Path(f"{self._output_prefix}.staar_results.tsv")
        outputs = bgzip_and_tabix(output_tsv, skip_row=1, sequence_row=2, begin_row=3,
                                  end_row=4)

        return outputs


@dxpy.entry_point('multithread_gene_model')
def multithread_gene_model(null_model, pheno_name, tarball_prefix, chromosome, genes, chunk_file, bgen, index, sample,
                           staar_samples, staar_variants, tarball_type, transcripts_table) -> Path:
    """
    Run a STAAR gene model in a multithreaded way

    :param null_model: a path to the null model RDS file
    :param pheno_name: the phenotype name
    :param tarball_prefix: the tarball prefix to work with
    :param chromosome: the chromosome chunk to work with
    :param genes: list of genes to run
    :param chunk_file: a path to the chunk JSON file (contains the genetic coordinates & variants)
    :param bgen: InputFileHandler for the bgen file
    :param index: InputFileHandler for the bgen index file
    :param sample: InputFileHandler for the bgen sample file
    :param staar_samples: STAAR samples table for the chunk we are working with
    :param staar_variants: STAAR variants table for the chunk we are working with
    :param tarball_type: the tarball type (TarballType enum)
    :param transcripts_table: a path to the transcripts table
    :return: Path to the output STAAR results TSV file (post-annotation)
    """

    # load our VM environment
    cmd_executor = build_default_command_executor()
    null_model = InputFileHandler(null_model).get_file_handle()
    staar_samples = InputFileHandler(staar_samples).get_file_handle()
    staar_variants = InputFileHandler(staar_variants).get_file_handle()
    chunk_file = InputFileHandler(chunk_file).get_file_handle()
    transcripts_table = InputFileHandler(transcripts_table).get_file_handle()

    with open(chunk_file, "r") as f:
        staar_data = json.load(f)

    # download our bgen files
    bgen_path = bgen.get_file_handle()
    index_path = index.get_file_handle()
    sample_path = sample.get_file_handle()

    thread_utility = ThreadUtility()

    for gene in genes:
        # generate a csr matrix from the bgen files
        matrix, summary_dict = generate_csr_matrix_from_bgen(
            bgen_path=bgen_path,
            sample_path=sample_path,
            variant_filter_list=staar_data[chromosome][gene]['vars'],
            chromosome=staar_data[chromosome][gene]['chrom'],
            start=staar_data[chromosome][gene]['min'],
            end=staar_data[chromosome][gene]['max'],
            should_collapse_matrix=False
        )

        # export matrix to file
        mmwrite(f"{tarball_prefix}.{chromosome}.STAAR.mtx", matrix)

        thread_utility.launch_job(
            function=staar_genes,
            inputs={
                'staar_null_path': null_model,
                'pheno_name': pheno_name,
                'gene': gene,  # single ENST ID string
                'mask_name': tarball_prefix,
                'staar_matrix': f"{tarball_prefix}.{chromosome}.STAAR.mtx",
                'staar_samples': staar_samples,
                'staar_variants': staar_variants,
                'out_dir': Path('.'),
            },
            outputs=['staar_result']
        )
    thread_utility.submit_and_monitor()
    # Print a preliminary STAAR output
    completed_staar_files = []
    # And gather the resulting futures
    for result in thread_utility:
        # Each result is a dict with {'staar_result': STAARModelResult(...)}
        staar_result = result["staar_result"]
        completed_staar_files.append(staar_result)

    # Annotate STAAR output
    transcript = pd.read_csv(transcripts_table, sep='\t', index_col=0)
    output_model = Path(f'{chromosome}.staar_results.tsv')
    process_model_outputs(input_models=completed_staar_files,
                          output_path=output_model,
                          tarball_type=tarball_type,
                          transcripts_table=transcript)

    return output_model
