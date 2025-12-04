import csv
import json
from pathlib import Path
from typing import List, Tuple

import dxpy
import pandas as pd
import pandas.core.series
from general_utilities.association_resources import bgzip_and_tabix
from general_utilities.bgen_utilities.genotype_matrix import generate_csr_matrix_from_bgen
from general_utilities.import_utils.file_handlers.export_file_handler import ExportFileHandler
from general_utilities.import_utils.file_handlers.input_file_handler import InputFileHandler
from general_utilities.job_management.command_executor import build_default_command_executor
from general_utilities.job_management.joblauncher_factory import joblauncher_factory
from general_utilities.job_management.thread_utility import ThreadUtility
from general_utilities.linear_model import linear_model
from general_utilities.import_utils.import_lib import LOGGER
from general_utilities.linear_model.linear_model import LinearModelResult
from general_utilities.linear_model.proccess_model_output import merge_glm_staar_runs, process_model_outputs
from general_utilities.linear_model.staar_model import staar_null, staar_genes, load_staar_genetic_data
from general_utilities.mrc_logger import MRCLogger
from scipy.io import mmwrite

from extract.extract_association_pack import ExtractAssociationPack


class GeneExtractionPipeline:
    def __init__(self, output_prefix: str, association_pack: ExtractAssociationPack,
                 gene_infos: list, gene_chunk_map: list, chromosomes: set,
                 transcripts_table: pd.DataFrame):
        self._logger = MRCLogger(__name__).get_logger()
        self._output_prefix = output_prefix
        self._association_pack = association_pack
        self._gene_infos = gene_infos
        self._gene_chunk_map = gene_chunk_map
        self._chromosomes = chromosomes
        self._transcripts_table = transcripts_table
        self._outputs: List[Path] = []

    def run(self) -> List[Path]:
        self._logger.info("Loading VEP annotations...")
        self._download_vep()

        self._logger.info("Filtering variant files to appropriate individuals...")
        self._filter_all_individuals()

        self._logger.info("Extracting variant information...")
        self._collect_variant_information()

        self._logger.info("Running linear models...")
        self._run_linear_models()
        self._run_staar_models()

        merged_outputs = merge_glm_staar_runs(self._output_prefix,
                                              self._association_pack.is_snp_tar,
                                              self._association_pack.is_gene_tar)
        # Flatten in case it returns [(tuple)] instead of [item1, item2]
        for item in merged_outputs:
            if isinstance(item, tuple):
                self._outputs.extend(item)
            else:
                self._outputs.append(item)

        self._outputs.append(self._export_phenotypes())

        return self._outputs

    def _download_vep(self) -> None:
        for chunk_key in self._chromosomes:
            vep_dx = self._association_pack.bgen_dict[chunk_key]['vep']
            output_filename = f"{chunk_key}.filtered.vep.tsv.gz"
            vep_dx.get_file_handle().rename(output_filename)

    def _filter_all_individuals(self) -> None:
        thread_utility = ThreadUtility(self._association_pack.threads, thread_factor=4)

        if self._association_pack.is_snp_tar:
            jobs = [('SNP', prefix) for prefix in self._association_pack.tarball_prefixes]
        elif self._association_pack.is_gene_tar:
            jobs = [('GENE', prefix) for prefix in self._association_pack.tarball_prefixes]
        else:
            jobs = [
                (chromosome, prefix)
                for chromosome in self._chromosomes
                for prefix in self._association_pack.tarball_prefixes
            ]

        for chromosome, tarball_prefix in jobs:
            thread_utility.launch_job(
                function=self._filter_individuals,
                inputs={
                    'tarball_prefix': tarball_prefix,
                    'chromosome': chromosome
                }
            )

        thread_utility.submit_and_monitor()

    def _collect_variant_information(self) -> None:
        thread_utility = ThreadUtility(self._association_pack.threads, thread_factor=2)

        for gene_info, chromosome in self._gene_chunk_map:
            for tarball_prefix in self._association_pack.tarball_prefixes:
                thread_utility.launch_job(
                    function=self._annotate_variants,
                    inputs={
                        'tarball_prefix': tarball_prefix,
                        'gene_info': gene_info,
                        'chromosomes': chromosome
                    },
                    outputs=['variant_file', 'carriers_file']
                )

        thread_utility.submit_and_monitor()
        for result in thread_utility:
            for key in ['variant_file', 'carriers_file']:
                value = result.get(key)
                if value:
                    if isinstance(value, (list, tuple)):
                        self._outputs.extend(value)
                    else:
                        self._outputs.append(value)

    def _filter_individuals(self, tarball_prefix: str, chromosome: str) -> None:
        cmd = (
            f"plink2 "
            f"--bgen {tarball_prefix}.{chromosome}.BOLT.bgen ref-first "
            f"--sample {tarball_prefix}.{chromosome}.BOLT.sample "
            f"--keep SAMPLES_Include.txt "
            f"--make-bed "
            f"--out {tarball_prefix}.{chromosome}.saige_input"
        )
        self._association_pack.cmd_executor.run_cmd_on_docker(cmd)

        cmd_vcf = (
            f"plink2 "
            f"--bfile {tarball_prefix}.{chromosome}.saige_input "
            f"--recode vcf bgz "
            f"--out {tarball_prefix}.{chromosome}.saige_input"
        )
        self._association_pack.cmd_executor.run_cmd_on_docker(cmd_vcf)

        cmd_bcf = (
            f"bcftools view "
            f"-O b "
            f"-o {tarball_prefix}.{chromosome}.saige_input.bcf "
            f"{tarball_prefix}.{chromosome}.saige_input.vcf.gz"
        )
        self._association_pack.cmd_executor.run_cmd_on_docker(cmd_bcf)

    def _annotate_variants(self, tarball_prefix: str, gene_info: pandas.core.series.Series,
                           chromosomes: set) -> List[Path]:
        variant_index = pd.read_csv(f'{chromosomes}.filtered.vep.tsv.gz',
                                    sep="\t",
                                    dtype={'SIFT': str, 'POLYPHEN': str})

        with Path(f'{tarball_prefix}.{chromosomes}.SAIGE.groupFile.txt').open('r') as saige_group_file, \
                Path(f'{tarball_prefix}.{gene_info["SYMBOL"]}.variants.txt').open('w') as var_file:
            var_ids = []
            found_gene = False
            for line in saige_group_file:
                data = line.rstrip().split("\t")
                if data[0] == gene_info.name:
                    found_gene = True

                    for i in range(1, len(data)):
                        raw_id = data[i]

                        # Skip the header column 'var' if present
                        if raw_id == 'var':
                            continue

                        # Store both formats for matching
                        # Format 1: Colon-separated (e.g., chr7:44145173:GC:G) - for output file
                        curr_id_colon = raw_id.replace('_', ':').replace('/', ':')
                        var_file.write(curr_id_colon + "\n")

                        # Format 2: Underscore-separated (e.g., chr7_44145173_GC_G) - for VEP matching
                        curr_id_underscore = raw_id.replace(':', '_').replace('/', '_')

                        # Add both formats to var_ids for matching
                        var_ids.append(curr_id_colon)
                        var_ids.append(curr_id_underscore)

                    break

            if not found_gene:
                self._logger.warning(f"Gene {gene_info.name} NOT found in SAIGE group file!")

        relevant_vars = variant_index[variant_index['varID'].isin(var_ids)]

        variants_txt = f"{tarball_prefix}.{gene_info['SYMBOL']}.variants.txt"
        regions_txt = f"{tarball_prefix}.{gene_info['SYMBOL']}.regions.txt"

        with open(variants_txt) as fin, open(regions_txt, "w") as fout:
            for line in fin:
                stripped = line.strip()
                if not stripped or ":" not in stripped:
                    continue
                parts = stripped.split(":")
                if len(parts) < 2:
                    continue
                chrom, pos = parts[0], parts[1]
                # Strip 'chr' prefix to make it work with both formats
                chrom = chrom.replace('chr', '')
                fout.write(f"{chrom}\t{pos}\t{pos}\n")

        index_cmd = f"bcftools index -f {tarball_prefix}.{chromosomes}.saige_input.bcf"
        self._association_pack.cmd_executor.run_cmd_on_docker(index_cmd)

        cmd = (
            f"bcftools view --threads 2 "
            f"-R {regions_txt} "
            f"-Ob -o {tarball_prefix}.{gene_info['SYMBOL']}.variant_filtered.bcf "
            f"{tarball_prefix}.{chromosomes}.saige_input.bcf"
        )
        self._association_pack.cmd_executor.run_cmd_on_docker(cmd)

        cmd = (
            f"bcftools +fill-tags --threads 4 -Ob "
            f"-o {tarball_prefix}.{gene_info['SYMBOL']}.final.bcf "
            f"{tarball_prefix}.{gene_info['SYMBOL']}.variant_filtered.bcf"
        )
        self._association_pack.cmd_executor.run_cmd_on_docker(cmd)

        cmd = (
            f"bcftools query -f '%ID\\t%MAF\\t%AC\\t%AC_Het\\t%AC_Hom\\n' "
            f"-o {tarball_prefix}.{gene_info['SYMBOL']}.annotated_vars.txt "
            f"{tarball_prefix}.{gene_info['SYMBOL']}.final.bcf"
        )
        self._association_pack.cmd_executor.run_cmd_on_docker(cmd)

        cmd = (
            f"bcftools query -i \"GT='alt'\" -f "
            f"'[%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%SAMPLE\\t%GT\\n]' "
            f"-o {tarball_prefix}.{gene_info['SYMBOL']}.carriers.txt "
            f"{tarball_prefix}.{gene_info['SYMBOL']}.final.bcf"
        )
        self._association_pack.cmd_executor.run_cmd_on_docker(cmd)

        geno_table = pd.read_csv(f"{tarball_prefix}.{gene_info['SYMBOL']}.annotated_vars.txt",
                                 sep="\t",
                                 names=['varID', 'MAF_tested', 'AC_tested', 'AC_tested_Het', 'AC_tested_Hom'])
        geno_table = pd.merge(relevant_vars, geno_table, on='varID', how="left")

        carriers_table = pd.read_csv(f"{tarball_prefix}.{gene_info['SYMBOL']}.carriers.txt",
                                     sep="\t",
                                     names=['CHROM', 'POS', 'varID', 'REF', 'ALT', 'IID', 'GT'])

        variant_file = Path(f'{self._output_prefix}.{tarball_prefix}.{gene_info["SYMBOL"]}.variant_table.tsv')
        carriers_file = Path(f'{self._output_prefix}.{tarball_prefix}.{gene_info["SYMBOL"]}.carriers_formatted.tsv')
        geno_table.to_csv(path_or_buf=variant_file, index=False, sep="\t", na_rep='NA')
        carriers_table.to_csv(path_or_buf=carriers_file, index=False, sep="\t", na_rep='NA')

        return [variant_file, carriers_file]

    def _run_linear_models(self) -> None:
        null_model = linear_model.linear_model_null(
            phenotype=self._association_pack.pheno_names[0],
            phenofile=self._association_pack.final_covariates,
            is_binary=self._association_pack.is_binary,
            ignore_base=self._association_pack.ignore_base_covariates,
            found_quantitative_covariates=self._association_pack.found_quantitative_covariates,
            found_categorical_covariates=self._association_pack.found_categorical_covariates
        )

        thread_utility = ThreadUtility(self._association_pack.threads, thread_factor=2)
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
            genotype_packs[tarball_prefix] = result['genetic_data']

        thread_utility = ThreadUtility(self._association_pack.threads, thread_factor=1)
        for model in genotype_packs:
            for gene_info in self._gene_infos:
                thread_utility.launch_job(
                    function=linear_model.run_linear_model,
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

        fieldnames = ['ENST', 'mask_name', 'pheno_name', 'p_val_init', 'n_car', 'cMAC', 'n_model',
                      'p_val_full', 'effect', 'std_err']
        if self._association_pack.is_binary:
            fieldnames.extend(['n_noncar_affected', 'n_noncar_unaffected', 'n_car_affected', 'n_car_unaffected'])

        lm_stats_path = Path(f'{self._output_prefix}.lm_stats.tmp')
        with lm_stats_path.open('w') as lm_stats_file:
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

        process_model_outputs(input_models=finished_genes,
                              output_path=Path(f'{self._output_prefix}.genes.glm.stats.tsv'),
                              tarball_type=self._association_pack.tarball_type,
                              transcripts_table=self._transcripts_table)

    def _run_staar_null_wrapper(self, phenoname: str, phenofile: Path, is_binary: bool,
                                ignore_base: bool, found_quantitative_covariates: List[str],
                                found_categorical_covariates: List[str], sex: int,
                                sparse_kinship_file: Path, sparse_kinship_samples: Path) -> Tuple[str, Path]:
        """
        Runs the STAAR null model for the specified phenotype. This function serves as a wrapper
        for the STAAR null model computation, incorporating phenotype data, covariates, and
        kinship information to build the model. It facilitates binary or quantitative trait modeling
        based on the input parameters.

        :param phenoname: The name of the phenotype to be processed.
        :param phenofile: Path to the file containing phenotype data.
        :param is_binary: Indicates whether the phenotype is binary or continuous.
        :param ignore_base: Determines if the base covariates should be ignored.
        :param found_quantitative_covariates: List of found quantitative covariates to be included in the model.
        :param found_categorical_covariates: List of found categorical covariates to be included in the model.
        :param sex: Sex information for covariate adjustments.
        :param sparse_kinship_file: Path to the file containing sparse kinship matrix.
        :param sparse_kinship_samples: Path to the file containing sample information for sparse kinship matrix.

        :return: Dictionary containing:
            - 'phenotype': The name of the phenotype.
            - 'null_model': Path to the generated null model file.
        """
        self._logger.info(f"Running STAAR null model for phenotype {phenoname}")
        null_model_path = staar_null(
            phenofile=phenofile,
            phenotype=phenoname,
            is_binary=is_binary,
            ignore_base=ignore_base,
            found_quantitative_covariates=found_quantitative_covariates,
            found_categorical_covariates=found_categorical_covariates,
            sex=sex,
            sparse_kinship_file=sparse_kinship_file,
            sparse_kinship_samples=sparse_kinship_samples
        )
        return phenoname, null_model_path

    def _run_staar_models(self) -> None:
        """
        Run STAAR null models for each phenotype and chromosome combination.

        Creates a merged covariates file for STAAR null model and runs the null model
        for each phenotype and chromosome combination.
        """
        valid_staar_chromosomes = set()
        if self._gene_infos:
            with open('staar.gene_list', 'w') as gene_list_file:
                for gene_info in self._gene_infos:
                    gene_list_file.write(gene_info.name + '\n')
                    if gene_info['chrom'] not in valid_staar_chromosomes:
                        valid_staar_chromosomes.add(gene_info['chrom'])

            # Create merged covariates file matching STAAR approach
        self._logger.info("Creating merged covariates file for STAAR null model...")

        first_chrom = next(iter(self._chromosomes))
        sample_path = self._association_pack.bgen_dict[first_chrom]["sample"].get_file_handle()

        # Read sample file
        sample = pd.read_csv(sample_path, sep=r"\s+", header=0, dtype={'ID_2': str})
        sample = sample.drop(columns=["sex"], errors="ignore")
        sample = sample.iloc[1:].reset_index(drop=True)

        # Read covariates
        covar = pd.read_csv(self._association_pack.final_covariates, sep=' ', header=0, dtype={'IID': str})

        # Merge
        merged = sample.merge(covar, how="left", left_on="ID_2", right_on="IID")
        merged = merged.drop(columns=["ID_1", "missing"], errors="ignore")
        merged = merged.sort_values("ID_2")

        # Remove rows with missing phenotype
        for phenoname in self._association_pack.pheno_names:
            if phenoname in merged.columns:
                merged = merged.dropna(subset=[phenoname])

        # Remove rows with missing covariates
        required_cols = ['age', 'age_squared', 'batch']
        if self._association_pack.sex == 2:
            required_cols.append('sex')
        for i in range(1, 11):
            required_cols.append(f'PC{i}')
        required_cols.extend(self._association_pack.found_quantitative_covariates)
        required_cols.extend(self._association_pack.found_categorical_covariates)

        required_cols = [col for col in required_cols if col in merged.columns]
        merged = merged.dropna(subset=required_cols)

        self._logger.info(f"After removing samples with missing data: {len(merged)} samples remain")

        # Rename and save
        # Drop IID since we're using FID
        merged = merged.drop(columns=['IID', 'FID'], errors='ignore')
        merged = merged.rename(columns={'ID_2': 'FID'})
        merged_cov_path = Path("merged_covariates_for_staar.tsv")
        merged.to_csv(merged_cov_path, sep="\t", index=False)

        thread_utility = ThreadUtility(self._association_pack.threads, thread_factor=1)
        for phenoname in self._association_pack.pheno_names:
            thread_utility.launch_job(
                function=self._run_staar_null_wrapper,
                inputs={
                    'phenoname': phenoname,
                    'phenofile': merged_cov_path,  # Use the new merged file
                    'is_binary': self._association_pack.is_binary,
                    'ignore_base': self._association_pack.ignore_base_covariates,
                    'found_quantitative_covariates': self._association_pack.found_quantitative_covariates,
                    'found_categorical_covariates': self._association_pack.found_categorical_covariates,
                    'sex': self._association_pack.sex,
                    'sparse_kinship_file': self._association_pack.sparse_grm,
                    'sparse_kinship_samples': self._association_pack.sparse_grm_sample
                },
                outputs=['phenotype', 'null_model']
            )
        thread_utility.submit_and_monitor()

        # Collect null models
        null_models = {}
        for result in thread_utility:
            phenotype = result['phenotype']
            null_model = result['null_model']

            # If null_model is a tuple, extract the Path (second element)
            if isinstance(null_model, tuple):
                null_model = null_model[1]  # Get PosixPath from ('T2D', PosixPath(...))

            null_models[phenotype] = null_model

        self._logger.info(f"Collected {len(null_models)} null models: {list(null_models.keys())}")

        # Filter STAAR samples tables to match the samples in the null model
        # Read the samples that were used in the null model
        final_covars = pd.read_csv(merged_cov_path, sep='\t')
        null_model_samples = set(final_covars['FID'].astype(str))

        # Filter each STAAR samples table
        for tarball_prefix in self._association_pack.tarball_prefixes:
            for chromosome in self._chromosomes:
                samples_path = Path(f"{tarball_prefix}.{chromosome}.STAAR.samples_table.tsv")

                if samples_path.exists():
                    # Read the original STAAR samples table
                    staar_samples_df = pd.read_csv(samples_path, sep='\t')
                    original_count = len(staar_samples_df)

                    # Filter to only keep samples in the null model
                    staar_samples_df = staar_samples_df[
                        staar_samples_df['sampID'].astype(str).isin(null_model_samples)
                    ]
                    filtered_count = len(staar_samples_df)

                    # Overwrite with filtered version
                    staar_samples_df.to_csv(samples_path, sep='\t', index=False)

                    self._logger.info(
                        f"Filtered {samples_path.name}: {original_count} → {filtered_count} samples"
                    )

        launcher = joblauncher_factory(download_on_complete=True)

        for phenoname in self._association_pack.pheno_names:
            for tarball_prefix in self._association_pack.tarball_prefixes:
                for chromosome in self._chromosomes:
                    staar_data = load_staar_genetic_data(
                        tarball_prefix=tarball_prefix,
                        bgen_prefix=chromosome
                    )

                    valid_gene_ids = {gene.name for gene in self._gene_infos}
                    genes_per_chunk = {
                        chunk: [gene for gene in genes.keys() if gene in valid_gene_ids]
                        for chunk, genes in staar_data.items()
                    }

                    for chunk, gene_list in genes_per_chunk.items():
                        if not gene_list:
                            continue

                        subset_staar_data = {chunk: staar_data[chunk]}
                        chunk_json_path = Path(f"{tarball_prefix}.{chunk}.staar_chunk.json")
                        with chunk_json_path.open("w") as f:
                            json.dump(subset_staar_data, f, default=lambda o: list(o) if isinstance(o, set) else o)

                        working_chunk = self._association_pack.bgen_dict[chromosome]

                        exporter = ExportFileHandler(delete_on_upload=False)
                        # null model
                        null_model = exporter.export_files(null_models[phenoname])
                        # staar samples
                        staar_samples = exporter.export_files(f'{tarball_prefix}.{chromosome}.STAAR.samples_table.tsv')
                        # variant table
                        variants_table = exporter.export_files(
                            f'{tarball_prefix}.{chromosome}.STAAR.variants_table.tsv')
                        # chunk file
                        chunk_file = exporter.export_files(chunk_json_path)
                        # transcripts table
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
                                'bgen': working_chunk['bgen'].get_input_str(),
                                'index': working_chunk['index'].get_input_str(),
                                'sample': working_chunk['sample'].get_input_str(),
                                'staar_samples': staar_samples,
                                'staar_variants': variants_table,
                                'tarball_type': str(self._association_pack.tarball_type),
                                'transcripts_table': transcripts_table
                            },
                            outputs=['output_model']
                        )

        launcher.submit_and_monitor()

        completed_staar_chunks = []
        for result in launcher:
            result_file = InputFileHandler(result['output_model']).get_file_handle()
            df = pd.read_csv(result_file, sep='\t', index_col=0)
            completed_staar_chunks.append(df)

        if completed_staar_chunks:
            combined_staar = pd.concat(completed_staar_chunks, axis=0)
            # Merge with transcripts_table to get genomic coordinates
            combined_staar = combined_staar.merge(
                self._transcripts_table[['chrom', 'start', 'end']],
                left_index=True,
                right_index=True,
                how='left'
            )
            combined_staar = combined_staar.sort_values(by='start')
            combined_staar.to_csv(f'{self._output_prefix}.genes.STAAR.stats.tsv', sep='\t', index=True)
            output_tsv = Path(f"{self._output_prefix}.genes.STAAR.stats.tsv")
            outputs = bgzip_and_tabix(output_tsv, skip_row=1, sequence_row=14, begin_row=15, end_row=16)
            self._outputs.extend(outputs)

    def _export_phenotypes(self) -> Path:
        out_pheno_path = Path(f'{self._output_prefix}.phenotypes_covariates.formatted.tsv')
        Path('phenotypes_covariates.formatted.txt').rename(out_pheno_path)
        return out_pheno_path


@dxpy.entry_point('multithread_gene_model')
def multithread_gene_model(null_model, pheno_name, tarball_prefix, chromosome, genes, chunk_file, bgen: str, index: str,
                           sample: str,
                           staar_samples, staar_variants, tarball_type, transcripts_table) -> dict:
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
    bgen_path = InputFileHandler(bgen, download_now=True).get_file_handle()
    index_path = InputFileHandler(index, download_now=True).get_file_handle()
    sample_path = InputFileHandler(sample, download_now=True).get_file_handle()

    thread_utility = ThreadUtility()

    for gene in genes:
        # Check if gene exists in staar_data
        if gene not in staar_data[chromosome]:
            LOGGER.warning(f"Gene {gene} not found in staar_data for {chromosome}")
            continue

        gene_info = staar_data[chromosome][gene]
        n_variants = len(gene_info.get('vars', []))

        # Skip genes with fewer than 2 variants
        if n_variants < 2:
            LOGGER.warning(f"Skipping {gene}: fewer than 2 variants ({n_variants})")
            continue

        LOGGER.info(f"Processing gene {gene} with {n_variants} variants")

        # Read the filtered samples table
        staar_samples_df = pd.read_csv(staar_samples, sep='\t')
        keep_rows = staar_samples_df['row'].values

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
        # Subset matrix to filtered samples
        matrix = matrix[keep_rows, :]
        # Then save and use
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

    exporter = ExportFileHandler()
    uploaded_file = exporter.export_files(output_model)

    return {"output_model": uploaded_file}
