from extract.extract_association_pack import ExtractAssociationPack, ExtractProgramArgs
from general_utilities.import_utils.import_lib import ingest_wes_bgen, ingest_tarballs, TarballType
from general_utilities.import_utils.module_loader.ingest_data import IngestData


class ExtractIngestData(IngestData):
    """
    Ingest and validate data for the extract module.
    """

    def __init__(self, parsed_options: ExtractProgramArgs):
        super().__init__(parsed_options)

        if parsed_options.gene_ids is None:
            raise ValueError('--gene_ids is required for extraction')

        if parsed_options.association_tarballs is None:
            raise ValueError('--association_tarballs is required for gene extraction')

        # Process tarball information if provided
        tarball_type = None
        tarball_prefixes = None

        if parsed_options.association_tarballs is not None:
            tarball_type, tarball_prefixes = ingest_tarballs(parsed_options.association_tarballs)

        # Load BGEN file information
        bgen_dict = ingest_wes_bgen(parsed_options.bgen_index)

        # Handle optional files
        sparse_grm = (
            parsed_options.sparse_grm.get_file_handle()
            if parsed_options.sparse_grm else None
        )
        sparse_grm_sample = (
            parsed_options.sparse_grm_sample.get_file_handle()
            if parsed_options.sparse_grm_sample else None
        )

        # Create association pack
        self.set_association_pack(ExtractAssociationPack(
            association_pack=self.get_association_pack(),
            tarball_type=tarball_type,
            is_snp_tar=tarball_type == TarballType.SNP if tarball_type else False,
            is_gene_tar=tarball_type == TarballType.GENE if tarball_type else False,
            tarball_prefixes=tarball_prefixes if tarball_prefixes else [],
            bgen_dict=bgen_dict,
            gene_ids=parsed_options.gene_ids,
            sparse_grm=sparse_grm,
            sparse_grm_sample=sparse_grm_sample,
        ))
