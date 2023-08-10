import dxpy

from extract.extract_association_pack import ExtractAssociationPack, ExtractProgramArgs
from general_utilities.import_utils.genetics_loader import GeneticsLoader
from general_utilities.import_utils.import_lib import ingest_wes_bgen, ingest_tarballs
from general_utilities.import_utils.module_loader.ingest_data import IngestData


class ExtractIngestData(IngestData):

    def __init__(self, parsed_options: ExtractProgramArgs):
        super().__init__(parsed_options)

        # Put additional options/covariate processing required by this specific package here
        is_snp_tar, is_gene_tar, named_prefix, tarball_prefixes = ingest_tarballs(parsed_options.association_tarballs)
        bgen_dict = ingest_wes_bgen(parsed_options.bgen_index)

        if is_snp_tar is False and is_gene_tar is False and parsed_options.gene_ids is None:
            raise dxpy.AppError('Must provide gene IDs when NOT using a SNP/GENE tarball!')

        GeneticsLoader.ingest_sparse_matrix(parsed_options.sparse_grm,
                                            parsed_options.sparse_grm_sample)

        # Put additional covariate processing specific to this module here
        self.set_association_pack(ExtractAssociationPack(self.get_association_pack(),
                                                         is_snp_tar, is_gene_tar, tarball_prefixes,
                                                         bgen_dict, parsed_options.gene_ids))
