import dxpy
from general_utilities.import_utils.file_handlers.input_file_handler import InputFileHandler

from extract.extract_association_pack import ExtractAssociationPack, ExtractProgramArgs
from general_utilities.import_utils.genetics_loader import GeneticsLoader
from general_utilities.import_utils.import_lib import ingest_wes_bgen, ingest_tarballs, TarballType
from general_utilities.import_utils.module_loader.ingest_data import IngestData


class ExtractIngestData(IngestData):

    def __init__(self, parsed_options: ExtractProgramArgs):
        super().__init__(parsed_options)

        # Put additional options/covariate processing required by this specific package here
        tarball_type, tarball_prefixes = ingest_tarballs(parsed_options.association_tarballs)
        bgen_dict = ingest_wes_bgen(parsed_options.bgen_index)

        if tarball_type != TarballType.SNP and tarball_type != TarballType.GENE and parsed_options.gene_ids is None:
            raise dxpy.AppError('Must provide gene IDs when NOT using a SNP/GENE tarball!')

        sparse_grm = parsed_options.sparse_grm.get_file_handle()
        sparse_grm_sample = parsed_options.sparse_grm_sample.get_file_handle()

        # Put additional covariate processing specific to this module here
        self.set_association_pack(ExtractAssociationPack(association_pack=self.get_association_pack(),
                                                         tarball_type=tarball_type,
                                                         is_snp_tar=tarball_type==TarballType.SNP,
                                                         is_gene_tar=tarball_type==TarballType.GENE,
                                                         tarball_prefixes=tarball_prefixes,
                                                         bgen_dict=bgen_dict, gene_ids=parsed_options.gene_ids,
                                                         sparse_grm=sparse_grm,
                                                         sparse_grm_sample=sparse_grm_sample))
