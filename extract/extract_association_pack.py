from dataclasses import dataclass
from typing import Dict, List

import dxpy

from general_utilities.import_utils.import_lib import BGENInformation
from general_utilities.import_utils.module_loader.association_pack import AssociationPack, ProgramArgs


@dataclass
class ExtractProgramArgs(ProgramArgs):
    association_tarballs: dxpy.DXFile
    bgen_index: dxpy.DXFile
    gene_ids: List[str]
    sparse_grm: dxpy.DXFile
    sparse_grm_sample: dxpy.DXFile

    def __post_init__(self):
        """@dataclass automatically calls this method after calling its own __init__().

        This is required in the subclass because dataclasses do not call the __init__ of their super o.0

        """
        self._check_opts()

    def _check_opts(self):
        pass


class ExtractAssociationPack(AssociationPack):

    def __init__(self, association_pack: AssociationPack,
                 is_snp_tar: bool, is_gene_tar: bool, tarball_prefixes: List[str],
                 bgen_dict: Dict[str, BGENInformation], gene_ids: List[str]):

        super().__init__(association_pack.is_binary, association_pack.sex, association_pack.threads,
                         association_pack.pheno_names, association_pack.ignore_base_covariates,
                         association_pack.found_quantitative_covariates, association_pack.found_categorical_covariates,
                         association_pack.cmd_executor)

        self.is_snp_tar = is_snp_tar
        self.is_gene_tar = is_gene_tar
        self.is_non_standard_tar = is_snp_tar or is_gene_tar
        self.tarball_prefixes = tarball_prefixes
        self.bgen_dict = bgen_dict
        self.gene_ids = gene_ids
