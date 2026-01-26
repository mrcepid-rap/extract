from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

from general_utilities.import_utils.file_handlers.input_file_handler import InputFileHandler
from general_utilities.import_utils.import_lib import BGENInformation, TarballType
from general_utilities.import_utils.module_loader.association_pack import AssociationPack, ProgramArgs


@dataclass
class ExtractProgramArgs(ProgramArgs):
    association_tarballs: InputFileHandler
    bgen_index: InputFileHandler
    gene_ids: List[str]
    sparse_grm: InputFileHandler
    sparse_grm_sample: InputFileHandler

    def __post_init__(self):
        """@dataclass automatically calls this method after calling its own __init__().

        This is required in the subclass because dataclasses do not call the __init__ of their super o.0

        """
        self._check_opts()

    def _check_opts(self):
        pass


class ExtractAssociationPack(AssociationPack):

    def __init__(self, association_pack: AssociationPack, tarball_type: TarballType,
                 tarball_prefixes: List[Path],
                 bgen_dict: Dict[str, BGENInformation], gene_ids: List[str], sparse_grm: Path,
                 sparse_grm_sample: Path):
        super().__init__(association_pack.is_binary, association_pack.sex, association_pack.threads,
                         association_pack.pheno_names,
                         association_pack.found_quantitative_covariates, association_pack.found_categorical_covariates,
                         association_pack.cmd_executor, association_pack.final_covariates,
                         association_pack.inclusion_samples,
                         association_pack.exclusion_samples, association_pack.transcript_index)

        self.tarball_type = tarball_type
        self.tarball_prefixes = tarball_prefixes
        self.bgen_dict = bgen_dict
        self.gene_ids = gene_ids
        self.sparse_grm = sparse_grm
        self.sparse_grm_sample = sparse_grm_sample
