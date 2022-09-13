from typing import Type

import dxpy
from extract import extract_ingester
from extract.extract_association_pack import ExtractProgramArgs, ExtractAssociationPack
from extract.extract_variants import ExtractVariants
from runassociationtesting.module_loader import ModuleLoader


class LoadModule(ModuleLoader):

    def __init__(self, output_prefix: str, input_args: str):

        super().__init__(output_prefix, input_args)

    def start_module(self) -> None:

        # Start the extract variants tool
        extract_tool = ExtractVariants(self.output_prefix, self.association_pack)
        extract_tool.run_tool()

        # Retrieve outputs – all tools _should_ append to the outputs object so they can be retrieved here.
        self.set_outputs(extract_tool.get_outputs())

    def _load_module_options(self) -> None:

        example_dxfile = 'file-123...'

        self._parser.add_argument('--association_tarballs',
                                  help="Path or hash to list file / single tarball of masks from "
                                       "'mergecollapsevariants'",
                                  type=self.dxfile_input, dest='association_tarballs', required=True,
                                  metavar=example_dxfile)
        self._parser.add_argument('--bgen_index',
                                  help="list of bgen files and associated index/annotation",
                                  type=self.dxfile_input, dest='bgen_index', required=True,
                                  metavar=example_dxfile)
        self._parser.add_argument('--gene_ids',
                                  help="A valid ENST Gene ID OR Gene Symbol, or space-separated list of Symbols/IDs "
                                       "to extract carriers and phenotype/covariate information for. These ID(s) MUST "
                                       "exist in the file provided to --transcript_index",
                                  type=str, dest='gene_ids', required=True, nargs='+',
                                  metavar=("GENE_ID1", "GENE_ID2"))
        self._parser.add_argument('--sparse_grm',
                                  help="A sparse GRM generated by SAIGE step0",
                                  type=self.dxfile_input, dest='sparse_grm', required=True,
                                  metavar=example_dxfile)
        self._parser.add_argument('--sparse_grm_sample',
                                  help="Sparse GRM sample file with row/column names for 'sparse_grm'",
                                  type=self.dxfile_input, dest='sparse_grm_sample', required=True,
                                  metavar=example_dxfile)

    def _parse_options(self) -> ExtractProgramArgs:
        return ExtractProgramArgs(**vars(self._parser.parse_args(self._input_args.split())))

    def _ingest_data(self, parsed_options: ExtractProgramArgs) -> ExtractAssociationPack:
        ingested_data = extract_ingester.ExtractIngestData(parsed_options)
        return ingested_data.get_association_pack()
