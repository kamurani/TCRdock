"""CLI command for setting up AlphaFold modeling of TCR:pMHC complexes."""
import click as ck 
import pandas as pd
import os
import tcrdock
import sys
import argparse

from pathlib import Path

required_columns = "organism mhc_class mhc peptide va ja cdr3a vb jb cdr3b".split()


@ck.command()
@ck.option(
    "-t", "--targets-tsvfile", 
    #"--targets_tsvfile", 
    required=True,
    help="TSV formatted file with info on modeling targets",
    type=ck.Path(exists=True, file_okay=True, readable=True,
        resolve_path=True, path_type=Path), 
)
def setup(
    targets_tsvfile, 
):
    """CLI command for setting up AlphaFold modeling of TCR:pMHC complexes.
    
    """
    # load
    targets = pd.read_table(targets_tsvfile)
    print(targets.head())


    exit()


    if args.benchmark:
        exclude_self_peptide_docking_geometries = True
        min_single_chain_tcrdist = 36
        min_pmhc_peptide_mismatches = 3
        min_dgeom_peptide_mismatches = 3
        min_dgeom_paired_tcrdist = 48.5
        min_dgeom_singlechain_tcrdist = 0.5
    else:
        exclude_self_peptide_docking_geometries = False
        min_single_chain_tcrdist = -1
        min_pmhc_peptide_mismatches = -1
        min_dgeom_peptide_mismatches = -1
        min_dgeom_paired_tcrdist = -1
        min_dgeom_singlechain_tcrdist = -1

    # load the modeling targets
    targets = pd.read_table(targets_tsvfile)
    missing = [col for col in required_columns if col not in targets.columns]
    if missing:
        print('ERROR --targets_tsvfile is missing required columns:', missing)
        print('see --help message for details on required columns')
        exit()

    # check the gene names
    ok = tcrdock.sequtil.check_genes_for_modeling(targets)

    if not ok:
        print(f'ERROR some of the genes in {args.targets_tsvfile} are problematic,'
            ' see error messages above')
        sys.exit()

    # make the output dir
    output_dir = args.output_dir
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    if not args.maintain_relative_paths:
        # since we are saving filenames to the targets/templates files, we want these
        # to work from wherever AlphaFold is eventually run, so try to make them
        # absolute file paths
        output_dir = str(Path(output_dir).resolve())

    if not output_dir.endswith('/'):
        output_dir += '/' # this will cause issues on windows but code later assumes it...

    num_runs = 1 if args.new_docking else args.num_runs

    # run the setup code
    tcrdock.sequtil.setup_for_alphafold(
        targets, output_dir, clobber=True,
        min_single_chain_tcrdist=min_single_chain_tcrdist,
        exclude_self_peptide_docking_geometries=exclude_self_peptide_docking_geometries,
        min_pmhc_peptide_mismatches=min_pmhc_peptide_mismatches,
        num_runs = num_runs,
        min_dgeom_peptide_mismatches=min_dgeom_peptide_mismatches, #
        min_dgeom_paired_tcrdist = min_dgeom_paired_tcrdist, #
        min_dgeom_singlechain_tcrdist = min_dgeom_singlechain_tcrdist,
        exclude_pdbids_column = args.exclude_pdbids_column,
        use_opt_dgeoms = args.new_docking,
    )
