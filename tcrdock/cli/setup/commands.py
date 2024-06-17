"""CLI command for setting up AlphaFold modeling of TCR:pMHC complexes."""
import click as ck 
import pandas as pd
import os
import tcrdock
import sys
import argparse

from pathlib import Path

required_columns = "organism mhc_class mhc peptide va ja cdr3a vb jb cdr3b".split()
epilog = """
    The --targets-tsvfile file should contain the columns

    organism = 'mouse' or 'human'
    mhc_class = 1 or 2
    mhc = the MHC allele, e.g. 'A*02:01' or 'H2Db'
    peptide = the peptide sequence, for MHC class 2 should have length exactly 11
              (the 9 residue core plus 1 residue on either side)
    va = V-alpha gene
    ja = J-alpha gene
    cdr3a = CDR3-alpha sequence, starts with C, ends with the F/W/etc right before the
            GXG sequence in the J gene
    vb = V-beta gene
    jb = J-beta gene
    cdr3b = CDR3-beta sequence, starts with C, ends with the F/W/etc right before the
            GXG sequence in the J gene

    Check out the input files in the examples below.

Examples:

tcrdock setup --targets-tsvfile examples/benchmark/single_target.tsv \\
    --output-dir test_setup_single

or

tcrdock setup --targets-tsvfile examples/benchmark/full_benchmark.tsv \\
    --output-dir test_setup_full_benchmark --benchmark
"""

from tcrdock.cli import context_settings 

context_settings.update(
    show_default=True,
    help_option_names=["-h", "--help"],
)

@ck.command(
        epilog=epilog,
        context_settings=context_settings, 
)
@ck.option(
    "-t", "--targets-tsvfile", 
    required=True,
    help="TSV formatted file with info on modeling targets.",
    type=ck.Path(exists=True, file_okay=True, readable=True,
        resolve_path=True, path_type=Path), 
)
@ck.option(
    "-o", "--output-dir", 
    required=True,
    help="directory where the outputs will go; preferably empty "
        "or nonexistent (will create).",
    type=ck.Path(exists=False, file_okay=False, writable=True, 
        resolve_path=True, path_type=Path), 
)
@ck.option(
    "-n", "--num-runs", 
    default=3,
    help="Number of AlphaFold runs per target.",
    type=ck.INT,
    show_default=True,
)
@ck.option(
    "--benchmark/--no-benchmark", "-B", default=False, 
    help="Exclude sequence-similar templates (for benchmarking).",
    type=ck.BOOL,
)
@ck.option(
    "--maintain-relative-paths/--absolute-paths", "-R", default=False, 
    help="Keep the file paths in the targets and alignments files "
        "relative (assuming --output_dir is a relative path), rather "
        "than trying to resolve to absolute paths.",
    type=ck.BOOL,
)
@ck.option(
    "--exclude-pdbids-column", 
    help="Column in the --targets-tsvfile file with comma-separated "
        "lists of PDB files to exclude from modeling.",
    type=ck.STRING,
)
@ck.option(
    "--new-docking/--old-docking", "-N", default=False, 
    help="Use a new, unpublished approach for constructing the "
        "TCR:pMHC docking geometry in the AlphaFold templates. "
        "This only requires 1 run per target (ie, it is 3x faster than "
        "the default) and preliminary tests of the fine-tuned version "
        "are promising.",
    type=ck.BOOL,
)
def setup(
    targets_tsvfile: Path, 
    output_dir: Path,
    benchmark: bool,
    maintain_relative_paths: bool,
    exclude_pdbids_column: str,
    new_docking: bool,
):
    """
    Create <OUTPUT_DIR>/targets.tsv file and associated input files for AlphaFold modeling with the `tcrdock run` command.
    
    \f 

    Parameters
    ----------
    targets_tsvfile : Path
        Path to the TSV formatted file with info on modeling targets.
    output_dir : Path
        Path to the directory where the outputs will go.
    benchmark : bool
        Exclude sequence-similar templates (for benchmarking).
    maintain_relative_paths : bool
        Keep the file paths in the targets and alignments files relative (assuming --output_dir is a relative path), rather than trying to resolve to absolute paths.
    exclude_pdbids_column : str
        Column in the --targets_tsvfile file with comma-separated lists of PDB files to exclude from modeling.
    new_docking : bool
        Use a new, unpublished approach for constructing the TCR:pMHC docking geometry in the AlphaFold templates.

    Returns
    -------
    None
    """
    # load
    print(type(targets_tsvfile))
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
