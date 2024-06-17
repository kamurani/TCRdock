parser = argparse.ArgumentParser(
    description = "Create <output_dir>/targets.tsv file and associated input files "
    "for AlphaFold modeling with the run_alphafold_predictions.py script",
    epilog = f'''The --targets_tsvfile file should contain the columns

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

    

Example usage:

tcrdock setup --targets-tsvfile examples/benchmark/single_target.tsv \\
    --output-dir test_setup_single

or

python setup_for_alphafold.py --targets_tsvfile examples/benchmark/full_benchmark.tsv \\
    --output_dir test_setup_full_benchmark --benchmark

''',
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

parser.add_argument('--targets_tsvfile', required=True,
                    help='TSV formatted file with info on modeling targets')
parser.add_argument('--output_dir', required=True,
                    help='directory where the outputs will go; preferably empty '
                    'or nonexistent (will create)')
parser.add_argument('--num_runs', type=int, default=3,
                    help='Number of alphafold runs per target (default is 3)')
parser.add_argument('--benchmark', action='store_true',
                    help='Exclude sequence-similar templates (for benchmarking)')
parser.add_argument('--maintain_relative_paths', action='store_true',
                    help='Keep the file paths in the targets and alignments files '
                    'relative (assuming --output_dir is a relative path), rather '
                    'than trying to resolve to absolute paths')
parser.add_argument('--exclude_pdbids_column',
                    help='Column in the --targets_tsvfile file with comma-separated '
                    'lists of pdbfiles to exclude from modeling')
parser.add_argument('--new_docking', action='store_true',
                    help='Use a new, unpublished approach for constructing the '
                    'TCR:pMHC docking geometry in the AlphaFold templates. '
                    'This only requires 1 run per target (ie, it is 3x faster than '
                    'the default) and preliminary tests of the fine-tuned version '
                    'are promising.')

args = parser.parse_args()
