"""Script to test MHC sequence validity before running `setup_for_alphafold.py`."""

import click as ck 

from tcrdock.sequtil import get_mhc_class_1_alseq 

@ck.command()
@ck.option('--allele', '-a', help='MHC allele sequence')
def main(allele: str):
    """Check if the MHC allele sequence is valid."""
    if not allele:
        raise ValueError('Please provide a valid MHC allele sequence.')
    print(f'MHC allele sequence: {allele}')
    mhc_class_1 = get_mhc_class_1_alseq(allele) 
    print(f'MHC class 1 allele sequence: \n{mhc_class_1}')


    # trg_mhc_seq = trg_mhc_alseq.replace(ALL_GENES_GAP_CHAR,'')



if __name__ == '__main__':
    main()

