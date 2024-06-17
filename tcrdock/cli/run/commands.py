"""CLI command for running AlphaFold inference on a set of targets."""
from typing import List
import click as ck 
import numpy as np
import pandas as pd
import os
import sys
import itertools
import tcrdock

from os.path import exists
from pathlib import Path

epilog = """ 
    This script is borrowed from the original here:
    https://github.com/phbradley/alphafold_finetune/blob/main/run_prediction.py

    Examples:

    \b
    # The test_setup_single/ directory would first be made by running the
    # `tcrdock setup` command.
    
    \b
    tcrdock run --targets test_setup_single/targets.tsv 
        --outfile-prefix test_run_single --model-names model_2_ptm 
        --data-dir $ALPHAFOLD_DATA_DIR
"""
from tcrdock.cli import context_settings
context_settings.update(
    show_default=True,
    help_option_names=["-h", "--help"],
)
@ck.command(epilog=epilog, context_settings=context_settings)
@ck.option(
    "-t", "--targets", 
    required=True,
    help="File listing the targets to "
        "be modeled. See description of file format in the README "
        "https://github.com/phbradley/alphafold_finetune/blob/main/README.md"
        "and also examples in that repo's examples/*/*tsv",
    type=ck.Path(exists=True, file_okay=True, readable=True,
        resolve_path=True, path_type=Path), 
)
@ck.option(
    "--outfile-prefix",
    help="Prefix that will be prepended to the output filenames.",
    type=ck.STRING,
)
@ck.option(
    "--final-outfile-prefix",
    help="Prefix that will be prepended to the final output .tsv filename.",
    type=ck.STRING,
)
@ck.option(
    "--data-dir",
    help="Location of AlphaFold params/ folder.",
    type=ck.Path(exists=True, file_okay=False, readable=True,
        resolve_path=True, path_type=Path),
)
@ck.option(
    "--model-names", 
    multiple=True, help="Model names to use.", type=ck.STRING,
    default=['model_2_ptm'],
)
@ck.option(
    "--model-params-files", 
    multiple=True, help="Only needed if running with fine-tuned parameters or "
        "parameters in a non-default location (ie, not in the params/ "
        "folder in --data_dir).", type=ck.STRING,
)
@ck.option(
    "--verbose", is_flag=True, help="Print verbose output.",
)
@ck.option(
    "--ignore-identities", is_flag=True, help="Ignore the sequence identities column in the templates "
        "alignment files. Useful when modeling many different peptides "
        "using the same alignment file.",
)
@ck.option(
    "--write-pdbs/--no-pdbs", default=True, help="Write out pdbs.",
)
@ck.option(
    "--terse", is_flag=True, help="Dont write out pdbs or "
        "matrices with AlphaFold confidence values.",
)
@ck.option(
    "--resample-msa/--no-resample-msa", default=True, help="Randomly "
        "resample from the MSA during recycling. Perhaps useful for "
        "testing...",
)
def run(
    targets: Path,
    outfile_prefix: str,
    final_outfile_prefix: str,
    data_dir: Path,
    model_names: List, 
    model_params_files: List,
    verbose: bool,
    ignore_identities: bool,
    write_pdbs: bool,
    terse: bool,
    resample_msa: bool,
    max_content_width: int = 120,

    **params, 
) -> None:
    """
    Run simple template-based AlphaFold inference.

    \f
    Parameters
    ----------
    max_content_width : int
        Maximum width of the terminal rewrapping. 

    Returns
    -------
    None
    """
    import predict_utils

    targets = pd.read_table(targets)
    lens = [len(x.target_chainseq.replace('/',''))
            for x in targets.itertuples()]
    crop_size = max(lens)

    if verbose:
        import jax
        from os import popen # just to get hostname for logging, not necessary
        # print some logging info
        platform = jax.local_devices()[0].platform
        hostname = popen('hostname').readlines()[0].strip()

        print('cmd:', ' '.join(sys.argv))
        print('local_device:', platform, 'hostname:', hostname, 'num_targets:',
            targets.shape[0], 'max_len=', crop_size)

    sys.stdout.flush()

    model_runners = predict_utils.load_model_runners(
        model_names,
        crop_size,
        data_dir,
        model_params_files=model_params_files,
        resample_msa_in_recycling = resample_msa,
    )

    final_dfl = []
    for counter, targetl in targets.iterrows():
        print('START:', counter, 'of', targets.shape[0])

        alignfile = targetl.templates_alignfile
        assert exists(alignfile)

        query_chainseq = targetl.target_chainseq
        if 'outfile_prefix' in targetl:
            outfile_prefix = targetl.outfile_prefix
        else:
            assert outfile_prefix is not None
            if 'targetid' in targetl:
                outfile_prefix = outfile_prefix+'_'+targetl.targetid
            else:
                outfile_prefix = f'{outfile_prefix}_T{counter}'

        query_sequence = query_chainseq.replace('/','')
        num_res = len(query_sequence)

        data = pd.read_table(alignfile)
        cols = ('template_pdbfile target_to_template_alignstring identities '
                'target_len template_len'.split())
        template_features_list = []
        for tnum, row in data.iterrows():
            #(template_pdbfile, target_to_template_alignstring,
            # identities, target_len, template_len) = line[cols]

            assert row.target_len == len(query_sequence)
            target_to_template_alignment = {
                int(x.split(':')[0]) : int(x.split(':')[1]) # 0-indexed
                for x in row.target_to_template_alignstring.split(';')
            }

            template_name = f'T{tnum:03d}' # dont think this matters
            template_features = predict_utils.create_single_template_features(
                query_sequence, row.template_pdbfile, target_to_template_alignment,
                template_name, allow_chainbreaks=True, allow_skipped_lines=True,
                expected_identities = None if ignore_identities else row.identities,
                expected_template_len = row.template_len,
            )
            template_features_list.append(template_features)

        all_template_features = predict_utils.compile_template_features(
            template_features_list)

        msa=[query_sequence]
        deletion_matrix=[[0]*len(query_sequence)]

        all_metrics = predict_utils.run_alphafold_prediction(
            query_sequence=query_sequence,
            msa=msa,
            deletion_matrix=deletion_matrix,
            chainbreak_sequence=query_chainseq,
            template_features=all_template_features,
            model_runners=model_runners,
            out_prefix=outfile_prefix,
            crop_size=crop_size,
            dump_pdbs = not (not write_pdbs or terse),
            dump_metrics = not terse,
        )


        outl = targetl.copy()
        for model_name, metrics in all_metrics.items():
            plddts = metrics['plddt']
            paes = metrics.get('predicted_aligned_error', None)
            filetags = 'pdb plddt ptm predicted_aligned_error'.split()
            for tag in filetags:
                fname = metrics.get(tag+'file', None)
                if fname is not None:
                    outl[f'{model_name}_{tag}_file'] = fname

            cs = query_chainseq.split('/')
            chain_stops = list(itertools.accumulate(len(x) for x in cs))
            chain_starts = [0]+chain_stops[:-1]
            nres = chain_stops[-1]
            assert nres == num_res
            outl[model_name+'_plddt'] = np.mean(plddts[:nres])
            if paes is not None:
                outl[model_name+'_pae'] = np.mean(paes[:nres,:nres])
            for chain1,(start1,stop1) in enumerate(zip(chain_starts, chain_stops)):
                outl[f'{model_name}_plddt_{chain1}'] = np.mean(plddts[start1:stop1])

                if paes is not None:
                    for chain2 in range(len(cs)):
                        start2, stop2 = chain_starts[chain2], chain_stops[chain2]
                        pae = np.mean(paes[start1:stop1,start2:stop2])
                        outl[f'{model_name}_pae_{chain1}_{chain2}'] = pae
        final_dfl.append(outl)

    if final_outfile_prefix:
        outfile_prefix = final_outfile_prefix
    elif outfile_prefix:
        outfile_prefix = outfile_prefix
    elif 'outfile_prefix' in targets.columns:
        outfile_prefix = targets.outfile_prefix.iloc[0]
    else:
        outfile_prefix = None

    if outfile_prefix:
        outfile = f'{outfile_prefix}_final.tsv'
        pd.DataFrame(final_dfl).to_csv(outfile, sep='\t', index=False)
        print('made:', outfile)