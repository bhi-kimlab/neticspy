__version__ = '0.1.3'

import argparse
import cleanlog
import warnings
import sys
import time
import numpy as np
import pandas as pd

from . import util
from . import diffusion

from collections import defaultdict
from .const import DIRECTION

from scipy.stats import beta
from scipy.stats import combine_pvalues

warnings.filterwarnings(action='ignore', category=RuntimeWarning)

logger = cleanlog.ColoredLogger('NetICS')
logger.setLevel(cleanlog.INFO)

def parse_args():
    parser = argparse.ArgumentParser(description='Python implementation of NetICS')

    subparsers = parser.add_subparsers(dest='command', help='Subcommands.')

    subparser_diffuse = subparsers.add_parser('diffuse', help='Prepare diffusion matrix for the network.')
    subparser_diffuse.add_argument(
        '-j',
        '--adj',
        required=True,
        type=str,
        help='Adjacency matrix of the directed interaction network.'
    )
    subparser_diffuse.add_argument(
        '-b',
        '--beta',
        type=float,
        default=0.4,
        help='Restart probability for the insulated diffusion. Default: 0.4 (For the network from Wu et al., 2010)'
    )
    subparser_diffuse.add_argument(
        '-o',
        '--output',
        required=True,
        help='Output filename for diffusion matrix in .npz format.'
    )

    subparser_rank = subparsers.add_parser('rank', help='Run NetICS algorithm and rank genes.')
    subparser_rank.add_argument(
        '-a',
        '--aberration',
        required=True,
        type=str,
        help='Input two-column table (without headers) containing genetically aberrant genes for each sample. It contain two columns that map every gene (1st column) to the samples that it it genetically aberrant (2nd column).'
    )
    subparser_rank.add_argument(
        '-f',
        '--diffusion-matrix',
        required=True,
        help='Path to .npz file for diffusion matrix.'
    )
    subparser_rank.add_argument(
        '-n',
        '--network',
        required=True,
        help='Input file (without headers) that contains the list of the genes that are present in the network.\nThey should be in the same order as in the rows of the adjacency matrix.'
    )
    subparser_rank.add_argument(
        '-d',
        '--degs',
        default=None,
        help='Two-column table (without headers) with the names of predefined differentially expressed genes and the corresponding samples.',
    )
    subparser_rank.add_argument(
        '-o',
        '--output-prefix',
        required=True,
        help='Prefix of the output file to save raw NetICS result and aggregated ranks.'
    )
    subparser_rank.add_argument(
        '-b',
        '--beta',
        type=float,
        default=0.4,
        help='Restart probability for the insulated diffusion. Default: 0.4 (For the network from Wu et al., 2010)'
    )
    # parser.add_argument(
    # 	'-r',
    # 	'--rank',
    # 	default='SUM',
    # 	help="'MEDIAN' uses the median of the sample-specific ranks.\n'RRA' uses the Robust Rank Aggregation method to integrate sample-specific ranked lists.\n'SUM' uses the sum of the sample-specific ranks."
    # )
    subparser_rank.add_argument(
        '-v',
        '--verbose',
        default=False,
        action='store_true',
        help='Print debug messages'
    )

    return parser.parse_args()

def diffuse(filename_adj, restart_prob, output):
    logger.info('Making diffusion matrix...')

    logger.info('Started making forward diffusion matrix...')
    adj = np.loadtxt(open(filename_adj), delimiter='\t')
    F = diffusion.insulated_diff(util.row_normalize(adj), restart_prob)
    logger.info('Done!')

    logger.info('Started making backward diffusion matrix...')
    F_opposite = diffusion.insulated_diff(util.row_normalize(adj.conj().transpose()), restart_prob)
    logger.info('Done!')

    np.savez(output, forward=F, backward=F_opposite)
    logger.info(f'Successfully saved diffusion matrix to {output}.')

def netics_fun(
        filename_aberration,
        filename_genes,
        output,
        filename_diffusion_matrix,
        filename_deg_list=None,
        verbose=False,
        # rank_method='SUM',
    ):
    if verbose:
        logger.setLevel(cleanlog.DEBUG)
    # Accept arguments and check input.
    # if rank_method not in ['SUM', 'MEDIAN', 'RRA']:
    # 	print('Wrong rank method: %s' % rank_method)
    # 	print('Rank method should be MEDIAN, RRA or SUM')
    # 	sys.exit(1)

    # unique_samples, mutation_data = read_mutations(filename_aberration)

    # Read network genes, line by line.
    network_genes = [l.strip().upper() for l in open(filename_genes).readlines()]
    gene2idx = {g:i for i, g in enumerate(network_genes)}
    network_gene_set = set(network_genes)

    # Aberrations.
    mutation_df = read_mutations(filename_aberration)
    mutation_df['idx'] = mutation_df.gene.map(gene2idx)
    mutation_df = mutation_df.dropna()

    # DEGs.
    deg_df = pd.read_csv(filename_deg_list, sep='\t', names=['gene', 'sample'])
    deg_df = deg_df[deg_df.gene.isin(network_gene_set)]
    deg_df['idx'] = deg_df.gene.map(gene2idx)
    deg_df = deg_df.dropna()

    # Determine the direction of the diffusion.
    choose_mut_diff = DIRECTION.DOWN if filename_deg_list is None else DIRECTION.BOTH

    # Load or compute diffusion matrix.
    diffusion_matrices = np.load(filename_diffusion_matrix)
    F, F_opposite = diffusion_matrices['forward'], diffusion_matrices['backward']

    logger.info('Running NetICS...')
    final_result = []
    for sample in mutation_df['sample'].unique():
        mutation_df_per_sample = mutation_df[mutation_df['sample'] == sample]
        deg_df_per_sample = deg_df[deg_df['sample'] == sample]

        result = prioritization(sample, mutation_df_per_sample, deg_df_per_sample, F, F_opposite, network_genes, choose_mut_diff)
        final_result.append(result)

    final_result = pd.concat(final_result)

    logger.info(f'Saving output to {output}.raw.txt...')
    final_result.to_csv(f'{output}.raw.csv', index=False)

    # Rank aggregation.
    rank_agg_result = final_result.pivot_table(values='rank', index='gene', aggfunc=['mean', 'median'])
    rank_agg_result.columns = ['rank_mean', 'rank_median']
    rank_agg_result.sort_values('rank_mean').to_csv(f'{output}.rank_aggregated.csv')

    return final_result

def read_mutations(filename):
    return pd.read_csv(filename, sep='\t', names=['gene', 'sample'])

def prioritization(sample, mutation_df, deg_df, F, F_opposite, network_genes, choose_up_down_flag):
    num_genes = len(network_genes)

    result = {
        'sample': [sample] * num_genes,
        'gene': network_genes,
    }
    aberrant_gene_idx, deg_idx = mutation_df.idx.astype(int).values, deg_df.idx.astype(int).values
    logger.debug(f'{len(aberrant_gene_idx)}, {len(deg_idx)}')

    if len(aberrant_gene_idx) == 0:
        flag = DIRECTION.UP
    else:
        flag = choose_up_down_flag

    logger.info(f'Computing diffusion scores for {sample}.')
    diffusion_score = diffusion.diffuse_all(flag, aberrant_gene_idx, deg_idx, F, F_opposite)[0]
    logger.debug(f'{sample}, {len(mutation_df)}, {len(deg_df)}, {flag}')

    result['diffusion_score'] = diffusion_score
    result = pd.DataFrame(result)
    result['rank'] = result['diffusion_score'].rank(ascending=False)

    return result.sort_values('rank')

def main():
    args = parse_args()

    if args.command == 'diffuse':
        diffuse(
            filename_adj=args.adj,
            restart_prob=args.beta,
            output=args.output
        )
    else:
        netics_fun(
            args.aberration,
            args.network,
            args.output_prefix,
            args.diffusion_matrix,
            args.degs,
            verbose=args.verbose,
        )

if __name__ == '__main__':
    main()

