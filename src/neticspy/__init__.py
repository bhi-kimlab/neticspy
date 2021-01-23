__version__ = '0.0.5'

import argparse
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

warnings.filterwarnings(action='ignore', category=RuntimeWarning)

def parse_args():
	parser = argparse.ArgumentParser(description='Python implementation of NetICS')

	parser.add_argument(
		'-a',
		'--aberration',
		required=True,
		type=str,
		help='Input two-column table (without headers) containing genetically aberrant genes for each sample. It contain two columns that map every gene (1st column) to the samples that it it genetically aberrant (2nd column).'
	)
	parser.add_argument(
		'-j',
		'--adj',
		required=True,
		type=str,
		help='Adjacency matrix of the directed interaction network.'
	)
	parser.add_argument(
		'-n',
		'--network',
		required=True,
		help='Input file that contains the list of the genes that are present in the network.\nThey should be in the same order as in the rows of the adjacency matrix adj. An example file is given that contains the gene names of the network described in (Wu et al., 2010).'
	)
	parser.add_argument(
		'-d',
		'--degs',
		default=None,
		help='List of the names of predefined differentially expressed genes.',
	)
	parser.add_argument(
		'-o',
		'--output',
		required=True,
		help='File to save NetICS result.'
	)
	parser.add_argument(
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
	parser.add_argument(
		'-f',
		'--diffusion-matrix',
		default=None,
		help='(Optional) Path to .npy file for diffusion matrix.'
	)
	parser.add_argument(
		'-g',
		'--opposite-diffusion-matrix',
		default=None,
		help='(Optional) Path to .npy file for diffusion matrix in opposite direction.'
	)

	return parser.parse_args()

def netics_fun(
		filename_aberration,
		filename_adj,
		filename_net,
		output,
		degs=None,
		filename_F=None,
		filename_F_opposite=None,
		restart_prob=0.4,
		# rank_method='SUM',
	):
	# Accept arguments and check input.
	# if rank_method not in ['SUM', 'MEDIAN', 'RRA']:
	# 	print('Wrong rank method: %s' % rank_method)
	# 	print('Rank method should be MEDIAN, RRA or SUM')
	# 	sys.exit(1)

	unique_samples, mutation_data = read_mutations(filename_aberration)
	adj = np.loadtxt(open(filename_adj), delimiter='\t')

	choose_mut_diff = DIRECTION.DOWN if degs is None else DIRECTION.BOTH

	# Read network genes, line by line.
	network_genes = [l.strip().upper() for l in open(filename_net).readlines()]
	network_gene_set = set(network_genes)

	# Copy data from input mutation_data.
	samples = defaultdict(list)
	for i in range(len(mutation_data)):
		samples['DNA'].append(mutation_data[i])
		samples['DNA_network'].append(np.char.upper(np.intersect1d(network_genes, samples['DNA'][i]).astype(str)))
	del mutation_data

	# Read DEG names.
	degs_df = pd.read_csv(degs, sep='\t', names=['gene', 'sample'])
	degs_df = degs_df[degs_df.gene.isin(network_gene_set)]

	diff_expr_genes = []
	for sample in unique_samples:
		diff_expr_genes.append(degs_df[degs_df['sample'] == sample].gene.str.upper().values)

	for i in range(len(samples['DNA'])):
		samples['RNA'].append(diff_expr_genes[i])
		# samples['RNA_network'].append(np.intersect1d(network_genes, samples['RNA'][i]).astype(str))
		samples['RNA_network'].append(diff_expr_genes[i])

		# print(len(np.intersect1d(network_genes, diff_expr_genes).astype(str)))
		# print(samples['RNA_network'][-1])

	# Load or compute diffusion matrix.
	if filename_F is not None and filename_F_opposite is not None:
		F = np.load(filename_F)
		F_opposite = np.load(filename_F_opposite)
	else:
		print('Computing diffused matrix...')
		F = diffusion.insulated_diff(util.row_normalize(adj), restart_prob)
		F_opposite = diffusion.insulated_diff(util.row_normalize(adj.conj().transpose()), restart_prob)

	print('Running NetICS...')
	# ranked_list_genes, scores = prioritization(samples, F, F_opp, network_genes, choose_mut_diff, rank_method, unique_samples)
	final_data = prioritization(samples, F, F_opposite, network_genes, choose_mut_diff, unique_samples)

	print(f'Saving output to {output}...')
	# res = pd.DataFrame({
		# 'Genes': ranked_list_genes,
		# 'Scores': scores,
	# })
	# res.to_csv(output, sep='\t', header=True, index=False)
	final_data.to_csv(output, index=False)
	return final_data

def read_mutations(filename):
	g = np.loadtxt(open(filename, 'r'), delimiter='\t', dtype=str)

	unique_samples = np.unique(g[:, 1])
	mutation_data = []
	for sample in unique_samples:
		mutation_data.append(np.char.upper(g[:, 0][np.isin(g[:, 1], sample)]))

	return unique_samples, mutation_data

def prioritization(samples, F, F_opposite, network_genes, choose_up_down_flag, unique_samples):
	num_samples, num_genes = len(samples['DNA_network']), len(network_genes)

	final_data = []
	for i, sample in enumerate(unique_samples): # Number of samples
		tmp_result = {
			'Sample': [sample] * num_genes,
			'Gene': network_genes,
		}
		aberrant_genes, de_genes = samples['DNA_network'][i], samples['RNA_network'][i]

		if len(aberrant_genes) == 0:
			flag = 2  # If there is no aberrant genes, just diffuse up from DE genes.
		else:
			flag = choose_up_down_flag

		diffusion_scores = diffusion.diffuse_all(flag, aberrant_genes, de_genes, network_genes, F, F_opposite)
		print(sample, len(samples['DNA_network'][i]), len(samples['RNA_network'][i]))
		print(diffusion_scores.flatten().sum())
		tmp_result['Diffusion_score'] = diffusion_scores.flatten()

		tmp_result = pd.DataFrame(tmp_result)
		tmp_result['Rank'] = tmp_result['Diffusion_score'].rank(ascending=False)

		final_data.append(tmp_result)

		# weights_all[i, :] = diffuse_all(choose_up_down_flag, samples['DNA_network'][i], samples['RNA_network'][i], network_genes, F, F_opp)
		# to_sort_sorted = sortrows(np.array([weights_all[i, :], np.arange(len(weights_all[i,:]))]).conj().transpose())
		# positions_sorted[:, i] = np.flip(to_sort_sorted[:, 1].conj().transpose())

	final_data = pd.concat(final_data)
	return final_data

	# for j in range(num_genes): # Number of genes
	# 	rank = np.zeros([num_samples, 1])

	# 	for i in range(num_samples): # Number of samples
	# 		#rank[i] = (positions_sorted[:,i] == j).nonzero()
	# 		rank[i] = (positions_sorted[:, i] == j).nonzero()[0] + 1

	# 	final_weights_median[j] = np.median(rank)
	# 	final_weights_sum[j] = rank.sum(axis=0)
	# 	final_weights_rho[j] = rhoScores(rank / max(rank))

	# # Choose the aggregation method
	# if rank_method == 'MEDIAN':
	# 	scores = final_weights_median.astype(np.int64)
	# elif rank_method == 'SUM':
	# 	scores = final_weights_sum.astype(np.int64)
	# else:
	# 	scores = final_weights_rho

	# # Sort scores
	# to_sort_sorted_rho = sortrows(np.array([scores, np.arange(len(scores))]).conj().transpose())
	# ranked_list_genes = network_genes[to_sort_sorted_rho[:,1].astype(np.int64)]
	# ranked_scores = scores[to_sort_sorted_rho[:,1].astype(np.int64)]
	# return ranked_list_genes, ranked_scores

def sortrows(M):
	if len(M.shape) != 2:
		raise ValueError('M must be 2d numpy.array')
	M_columns = tuple(M[:,c] for c in range(M.shape[1]-1, -1, -1))
	return M[np.lexsort(M_columns), :]

def rhoScores(r, topCutoff=None):
	r = r.reshape((1,-1))

	rho = np.empty(1)
	rho[:] = np.nan

	r1 = r[0,:]
	x = betaScores(r1)[0]
	# Correct using Bonferroni method.
	rho[0] = correctBetaPvalues(min(x), sum(~np.isnan(x)))
	return rho

def betaScores(r):
	n = sum(~np.isnan(r))
	p = np.empty((1,len(r)))
	p[:] = np.nan
	# Sort the values
	r.sort()
	# Get the order statistics and calculates p-values for each of the order statistics. These are based on their expected distribution under the null hypothesis of uniform distribution.
	p[:,:n] = beta.cdf(r[0:n], np.arange(1, n+1), np.arange(n, 0, -1))
	return p

def correctBetaPvalues(p, k):
	pval = beta.cdf(p, 1, k)
	return pval

def main():
	args = parse_args()
	netics_fun(args.aberration, args.adj, args.network, args.output, args.degs, args.diffusion_matrix, args.opposite_diffusion_matrix, args.beta)

if __name__ == '__main__':
	main()

