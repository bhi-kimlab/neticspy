__version__ = '0.0.0'

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from scipy.special import gammainc

def insulated_diff( W, b ):
	temp = np.identity(W.shape[1]) - (1-b)*W
	F = b*np.linalg.inv(temp)
	return F

def norm_adj( adf ):
	degree = adj.conj().transpose().sum(axis=0)
	W = np.zeros((adj[:,0].shape[1],adj[0,].shape[1]))
	for i in range(adj[:,0].shape[1]):
		for j in range(adj[0,].shape[1]):
			if adj[i,j] != 0:
				W[i,j] = 1/degree[i]
	return W

def np_sortrows(M):
	if len(M.shape) != 2:
		raise ValueError('M must be a 2d numpy.array')
	M_columns = tuple(M[:,c] for c in range(M.shape[1]-1, -1, -1))
	return M[np.lexsort(M_columns), :]

def prioritization( samples, F, F_opp, network_genes, choose_up_down_flag, rank_method_str):
	weights_all = np.zeros((samples_DNA_network.shape[1],network_genes.shape[1]))
	positions_sorted = np.zeros((network_genes.shape[1],samples_DNA_network.shape[1]))
	final_weights_median, final_weights_sum, final_weights_rho = np.zeros((1,network_genes.shape[1]))

	for i in range(samples_DNA_network.shape[1]): # number of samples
		weights_all[i,] = diffuse_all(choose_up_down_flag, samples_DNA_network[i], samples_RNA_network[i], network_genes, F, F_opp)
		to_sorted_sorted = np_sortrows(np.array([weights_all[i,],np.arange(weights_all[i,].shape[1])]).conj().transpose())
		positions_sorted[:,i] = np.fliplr(to_sort_sorted[:,1].conj().transpose()).conj().transpose()
	for j in range(network_genes.shape[1]): # number of genes
		rank = np.zeros((samples_DNA_network.shape[1],1))
		for i in range(samples_DNA_network.shape[1]): # number of samples
			rank[i] = (positions_sorted[:,i] == j).nonzero()
		final_weights_median[j] = np.median(rank)
		final_weights_sum[j] = rank.sum(axis=0)
		#final_weights_rho[j] = ?????

	# choose the aggregation method
	if rank_method_str == 'MEDIAN':
		scores = final_weights_median
	elif rank_method_str == 'SUM':
		scores = final_weights_sum
	#else:
		#scores = final_weights_rho

	# sort scores
	to_sort_sorted_rho = np_sortrows(np.array([scores, np.arange(scores.shape[1])]).conj().transpose())
	ranked_list_genes = network_genes[to_sort_sorted_rho[:,1]]
	ranked_scores = scores[to_sort_sorted_rho[:,1]]

	return ranked_list_genes, ranked_scores

def diffuse_all(choose_up_down_flag, sample_DNA, sample_RNA, network_genes, F, F_opp):
	# diffuse mutation labels
	if choose_up_down_flag == 1 or choose_up_down_flag == 3:
		mutation_weights = diffuse(sample_DNA.intersection(network_genes), network_genes, F)
	# diffuse DE labels
	if choose_up_down_flag == 2 or choose_up_down_flag == 3:
		DE_weights = diffuse(sample_RNA.intersection(network_genes), network_genes, F_opp)
	# combines scores by multiplication (use both mutations and DE)
	if choose_up_down_flag == 3:
		ret_weights_all = mutation_weights * DE_weights
	elif choose_up_down_flag == 1:
		ret_weights_all = mutation_weights
	else:
		ret_weights_all = DE_weights
	return ret_weights_all

def diffuse(inter_mut_sampl_network, network_genes, F):
	positions = (network_genes in inter_mut_sampl_network).ravel().nonzero()
	mutation_weights = (1/positions.shape[1])*np.ones((1,positions.shape[1]))*F[positions-1,:]
	return ret_weights_all

def np_isempty(M):
	if len(M) == 0:
		return 0
	else:
		return 1

def netics_fun(filenameMu=None, adj=None, restart_prob=0.4, rank_method_str='MEDIAN', filenameNet=None, filenameRNA=None, filenamePR=None):
	filenameRNA = ''
	filenamePR = ''
	ranked_list_genes = np.array([])

	# Accept arguments
	if None in [filenameMu, adj, filenameNet]:
		raise ValueError('Error: Missing input arguments.')
	else:
		mutation_data = read_mutations(filenameMu)
		choose_mut_diff = 1
	if filenameRNA is not None:
		choose_mut_diff = 3

	# read network genes
	fidnet = open(filenameNet, 'r')
	g = np.loadtxt(fidnet, dtype=str, delimiter='\n')
	network_genes = g
	fidnet.close()

	# copy data from input mutation_data
	samples = np.array([])
	for i in range(mutation_data.shape[1]):
		samples_DNA[i] = np.char.upper(mutation_data[i])
		samples_DNA_network[i] = np.char.upper(network_genes.intersection(samples_DNA[i]))
	del mutation_data

	# read differentially expressed genes
	diff_expr_genes = np.array([])
	if choose_mut_diff != 1:
		diff_expr_genes = read_diff_expr(filenameRNA, filenamePR)
		if np_isempty(diff_expr_genes):
			choose_mut_diff = 1
	else:
		print('No input files were provided for differentially expressed genes.')
	for i in range(samples_DNA.shape[1]):
		samples_RNA[i] = diff_expr_genes
		samples_RNA_networks[i] = np.char.upper(network_genes.intersection(samples_RNA[i]))
	
	# compute diffused matrix
	print('Computing diffused matrix...')
	F = insulated_diff(norm_adj(adj), restart_prob)
	F_opp = insulated_diff(norm_adj(adj.conj().transpose()), restart_prob)
	print('Running NetICS...')
	ranked_list_genes, scores = prioritization(samples, F, F_opp, network_genes, choose_mut_diff, rank_method_str)
	return ranked_list_genes, scores

def read_mutations(filename):
	fid = open(filename, 'r')
	g = np.loadtxt(fid, dtype=str, delimiter='\t')
	fid.close()
	unique_samples = np.unique(g[:,1])
	mutation_data = np.empty([unique_samples.shape[1],1])
	for i in range(uniq_samples.shape[1]):
		mutation_data[i] = g[:,0][g[:,1] in uniq_samples[i]]
	return mutation_data

def pchisq(x,a):
# PCHISQ	The chisquare distribution function
#	F = pchisq(x, DegreeOfFreedom)
	if a <= 0:
		raise ValueError('Degrees Of Freedom is wrong')
	F = pgamma(x/2, a*0.5)
	return F

def pgamma(x,a):
#PGAMMA	The gamma distribution function
#	F = pgamma(x,a)
	if a<= 0:
		raise ValueError('Parameter a is wrong')
	F = gammainc(a,x)
	I0 = (x<0).nonzero()
	F[I0] = np.zeros(I0.shape)
	return F

def read_diff_expr(filenameDE, filenamePR):
	diff_expr_genes = np.array([])
	RNA_names = np.array([])
	rppa_names = np.array([])

	if not np_isempty(filenamePR):
		rppa_names, rppa_pval = read_file(filenamePR)
	else:
		print('No Input file was provided for differentially expressed genes at the proteome level.')
	
	if not np_isempty(filenameDE):
		RNA_names, RNA_pval = read_file(filenameDE)
	else:
		print('No Input file was provided for differentially expressed genes at the RNA level.')
	
	# check for duplicate genes in the input files
	if (np.unique(RNA_names).shape[1] != RNA_names.shape[1]) or (np.unique(rppa_names).shape[1] != rppa_names.shape[1]):
		print('Input files for differentially expressed genes should contain only one entry per gene.')
		print('Input files for differentially expressed genes were ignored.')
		return None
	# check it the input files for differential expression are empty
	if np_isempty(rppa_names) and np_isempty(RNA_names):
		print('Only genetically aberrant genes are used for diffusion.')
		return None
	# only differential expressed genes at the RNA level are provided
	elif np_isempty(rppa_names) and not np_isempty(RNA_names):
		pval_all = multipletests(RNA_pval, alpha=0.05, method='fdr_bh').pvals_corrected
		diff_expr_genes = RNA_names[pval_all < 0.05]
	# only differential expressed genes at the protein level are provided
	elif np_isempty(RNA_names) and not np_isempty(rppa_names):
		pval_all = multipletests(rppa_pval, alpha=0.05, method='fdr_bh').pvals_corrected
		diff_expr_genes = rppa_names[pval_all < 0.05]
	else:
		names_intersection = RNA_names.intersection(rppa_names)
		names_setdiff_rna = RNA_names.difference(rppa_names)
		names_setdiff_pr = rppa_names.difference(RNA_names)

		rna_pval_names_intersection = RNA_pval[RNA_names in names_intersection]
		rppa_pval_names_intersection = rppa_pval[rppa_names in names_intersection]

		rna_pval_names_setdiff_rna = RNA_pval[RNA_names in names_setdff_rna]
		pr_pval_names_setdiff_pr = rppa_pval[rppa_names in names_setdiff_pr]

		all_names = np.concatenate((names_intersection.conj().transpose(), names_setdiff_rna.conj().transpose(), names_setdiff_pr.conj().transpose()))
		# Fishers method
		pval_all = (1-pchisq(-2*(np.log(rna_pval_names_intersection)+np.log(rppa_pval_names_intersection)),4))
		pval_all = np.concatenate((pval_all.conj().transpose(), rna_pval_names_setdiff_rna.conj().transpose(), pr_pval_names_setdiff_pr.conj().transpose()))

		pval_all = multipletests(pval_all, alpha=0.05, method='fdr_bh').pvals_corrected
		diff_expr_genes = all_names[pval_all < 0.05]
	return diff_expr_genes

def read_file(filenamePR):
	fid = open(filenamePR, 'r')
	g = np.loadtxt(fid, dtype=srt, delimiter='\t')
	fid.close()
	names = g[:,0]
	pval = g[:,1].astype(np.float64)
	return names, pval

def run():
	pass
