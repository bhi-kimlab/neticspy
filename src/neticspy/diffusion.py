import time
import numpy as np

from .const import DIRECTION

def insulated_diff(W, b):
	start = time.time()
	temp = np.identity(len(W)) - (1-b) * W
	F = b * np.linalg.pinv(temp)
	print(f'Insulated diffusion was done in {time.time() - start:.2f} seconds.')
	return F

def diffuse_all(diffusion_direction, aberrant_genes, degs, network_genes, F, F_opp):
	# Diffuse mutation labels
	if diffusion_direction == DIRECTION.DOWN or diffusion_direction == DIRECTION.BOTH:
		aberration_scores = diffuse(aberrant_genes, network_genes, F)

	# Diffuse DE labels
	if diffusion_direction == DIRECTION.UP or diffusion_direction == DIRECTION.BOTH:
		de_scores = diffuse(degs, network_genes, F_opp)

	# Combines scores by multiplication (use both mutations and DE)
	if diffusion_direction == DIRECTION.BOTH:
		ret_scores_all = aberration_scores * de_scores
	elif diffusion_direction == DIRECTION.DOWN:
		ret_scores_all = aberration_scores
	else:
		ret_scores_all = de_scores

	return ret_scores_all

def diffuse(target_genes, network_genes, F):
	positions = np.isin(network_genes, target_genes).nonzero()[0]
	scores = (1 / len(positions)) * np.matmul(np.ones((1, len(positions))), F[positions, :])
	return scores
