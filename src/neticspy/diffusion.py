import time
import numpy as np

from .const import DIRECTION

def insulated_diff(W, b):
	start = time.time()
	temp = np.identity(len(W)) - (1-b) * W
	F = b * np.linalg.pinv(temp)
	print(f'Insulated diffusion was done in {time.time() - start:.2f} seconds.')
	return F

def diffuse_all(diffusion_direction, aberrant_gene_idx, deg_idx, F, F_opp):
	# Diffuse mutation labels
	if diffusion_direction == DIRECTION.DOWN or diffusion_direction == DIRECTION.BOTH:
		aberration_scores = diffuse(aberrant_gene_idx, F)

	# Diffuse DE labels
	if diffusion_direction == DIRECTION.UP or diffusion_direction == DIRECTION.BOTH:
		de_scores = diffuse(deg_idx, F_opp)

	# Combines scores by multiplication (use both mutations and DE)
	if diffusion_direction == DIRECTION.BOTH:
		ret_scores_all = aberration_scores * de_scores
	elif diffusion_direction == DIRECTION.DOWN:
		ret_scores_all = aberration_scores
	else:
		ret_scores_all = de_scores

	return ret_scores_all

def diffuse(seed_idx, F):
	scores = 1 / len(seed_idx) * np.matmul(
		np.ones((1, len(seed_idx))),
		F[seed_idx, :]
	)
	return scores
