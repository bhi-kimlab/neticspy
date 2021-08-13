import time
import os
import numpy as np

from .const import DIRECTION
from . import util

def insulated_diff(W, b):
    #start = time.time()
    temp = np.identity(len(W)) - (1-b) * W
    F = b * np.linalg.pinv(temp)
    #print(f'Insulated diffusion was done in {util.time_format(time.time()-start)}.')
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

def diffuse_all_permutation(diffusion_direction, aberrant_gene_seeds, deg_seeds, F, F_opp):
    if diffusion_direction in [DIRECTION.DOWN, DIRECTION.BOTH]:
        aberration_scores = diffuse_many(aberrant_gene_seeds, F)

    # Diffuse DE labels
    if diffusion_direction in [DIRECTION.UP, DIRECTION.BOTH]:
        de_scores = diffuse_many(deg_seeds, F_opp)

    # Combines scores by multiplication (use both mutations and DE)
    if diffusion_direction == DIRECTION.BOTH:
        ret_scores_all = aberration_scores * de_scores
    elif diffusion_direction == DIRECTION.DOWN:
        ret_scores_all = aberration_scores
    else:
        ret_scores_all = de_scores

    return ret_scores_all

def diffuse(seed_idx, F):
    scores = F[seed_idx, :].mean(axis=0)
    return scores

def diffuse_many(seed_one_hot, F):
    """
    seed_one_hot: Seeds in one-hot representation. size (num_permutation x num_seeds)
    scores: size (num_permutation x num_gene)
    """
    num_seeds = seed_one_hot[0].sum()
    scores = 1 / num_seeds * np.matmul(seed_one_hot, F)

    return scores
