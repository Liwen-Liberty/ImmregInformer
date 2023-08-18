import pandas as pd
import numpy as np
import torch
from tqdm import tqdm
#-------------------------- Data load ----------------------------
def load_data(config):
    score_gene_num_vector=np.load("./data/score_gene_num_vector.npy")

    train = np.array(score_gene_num_vector[:,1:],dtype=int)
    targets = np.array(score_gene_num_vector[:,0],dtype=float)

    return train,targets

