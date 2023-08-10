import pandas as pd
import numpy as np
import torch
from tqdm import tqdm
#--------------------------加载数据----------------------------
def load_data(config):
    score_gene_num_vector=np.load("/home/jby2/SNN-master/code_by_jby/data/score_gene_num_vector.npy")

    # score_mutation=np.load("/home/jby2/SNN-master/code_by_jby/data/score_mutation.npy")

    train = np.array(score_gene_num_vector[:,1:],dtype=int)
    targets = np.array(score_gene_num_vector[:,0],dtype=float)
    # data_target =  pd.read_csv(r"/home/jby2/SNN-master/data/Solid_tumor_mutation_status_mat.csv",header=None)

    # gene_mutation_dia_1=np.load("/home/jby2/SNN-master/code_by_jby/data/gene_mutation_dia_1.npy")
    return train,targets

