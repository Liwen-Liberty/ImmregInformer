import pandas as pd
import numpy as np
def load_data_save():
    df = pd.read_csv(r"./data/Solid_tumor_CYT_score_df.csv")
    rf = pd.read_csv(r"./data/Solid_tumor_mutation_status_mat.csv")
    
    final=[]
    for i in range(len(df)):
        pair=[]
        PatientID=df.loc[i]["PatientID"]
        CYT_score=df.loc[i]["CYT_score"]
        PatientID_feature=rf.loc[:,PatientID].tolist() 
        pair.append(CYT_score)
        pair.extend(PatientID_feature)
        final.append(pair)
    
    final_array=np.array(final)
    np.save("./code_by_jby/data/score_mutation.npy",final_array)

def load_data_save2():
    df = pd.read_csv(r"./data/Solid_tumor_CYT_score_df.csv")
    rf = pd.read_csv(r"./data/Solid_tumor_mutation_status_mat.csv")
    
    n_gene=487
    final=[]
    for i in range(len(df)):
        pair=[]
        CYT_score=df.loc[i]["CYT_score"]
        pair.append(CYT_score)
        for j in range(n_gene):
            pair.append(j)
        final.append(pair)
    
    
    final_array=np.array(final)
    np.save("./data/score_gene_num_vector.npy",final_array)

def load_data_save3():
    df = pd.read_csv(r"./data/Solid_tumor_mutation_cooccurrence_mat.csv",sep=",",index_col=0)
    tem=df.values
    
    np.save("./data/gene_mutation_dia_1.npy",tem)
    
    row, col = np.diag_indices_from(tem)
    tem[row,col]=0

    np.save("./data/gene_mutation_dia_0.npy",tem)
if __name__ == "__main__":
    
    load_data_save3()
    
