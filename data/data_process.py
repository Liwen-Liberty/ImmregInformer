import pandas as pd
import numpy as np
def load_data_save():
    df = pd.read_csv(r"/home/jby2/SNN-master/data/Solid_tumor_CYT_score_df.csv")
    rf = pd.read_csv(r"/home/jby2/SNN-master/data/Solid_tumor_mutation_status_mat.csv")
    #两个pands匹配构建最终的数据
    final=[]#保存最终特征矩阵
    for i in range(len(df)):
        pair=[]
        PatientID=df.loc[i]["PatientID"]
        CYT_score=df.loc[i]["CYT_score"]
        PatientID_feature=rf.loc[:,PatientID].tolist() #对应PatientID的特征列
        pair.append(CYT_score)
        pair.extend(PatientID_feature)
        final.append(pair)
    #保存到本地
    final_array=np.array(final)
    np.save("/home/jby2/SNN-master/code_by_jby/data/score_mutation.npy",final_array)

def load_data_save2():
    df = pd.read_csv(r"/home/jby2/SNN-master/data/Solid_tumor_CYT_score_df.csv")
    rf = pd.read_csv(r"/home/jby2/SNN-master/data/Solid_tumor_mutation_status_mat.csv")
    #基因的个数
    n_gene=487
    final=[]
    for i in range(len(df)):
        pair=[]
        CYT_score=df.loc[i]["CYT_score"]
        pair.append(CYT_score)
        for j in range(n_gene):
            pair.append(j)
        final.append(pair)
    
    #保存到本地
    final_array=np.array(final)
    np.save("/home/jby2/SNN-master/code_by_jby/data/score_gene_num_vector.npy",final_array)

def load_data_save3():
    df = pd.read_csv(r"/home/jby2/SNN-master/data/Solid_tumor_mutation_cooccurrence_mat.csv",sep=",",index_col=0)
    tem=df.values
    #对角线元素为1
    np.save("/home/jby2/SNN-master/code_by_jby/data/gene_mutation_dia_1.npy",tem)
    #对角线元素为0，用于看哪种效果好
    row, col = np.diag_indices_from(tem)
    tem[row,col]=0

    np.save("/home/jby2/SNN-master/code_by_jby/data/gene_mutation_dia_0.npy",tem)
if __name__ == "__main__":
    
    # 我们需要的第一个文件是，列是label feature是01向量，表示基因是否突变
    # load_data_save()
    # 我们需要的第二个文件是，列是label（score） feature是gene的数字编码，1代表基因1,2代表基因2，类似一句话中的每个词，因为这里gene用的都是一样的
    # 所以这里的特征向量都是从1到基因个数的数值
    # load_data_save2()

    #我们需要的第三个文件，是gene与gene之间的突变相关性

    load_data_save3()
    
