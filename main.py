
from load_data import load_data
from Transformer_model import Model
from args import Config
#---------------------------------------------------
import pandas as pd
from collections import Counter
import pandas as pd
import torch.nn as nn
import torch.nn.functional as F
from torch.autograd import Variable
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
from sklearn.model_selection import train_test_split, GroupKFold, KFold
import numpy as np
import torch
from torch import autograd
import os
from tqdm import tqdm
from sklearn.metrics import r2_score
# import shap
# imports from captum library
from captum.attr import LayerConductance, LayerActivation, LayerIntegratedGradients
from captum.attr import IntegratedGradients, DeepLift, GradientShap, NoiseTunnel, FeatureAblation
config = Config()
#train:每个基因的数字表示，类似句子中的每个单词，用数字编码不同的单词
#target:每个样本的标签，类似每个句子的标签
#score_mutation:基因突变矩阵，用于做回归任务
#gene_mutation_dia_1:基因突变相关性矩阵，用于加到注意力中，相当于自己首先计算固定注意力



train,targets = load_data(config)#加载数据
# tag_data = tag_data.T
# tag_data = tag_data.drop([0],axis=1)

batch_size = config.batch_size

x_train = torch.from_numpy(train[:]).to(torch.long).to(config.device)
y_train = torch.from_numpy(targets[:]).to(torch.float).to(config.device)

model = Model(config)#调用transformer的编码器
model.to(config.device)
optimizer = torch.optim.Adam(model.parameters(),lr=config.learning_rate)
loss_func = nn.MSELoss()#回归的任务
model.train()#模型中有BN和Droupout一定要添加这个说明
print('开始迭代....')
#开始迭代
empty=torch.empty(0,config.n_vocab*config.embed).to(config.device)#用于保存最终的embedding
train_acc_temp=0
for step in range(config.num_epochs):
    print('step=',step+1)
    #-----------------训练内容------------------
    for batch in range(0,x_train.shape[0],batch_size):
        train_pre,out_emb,selfattention = model(x_train[batch:batch+batch_size])     # 喂给 model训练数据 x, 输出预测值
        train_loss = loss_func(train_pre.view(-1), y_train[batch:batch+batch_size])

        #----------- -----计算准确率----------------
        #回归模型评价指标R2_score
        # R2_score = 1，样本中预测值和真实值完全相等，没有任何误差，表示回归分析中自变量对因变量的解释越好。
        # R2_score = 0。此时分子等于分母，样本的每项预测值都等于均值。
        train_acc=r2_score(y_train[batch:batch+batch_size].data.cpu(),train_pre.view(-1).data.cpu())
        train_loss_print=float(train_loss.data.cpu())
        #-----------------反向传播更新---------------
        optimizer.zero_grad()   # 清空上一步的残余更新参数值
        train_loss.backward()         # 以训练集的误差进行反向传播, 计算参数更新值
        optimizer.step()        # 将参数更新值施加到 net 的 parameters 上
        # if step==config.num_epochs-1:
        #     #最后一个step开始保存embedding,一个batch一个batch保存
        empty=torch.cat((empty,out_emb),0)


    print("train_loss:",train_loss_print)
    print("train_acc:",train_acc)
    if train_acc>train_acc_temp:
        #保存embedding和label
        np.savetxt('/home/jby2/SNN-master/code_by_jby_for_attention/Results/embedding.csv',empty.detach().cpu().numpy(),fmt='%.10f',delimiter=',')
        np.savetxt('/home/jby2/SNN-master/code_by_jby_for_attention/Results/label.csv',y_train.detach().cpu().numpy(),fmt='%.10f',delimiter=',')
        
        #保存模型
        torch.save(model.state_dict(), '/home/jby2/SNN-master/code_by_jby_for_attention/Results/model_state_dict.pt')

        #保存权重
        finaltensor=torch.mm(model.fc2.weight.data,model.fc1.weight.data)
        np.savetxt('/home/jby2/SNN-master/code_by_jby_for_attention/Results/fc1weight.csv',finaltensor.detach().cpu().numpy(),fmt='%.10f',delimiter=',')


        #保存自注意力
        np.savetxt('/home/jby2/SNN-master/code_by_jby_for_attention/Results/selfattention.csv',selfattention.detach().cpu().numpy(),fmt='%.10f',delimiter=',')
        print('保存成功')
        train_acc_temp=train_acc
    empty=torch.empty(0,config.n_vocab*config.embed).to(config.device)


