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




train,targets = load_data(config)
# tag_data = tag_data.T
# tag_data = tag_data.drop([0],axis=1)

batch_size = config.batch_size

x_train = torch.from_numpy(train[:]).to(torch.long).to(config.device)
y_train = torch.from_numpy(targets[:]).to(torch.float).to(config.device)

model = Model(config)
model.to(config.device)
optimizer = torch.optim.Adam(model.parameters(),lr=config.learning_rate)
loss_func = nn.MSELoss()
model.train()
print('Start iteration....')

empty=torch.empty(0,config.n_vocab*config.embed).to(config.device) # final embedding
train_acc_temp=0
for step in range(config.num_epochs):
    print('step=',step+1)
    #----------------- Training ------------------
    for batch in range(0,x_train.shape[0],batch_size):
        train_pre,out_emb,selfattention = model(x_train[batch:batch+batch_size])   
        train_loss = loss_func(train_pre.view(-1), y_train[batch:batch+batch_size])

        #----------------- Check accuracy ----------------
        
        train_acc=r2_score(y_train[batch:batch+batch_size].data.cpu(),train_pre.view(-1).data.cpu())
        train_loss_print=float(train_loss.data.cpu())
        
        #----------------- Backpropagation ---------------
        optimizer.zero_grad()   
        train_loss.backward()        
        optimizer.step()        
        
        empty=torch.cat((empty,out_emb),0)


    print("train_loss:",train_loss_print)
    print("train_acc:",train_acc)
    if train_acc>train_acc_temp:
        # save embedding, label
        np.savetxt('./Results/embedding.csv',empty.detach().cpu().numpy(),fmt='%.10f',delimiter=',')
        np.savetxt('./Results/label.csv',y_train.detach().cpu().numpy(),fmt='%.10f',delimiter=',')
        
        # save model
        torch.save(model.state_dict(), './Results/model_state_dict.pt')

        # save weight
        finaltensor=torch.mm(model.fc2.weight.data,model.fc1.weight.data)
        np.savetxt('./Results/fc1weight.csv',finaltensor.detach().cpu().numpy(),fmt='%.10f',delimiter=',')


        # save correlation from self-attention 
        np.savetxt('./Results/selfattention.csv',selfattention.detach().cpu().numpy(),fmt='%.10f',delimiter=',')
        print('Output saved')
        train_acc_temp=train_acc
    empty=torch.empty(0,config.n_vocab*config.embed).to(config.device)


