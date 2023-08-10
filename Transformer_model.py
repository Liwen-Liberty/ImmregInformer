
from os import device_encoding
from random import random
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
import copy
import random
from args import Config
# from main import train_idx,expain_num
config = Config()

'''Attention Is All You Need'''


class Model(nn.Module):
    def __init__(self, config):
        super(Model, self).__init__()
        if config.embedding_pretrained is not None:
            self.embedding = nn.Embedding.from_pretrained(config.embedding_pretrained, freeze=False)
        else:
            # self.embedding = nn.Embedding(config.n_vocab, config.embed, padding_idx=config.n_vocab - 1)
            self.embedding = nn.Embedding(config.n_vocab, config.embed)
        self.postion_embedding = Positional_Encoding(config.embed, config.pad_size, config.dropout, config.device)
        self.encoder = Encoder(config.dim_model, config.num_head, config.hidden, config.dropout)
        self.encoders = nn.ModuleList([
            copy.deepcopy(self.encoder)
            # Encoder(config.dim_model, config.num_head, config.hidden, config.dropout)
            for _ in range(config.num_encoder)])

        # self.fc1 = nn.Linear(config.pad_size * config.dim_model, config.num_classes)
        # self.fc2 = nn.Linear(config.last_hidden, config.num_classes)
        self.my_score_mutation = My_score_mutation(config.score_mutation_path, config.batch_size,config.n_vocab,config.embed,config.device,config.start)
        
        self.fc1 = nn.Linear(config.n_vocab * config.dim_model, 256)
        self.relu1 = torch.nn.ReLU()
        self.droupout1 = torch.nn.Dropout(p=0)
        self.fc2 = nn.Linear(256, config.num_classes)
    def forward(self, x):
        out = self.embedding(x)
        #return out
        out = self.postion_embedding(out)
        for encoder in self.encoders:
            out,selfattention = encoder(out)#64,487,10
        out = out.view(out.size(0), -1)#64,4870
        # out = torch.mean(out, 1)
        out_emb = self.my_score_mutation(out)#对out进行进步处理，按照基因的突变信息进行选择
        out = self.fc1(out_emb)
        out = self.relu1(out)
        out = self.droupout1(out)
        out = self.fc2(out)
        return out,out_emb,selfattention
    
    # def forward(self,x):
    #     x = x.long()
    #     out = self.embedding(x)
    #     #return out
    #     out = self.postion_embedding(out)
    #     for encoder in self.encoders:
    #         out = encoder(out)#64,487,10
    #     out = out.view(out.size(0), -1)#64,4870
    #     # out = torch.mean(out, 1)
    #     index = train_idx
    #     out = self.my_score_mutation(out,index,expain_num,expain_num+1)#对out进行进步处理，按照基因的突变信息进行选择
    #     out = self.fc1(out)
    #     out = self.relu1(out)
    #     out = self.droupout1(out)
    #     out = self.fc2(out)
    #     out = out.double()
    #     return out


class My_score_mutation(nn.Module):
    #我们自己定义的类，目的是实现在最后分类之前，加入突变矩阵的信息
    def __init__(self, score_mutation_path, batch_size,n_vocab,embed,device,start):
        super(My_score_mutation, self).__init__()
        self.score_mutation_path = score_mutation_path
        self.batch_size = batch_size
        self.n_vocab = n_vocab
        self.embed = embed
        self.device = device
        self.start = start
    def forward(self, x):
        #x:batch_size,n_vocab*embed
        score_mutation=np.load(self.score_mutation_path)#突变数组
        #构造一个矩阵，如果有突变则保留该gene embedding，否则让它对应乘以0，不考虑该基因
        start=self.start
        final=[]
        for i in range(start,start+x.shape[0]):
            tem_score_mutation=score_mutation[i,1:]
            pair=[]
            for j in range(len(tem_score_mutation)):
                if tem_score_mutation[j]==0:
                    pair.extend([0]*self.embed)
                if tem_score_mutation[j]==1:
                    pair.extend([1]*self.embed)
            final.append(pair)
        final_tensor=torch.LongTensor(final).to(self.device)
        # final_tensor=torch.ones(final_tensor.shape).to(self.device)
        out = final_tensor.mul(x).to(self.device) #点积
        # out = (final_tensor+x).to(self.device) #相加
        if start+x.shape[0]==score_mutation.shape[0]:
            self.start=0
        else:
            self.start=start+x.shape[0]
        return out




class Encoder(nn.Module):
    def __init__(self, dim_model, num_head, hidden, dropout):
        super(Encoder, self).__init__()
        self.attention = Multi_Head_Attention(dim_model, num_head, dropout)
        self.feed_forward = Position_wise_Feed_Forward(dim_model, hidden, dropout)

    def forward(self, x):
        out,selfattention = self.attention(x)
        out = self.feed_forward(out)
        return out,selfattention


class Positional_Encoding(nn.Module):
    def __init__(self, embed, pad_size, dropout, device):
        super(Positional_Encoding, self).__init__()
        self.device = device
        # self.pe = torch.tensor([[pos / (10000.0 ** (i // 2 * 2.0 / embed)) for i in range(embed)] for pos in range(pad_size)])
        # self.pe[:, 0::2] = np.sin(self.pe[:, 0::2])
        # self.pe[:, 1::2] = np.cos(self.pe[:, 1::2])
        # self.dropout = nn.Dropout(dropout)

    def forward(self, x):
        # out = x + nn.Parameter(self.pe, requires_grad=False).to(self.device)
        # out = self.dropout(out)
        out = x.to(self.device) #这里没有位置编码
        return out


class Scaled_Dot_Product_Attention(nn.Module):
    '''Scaled Dot-Product Attention '''
    def __init__(self):
        super(Scaled_Dot_Product_Attention, self).__init__()
    def forward(self, Q, K, V, scale=None):
        '''
        Args:
            Q: [batch_size, len_Q, dim_Q]
            K: [batch_size, len_K, dim_K]
            V: [batch_size, len_V, dim_V]
            scale: 缩放因子 论文为根号dim_K
        Return:
            self-attention后的张量，以及attention张量
        '''
        attention = torch.matmul(Q, K.permute(0, 2, 1))

        add_our_attention=True
        if add_our_attention:
            gene_mutation_dia_1=np.load("/home/jby2/SNN-master/code_by_jby/data/gene_mutation_dia_0.npy")
            tem_tensor=torch.FloatTensor(gene_mutation_dia_1)
            tem_tensor=tem_tensor.unsqueeze(0).repeat(attention.shape[0], 1, 1).to(config.device)#因为是多头注意力，所以有多组需要加上自己定义的attention
            attention=attention+tem_tensor
        if scale:
            attention = attention * scale
        
        # if mask:  # TODO change this
        #     attention = attention.masked_fill_(mask == 0, -1e9)
        attention = F.softmax(attention, dim=-1)
        
        context = torch.matmul(attention, V)
        return context,attention[0]


class Multi_Head_Attention(nn.Module):
    def __init__(self, dim_model, num_head, dropout=0.0):
        super(Multi_Head_Attention, self).__init__()
        self.num_head = num_head
        assert dim_model % num_head == 0
        self.dim_head = dim_model // self.num_head
        self.fc_Q = nn.Linear(dim_model, num_head * self.dim_head)
        self.fc_K = nn.Linear(dim_model, num_head * self.dim_head)
        self.fc_V = nn.Linear(dim_model, num_head * self.dim_head)
        self.attention = Scaled_Dot_Product_Attention()
        self.fc = nn.Linear(num_head * self.dim_head, dim_model)
        self.dropout = nn.Dropout(dropout)
        self.layer_norm = nn.LayerNorm(dim_model)

    def forward(self, x):
        batch_size = x.size(0)
        Q = self.fc_Q(x)
        K = self.fc_K(x)
        V = self.fc_V(x)
        Q = Q.view(batch_size * self.num_head, -1, self.dim_head)
        K = K.view(batch_size * self.num_head, -1, self.dim_head)
        V = V.view(batch_size * self.num_head, -1, self.dim_head)
        # if mask:  # TODO
        #     mask = mask.repeat(self.num_head, 1, 1)  # TODO change this
        scale = K.size(-1) ** -0.5  # 缩放因子
        context,selfattention = self.attention(Q, K, V, scale)

        context = context.view(batch_size, -1, self.dim_head * self.num_head)
        out = self.fc(context)
        out = self.dropout(out)
        out = out + x  # 残差连接
        out = self.layer_norm(out)
        return out,selfattention


class Position_wise_Feed_Forward(nn.Module):
    def __init__(self, dim_model, hidden, dropout=0.0):
        super(Position_wise_Feed_Forward, self).__init__()
        self.fc1 = nn.Linear(dim_model, hidden)
        self.fc2 = nn.Linear(hidden, dim_model)
        self.dropout = nn.Dropout(dropout)
        self.layer_norm = nn.LayerNorm(dim_model)

    def forward(self, x):
        out = self.fc1(x)
        out = F.relu(out)
        out = self.fc2(out)
        out = self.dropout(out)
        out = out + x  # 残差连接
        out = self.layer_norm(out)
        return out
