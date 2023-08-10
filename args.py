import torch
class Config(object):

    """配置参数"""
    def __init__(self):
        self.model_name = 'Transformer'
        self.embedding_pretrained = None                                 # 预训练词向量
        self.device = torch.device('cuda:2' if torch.cuda.is_available() else 'cpu')   # 设备
        # self.device = torch.device('cpu')   # 设备
        self.dropout = 0                                              # 随机失活
        self.num_classes = 1                        # 类别数（回归任务最后输出一个得分即可）
        self.num_epochs = 300                                            # epoch数
        self.batch_size = 64                                           # mini-batch大小
        self.pad_size = 256                                              # 每句话处理成的长度(短填长切)（用不到）
        self.n_vocab = 487#词的个数，我们这里就是基因的个数487
        self.learning_rate = 5e-04                                       # 学习率
        self.embed = 1           # 词向量维度
        self.dim_model = 1
        self.hidden = 8
        self.last_hidden = 512#（用不到）
        self.num_head = 1#多头个数
        self.num_encoder = 1 #编码层个数
        self.n_splits = 5#k折交叉验证
        self.score_mutation_path="/home/jby2/SNN-master/code_by_jby/data/score_mutation.npy"
        self.start=0 #这里新加一个参数，用于score_mutation的batch的计数
    
