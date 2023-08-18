import torch
class Config(object):

    def __init__(self):
        self.model_name = 'Transformer'
        self.embedding_pretrained = None                                 
        self.device = torch.device('cuda:3' if torch.cuda.is_available() else 'cpu')
        self.dropout = 0                                              
        self.num_classes = 1                        
        self.num_epochs = 300                                           
        self.batch_size = 64                                           
        self.pad_size = 256                                             
        self.n_vocab = 487
        self.learning_rate = 5e-04                                       
        self.embed = 1           
        self.dim_model = 1
        self.hidden = 8
        self.last_hidden = 512
        self.num_head = 1
        self.num_encoder = 1 
        self.n_splits = 5
        self.score_mutation_path = "./data/score_mutation.npy"
        self.start = 0 
    