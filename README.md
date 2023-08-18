# ImmRegInformer: an interpretable Transformer-based method for prioritizing immune-regulatory cancer drivers
===========================================================================


[![license](https://img.shields.io/badge/python_-3.8.0_-blue)](https://www.python.org/)
[![license](https://img.shields.io/badge/torch_-1.12.0_-blue)](https://pytorch.org/)
[![license](https://img.shields.io/badge/scanpy_-1.9.0_-blue)](https://scanpy.readthedocs.io/en/stable/)
[![license](https://img.shields.io/badge/anndata_-0.8.0_-blue)](https://anndata-tutorials.readthedocs.io/en/latest/index.html/)
[![license](https://img.shields.io/badge/R_-4.2.2_-blue)](https://www.r-project.org/)

The identification and prioritization of immune-regulatory cancer driver mutations present a promising study for precision immunotherapy of cancer but remain considerable challenges. As there is no existing proven and efficient method to systematically explore the regulatory relationship between cancer driver mutations and immune response, we introduced a novel method, ImmRegInformer, that leverages the powerful self-attention mechanisms of the Transformer model and the lasso-regularised ordinal regression to address these challenges. In particular, our approach integrated the mutation co-occurrence information with the correlation value generated by the self-attention module to discern the underlying relationships between different driver mutations when regulating the immune cytolytic activity (CYT). Using ImmRegInformer, we identified 250 immune-regulating driver mutations out of 487 candidates in 8223 pan-cancer samples. The regulatory driver mutations had significantly higher mutation frequency and more cohesive interactions with the cytolytic signature genes, suggesting the accuracy of our identification. Based on the correlation analysis of the correlation value generated by the self-attention module and mutation co-occurrence odd ratio, we found their complementary roles in exploring the synergistic immune-regulatory mutation pairs. For model interpretability, the CYT-regulatory driver mutations were prioritized based on the integrated weight values. We verified the top-ranked driver mutations exhibited dominant associations with the CYT signature genes. In conclusion, this study paves the way for an advanced, interpretable machine learning framework that has potential implications in dissecting the immune regulatory factors in cancer. It also underscores the importance of employing deep learning methods like Transformer to unlock hidden insights into the biological complexities of cancer immunity, offering a new avenue for research on the immune regulatory mechanism and potential clinical applications. ImmRegInformer is freely available at https://github.com/Liwen-Liberty/ImmregInformer.

![Image text](https://github.com/Liwen-Liberty/ImmregInformer/blob/main/Figures/Figure1.png)


To fetch indicators of immune response regulation from tumor somatic mutation profile, we proposed a novel method called ImmRegInformer. The workflow was illustrated in Figure 1 for better understanding. At first, the mutation annotation format (MAF) file of somatic mutation was transformed into a binary matrix that indicated the mutation status of all candidate cancer driver genes. This matrix served as input data whose rows corresponded to cancer samples (observations) and columns to the driver mutation status (variables). Each driver mutation variable was then randomly initialized as an embedding vector, and fed to transformer encoder lay-ers. Then, to capture the mutation co-occurrence associations among cancer driver genes, the odd ratio matrix was integrated into the multi-head self-attention module. The dot product of the original mutation status profile and the embeddings generated by the ImmRegInformer encoder was calculated to concentrate the information from mutated genes. Subsequently, the modified embeddings encapsulate information regarding the interactions between various driver mutations was served as the final feature vectors. The final feature vector was further utilized in predicting the CYT score through a multi-layer perceptron (MLP) module. The weights in the MLP can be used to interpret the importance of the variables. Since there are challenges in defining an appropriate cut-off for variable weights, the lasso-regularized ordinal re-regression analysis was performed to discern driver mutations with significant immunoregulatory effects. Taking the CYT score as the dependent variable, we utilized the modified embeddings obtained from the Transformer encoder model as the independent variables.
## Table of Contents

- [Installation](#installation)
- [Quick start](#quick-start)
- [Contributing](#contributing)
- [Cite](#cite)
- [Contacts](#contacts)
- [License](#license)


## Installation

ImmRegInformer is tested to work under the:

```
* Python 3.8.0
* Torch 1.12.0
* Scanpy 1.9.0
* Anndata 0.8.0
* R 4.2.2
* Numpy 1.23.5
* Other basic Python and r toolkits
```
### Installation of other dependencies
* Install [R package glmnetcr](https://github.com/cran/glmnetcr) using ` devtools::install_github("cran/glmnetcr") ` in the R environment if you encounter any issue.


# Quick start
To reproduce our results:

## Data Preparation
- Calculation of the cytolytic activity score of tumor samples
- Calculation of the mutation co-occurrence odd ratios and p-values for each driver mutation pair


## 1. Application and validation in the identification of immune-regulators
```
python ./main.py
```

**Arguments**:

| **Arguments** | **Detail** |
| --- | --- |
| **Solid_tumor_mutation_status_mat.csv** | Mutation status profile of 8223 solid tumor samples with matched expression data. |
| **Solid_tumor_CYT_score_df.csv** | Cytolytic activity score of 8223 solid tumor samples. |
| **Solid_tumor_RNA_data.RData** | RNA-seq data of 8223 solid tumor samples. |
| **Solid_tumor_mutation_cooccurrence_mat.csv** | Co-occurrence analysis of each driver mutation pair based on Fisher's exact test. |

**Values**:

| **Output** | **Detail** |
| --- | --- |
| **model_state_dict.pt** | ImmRegInformer model. |
| **selfattention.csv** | Correlation matrix generated by the attention network in ImmregInformer. |
| **embedding.csv** | Embeddings generated by ImmregInformer. |
| **fc1weight.csv** | Variable  weights in the MLP. |




<div align="center">
  <img src="https://github.com/Liwen-Liberty/ImmregInformer/blob/main/Figures/Figure2.png" alt="Editor" width="500">
</div>

===========================================================================

<div align="center">
  <img src="https://github.com/Liwen-Liberty/ImmregInformer/blob/main/Figures/Figure3.png" alt="Editor" width="500">
</div>

## 2. Interplay of cancer drivers in regulating immune response

See Rscripts


Visualization of results:
<div align="center">
  <img src="https://github.com/Liwen-Liberty/ImmregInformer/blob/main/readmepic/Figure4.png" alt="Editor" width="500">
</div>



## 3. Prioritizing cancer drivers in order of contribution to immune regulation

See Rscripts

Visualization of results:


===========================================================================





# Contributing

liwen Xu, boya Ji, yijun Peng, wending Pi, shaoliang Peng*

# Cite
<p align="center">
  <a href="https://clustrmaps.com/site/1bpq2">
     <img width="200"  src="https://clustrmaps.com/map_v2.png?cl=ffffff&w=268&t=m&d=4hIDPHzBcvyZcFn8iDMpEM-PyYTzzqGtngzRP7_HkNs" />
   </a>
</p>


# Contacts
If you have any questions or comments, please feel free to email: xuliwen@hnu.edu.cn.

# License

[MIT ? Richard McRichface.](../LICENSE)
