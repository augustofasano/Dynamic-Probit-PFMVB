# Variational inference for the smoothing distribution in dynamic probit models

This repository is associated with the article [Fasano and Rebaudo (2021)](https://arxiv.org/abs/2104.07537). Variational inference for the smoothing distribution in dynamic probit models. The key contribution of the paper is outlined below.

> We develop a variational Bayes approach, extending the partially factorized mean-field variational approximation introduced by [Fasano, Durante and Zanella (2020)](https://arxiv.org/abs/1911.06743) for the static binary probit model to the dynamic setting exploting the theoretical results in [Fasano, Rebaudo, Durante and Petrone (2021)](https://link.springer.com/article/10.1007/s11222-021-10022-w).

This repository provides codes to implement the inference methods associated with such a new result. More precisely, we provide the `R` code to perform inference on the **smoothing distribution** under dynamic probit models with three different methods:
1. **i.i.d. sampling from the exact unified skew-normal distribution**
3. **partially factorized mean-field variational Bayes (PFM-VB) approximation**
4. **mean-field variational Bayes (MF-VB) approximation**

Structure of the repository:
* the functions to implement the above methods can be found in the `R` source file [`Functions.R`](https://github.com/augustofasano/Dynamic-Probit-PFMVB/blob/main/Functions.R)
* the financial dataset analyzed in the illustration can be found in [`Financial.Rdata`](https://github.com/augustofasano/Dynamic-Probit-PFMVB/blob/main/Financial.RData), while the original entire dataset is publicly available at [Yahoo Finance](https://finance.yahoo.com/)
* the code and a tutorial to reproduce the results in the paper is available at [`Illustration.md`](https://github.com/augustofasano/Dynamic-Probit-PFMVB/blob/main/Illustration.md)
