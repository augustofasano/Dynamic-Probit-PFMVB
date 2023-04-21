Introduction
============

As described in the [`README.md`](https://github.com/augustofasano/Dynamic-Probit-PFMVB/blob/main/README.md) file, this tutorial contains general guidelines and code to **perform the comparison for the financial applications** in the paper.

Preliminary operations
======================

Once the files [`Financial.RData`](https://github.com/augustofasano/Dynamic-Probit-PFMVB/blob/main/Financial.RData) and [`Functions.R`](https://github.com/augustofasano/Dynamic-Probit-PFMVB/blob/main/Functions.R) have been downloaded, set the working directory to the folder where they are located. Then, clean the workspace and load `Financial.RData` along with the source file `Functions.R` and other useful packages.

``` r
source("Functions.R")
seed=123
set.seed(seed)
library("TruncatedNormal"); library("mvtnorm"); library("Matrix")
library("ggplot2"); library("latex2exp")

load("Financial.RData")
```

Inference on the smoothing distribution
======================

We start our analysis setting the hyperparameters of the model and computing the parameters of the smoothing distribution
````r
p = ncol(matFF)
n = length(y)

SigmaEps = diag(c(0.01,0.01))
GG       = diag(p)
mu0      = matrix(0, nrow = p, ncol = 1)
Sigma0   = diag(rep(3,p))

smoothParams = getSmooothParams(y,
                                mu0 = mu0,
                                Sigma0 = Sigma0,
                                SigmaEps = SigmaEps,
                                matFF = matFF,
                                GG =GG)
Omega   = smoothParams$Omega

indSeq1 = seq(from = 1, by = 2, length.out = n)
indSeq2 = indSeq1 + 1
````
We perform inference via **i.i.d. sampling from the exact unified skew-normal posterior**.
``` r
nSim  = 1e4
smoot = smoothing(y=y,
                  mu0 = mu0,
                  Sigma0 = Sigma0,
                  SigmaEps = SigmaEps,
                  matFF = matFF,
                  GG = GG, seed = seed, nSim = nSim)

meanTheta1 = apply(X=smoot$smoothState$values[,,1], MARGIN=2, FUN = mean)
meanTheta2 = apply(X=smoot$smoothState$values[,,2], MARGIN=2, FUN = mean)
sdTheta1   = apply(X=smoot$smoothState$values[,,1], MARGIN=2, FUN = sd)
sdTheta2   = apply(X=smoot$smoothState$values[,,2], MARGIN=2, FUN = sd)
```

We perform inference via our proposed **partially factorized mean-field variational Bayes (PFM-VB)** approximation.
````r
paramsPFM = getParamsPFM(X = X,
                         y = y,
                         Omega = Omega,
                         moments = T,
                         tolerance = 1e-3)

meanTheta1_PFM = paramsPFM$postMoments.meanBeta[indSeq1]
meanTheta2_PFM = paramsPFM$postMoments.meanBeta[indSeq2]
sdTheta1_PFM   = sqrt(paramsPFM$postMoments.varBeta[indSeq1])
sdTheta2_PFM   = sqrt(paramsPFM$postMoments.varBeta[indSeq2])
````
We perform inference via the **mean-field variational Bayes (MF-VB)** approximation.
````r
paramsMF = getParamsMF(X = X,
                        y = y,
                        Omega = Omega,
                        tolerance = 1e-3)

meanTheta1_MF = paramsMF$meanBeta[indSeq1]
meanTheta2_MF = paramsMF$meanBeta[indSeq2]
sdTheta1_MF   = sqrt(paramsMF$diagV[indSeq1])
sdTheta2_MF   = sqrt(paramsMF$diagV[indSeq2])
````

Graphical comparison of the previous methods
======================

We plot the mean and the mean +/- standard deviation of the smoothing trajectories obtained with the three aforementioned methods.
````r
Mean   = c(meanTheta1,meanTheta2,
           meanTheta1_PFM,meanTheta2_PFM,
           meanTheta1_MF,meanTheta2_MF) 
Low    = c(meanTheta1 - sdTheta1,meanTheta2 - sdTheta2,
           meanTheta1_PFM - sdTheta1_PFM, meanTheta2_PFM - sdTheta2_PFM,
           meanTheta1_MF - sdTheta1_MF, meanTheta2_MF - sdTheta2_MF)
Upp    = c(meanTheta1 + sdTheta1,meanTheta2 + sdTheta2,
           meanTheta1_PFM + sdTheta1_PFM, meanTheta2_PFM + sdTheta2_PFM,
           meanTheta1_MF + sdTheta1_MF, meanTheta2_MF + sdTheta2_MF)
Par    = as.factor(rep(c(rep(1,n),rep(2,n)),3))
Method = as.factor(rep(c("IID","PMF-VB","MF-VB"),each=2*n))
Time   = rep(1:n,6)

Data_plot             = data.frame(Mean,Low,Upp,Par,Method,Time)
levels(Data_plot$Par) = c("1"=TeX("$\\theta_{1}$"), "2"=TeX("$\\theta_{2}$"))

Plot_smooth = ggplot(Data_plot,aes(x=Time,col=Method))+
  geom_line(aes(y=Mean),size=1.2)+
  geom_line(aes(y=Low),linetype = "dashed",size=1.2)+
  geom_line(aes(y=Upp),linetype = "dashed",size=1.2)+
  scale_colour_manual(values = alpha(c("black", "red","blue"),c(1,1,1)))+ theme_bw()+
  theme(axis.title=element_blank(),legend.title = element_blank(),axis.text=element_text(size=20),strip.text = element_text(size=30),legend.text = element_text(size=30))+
  facet_wrap(~Par,nrow=2,scales = "free_y",labeller=label_parsed,strip.position = "right")
````
![alt text](https://github.com/augustofasano/Dynamic-Probit-PFMVB/blob/main/ComparisonTrajectories.png)

Finally, we show the boxplot of the differences of the mean and the log standard deviations of the smoothing parameters obtained with the **PFM-VB** solution and the **MF-VB** solution, using the inferences obtained via **i.i.d. sampling from the exact unified skew-normal posterior** as benchmark.

````r
PFM_Diff_mean_Theta1  = meanTheta1 - meanTheta1_PFM
PFM_Diff_mean_Theta2  = meanTheta2 - meanTheta2_PFM
MF_Diff_mean_Theta1   = meanTheta1 - meanTheta1_MF
MF_Diff_mean_Theta2   = meanTheta2 - meanTheta2_MF
PFM_Diff_logsd_Theta1 = log(sdTheta1) - log(sdTheta1_PFM)
PFM_Diff_logsd_Theta2 = log(sdTheta2) - log(sdTheta2_PFM)
MF_Diff_logsd_Theta1  = log(sdTheta1) - log(sdTheta1_MF)
MF_Diff_logsd_Theta2  = log(sdTheta2) - log(sdTheta2_MF)

Diff_Mean    = c(PFM_Diff_mean_Theta1, PFM_Diff_mean_Theta2, MF_Diff_mean_Theta1, MF_Diff_mean_Theta2)
Diff_logsd   = c(PFM_Diff_logsd_Theta1, PFM_Diff_logsd_Theta2, MF_Diff_logsd_Theta1, MF_Diff_logsd_Theta2)

Par          = as.factor(rep(rep(c(rep(1,n),rep(2,n)),2),2))
Method       = as.factor(rep(rep(c("PFM-VB","MF-VB"),each=2*n),2))
Functional   = as.factor(rep(1:2,each=4*n))
Diff         = c(Diff_Mean,Diff_logsd)
Data_boxplot = data.frame(Diff,Method,Par,Functional)
col_meth     = c("red","blue")
levels(Data_boxplot$Par)        = c("1"=TeX("$\\theta_{1}$"), "2"=TeX("$\\theta_{2}$"))
levels(Data_boxplot$Functional) = c("1"=TeX("$\\Delta mean$"), "2"=TeX("$\\Delta \\log(sd)$"))
Plot_Diff = ggplot(Data_boxplot, aes(y=Diff, x=Method,col=Method)) + 
  geom_boxplot()+ scale_colour_manual(values = alpha(col_meth,c(1,1)))+
  theme_bw()+ geom_hline(yintercept=0, linetype="dashed", size=2)+
  theme(axis.title=element_blank(),axis.text=element_text(size=20) , axis.text.x = element_text(colour = col_meth), strip.text = element_text(size=30),legend.position = "none")+
  facet_grid(Functional~Par,scales = "free", labeller=label_parsed)
````
![alt text](https://github.com/augustofasano/Dynamic-Probit-PFMVB/blob/main/ComparisonBoxplot.png)
