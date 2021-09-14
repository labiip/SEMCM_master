# SEMCM_master
Predicting anti-cancer drug response through a self-expressive matrix completion model

Developer: Lin,Zhang(lin.zhang@cumt.edu.cn) from Engineering Research Center of Intelligent Control for Underground Space, School of Information and Control Engineering, China University of Mining and Technology

## **Requirement**

4GB memory

MATLAB R2015a or later

## **Related data information need to first load in SEMCM method** 

- /data/CCLE.mat
- /data/GDSC.mat

The first file CCLE.mat is a matrix of known drug sensitivity relationships between cell lines and drugs from Cancer Cell Line Encyclopedia. 
The second file GDSC.mat is a matrix of known drug sensitivity relationships between cell lines and drugs from Genomics of Drug Sensitivity in Cancer.

## **Run SEMCM to predict associations between cell line and drug**
*****************************************************************************
We provide two drug sensitivity relationships matrices.
To demonstrate SEMCM can achieve a good performance with better ave_PCC and ave_RMSE, we performing 10-fold cross validation to evaluate our method on predicting known cell line-drug associations by running the following code. 
(This process may need to take a few hours.)
	
	run SEMCM_CCLE.m
	run SEMCM_GDSC.m
*****************************************************************************

## **Contact**

Please feel free to contact us if you need any help: lin.zhang@cumt.edu.cn
