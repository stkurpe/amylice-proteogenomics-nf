  ##AMYPred-FRL

AMYPred-FRL is a novel approach to accurately predict amyloid proteins using feature representation learning. In this research study, we combined six well-known ML algorithms with ten different sequence descriptors to generate generate 60 probabilistic features (PFs), as opposed to state-of-the-art methods developed by the single feature-based approach.  Logistic Regression Recursive Characteristic Elimination (LR-RFE) method is used to find the optimal number of m feature from 60 PF to improve prediction efficiency. Finally, the 20 selected PFs were feed into the random forest method using a meta-predictor approach to construct the final model

###AMYPred-FRL uses the following dependencies:
Installation

Download AMYPred-FRL by

git clone https://github.com/saeed344/AMYPred-FRL

Installation has been tested in OS win 10 with Python 3.8.3

Since the package is written in python 3.8.3, python 3.8.3 with the pip tool must be installed first.
 AMYPred-FRL uses the following dependencies: numpy, scipy, scikit-learn, pandas, Xgboost 

You can install these packages first, by the following commands:

pip install numpy

pip install scipy

pip install scikit-learn

pip install pandas

pip install Xgboost

File Description
1. data
The ‘data’ file contains all datasets used in this work, including training set and independent set.
3. model
The ‘model’ file contains save model.
4. Blind_test.rar
The ‘Blind_test’ file contains the code for getting Blind dataset results.
5. model_test
The ‘model_test’ file contains code for evaluation the prediction ability of the model on existing indpendent data.

###Guiding principles: 

**The dataset file contains  TR_P_132.fasta, TR_N_305.fasta, TS_P_33.fasta, TS_N_77.fasta.

script test
#####################
To check whether the project can work normally, we can run

Stand_alone_AmyPredFRL using pycharm

To check blind dataset sequnces try Blind_test.py

 



