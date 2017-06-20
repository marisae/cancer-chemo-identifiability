# cancer-chemo-identifiability
Code for structural and practical identifiability of two cancer chemotherapy models. This is the source code for:

A Confidence Building Exercise in Data and Identifiability: Modeling Cancer Chemotherapy as a Case Study. Marisa C. Eisenberg and Harsh V. Jain (equal coauthors). In revision. 

Questions? Contact Harsh Jain (hjain@fsu.edu) or Marisa Eisenberg (marisae@umich.edu).

The code given here will simulate the taxol compartmental chemotherapy model described in the above manuscript (and the oxaliplatin model with minor modifications as well), and fit it to in vitro experimental data. The control data fitting code can be found in the `Control Data Code` folder. Code to run profile likelihoods to examine structural and practical identifiability is given in the `Profile Likelihood` folder (currently set up to profile 𝛼<sub>0</sub>, but can be used for any of the parameters). The current code is a bit rough, but we plan to clean it up once the paper is out. To run the code, start with `TaxolTx_fit.m`

Please cite the paper referenced above if you use our code or adapt it for your project.
