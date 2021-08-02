LDA_orginial - contains functions to run LDA as described in the Blei 2003 paper. This is also the only version which does not run the E steps in parallel.

LDA_orginial_par - same as above but the E-step is run in parallel

LDA_noalpha - since the alpha update rule in LDA_orginial is flawed, this version just treats alpha as a hyperparamter to be tuned

LDA_reshaped - this version is a further adaptation of LDA_orginial, which uses a count (document-term) matrix as an input rather than the document vectors

LDA_smoothed - this version assumes beta is a random variable, and is very close to the batch algorithm described in Hoffman 2010 
