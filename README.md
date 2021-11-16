# MSLDA

An R package of LDA implementations, built as part of my Compass Research Project.  

The package is so named as my project is about applying topic models such as LDA to mass spectrometry (MS) data.  

Specifically, the following functions are provided:  

* `lda_original` - as in the Blei 2003 paper, but edited to force $\alpha_i > 0$ - the only implementation that does not run in parallel

* `lda_original_par` - as above but the E-step is run in parallel - should give the same results with the same seed.

* `lda_noalpha` - since the alpha update rule in LDA_original is flawed, this version treats alpha as a hyperparameter which should be tuned.

* `lda_reshaped` - this version is a further adaptation of LDA_original, which uses a count (document-term) matrix as an input rather than the document vectors - this should give the same results as LDA_originial_par.

* `lda_reshaped_noalpha` - this is the reshaped algorithm above but with alpha as a hyperparamter - this should give the same results as LDA_noalpha.

* `lda_smoothed` - as in the Hoffman 2010 paper (batch LDA section), with edited equation for $\mathcal{L}$ - this version assumes beta is a random variable.

Note there are a number of underlying functions which are not exported and not well named/documented.
