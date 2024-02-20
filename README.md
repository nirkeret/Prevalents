# Prevalents

The "Prevalents" file contains RCPP functions that should be loaded to the R session.

"simulation data generation" - an R file that contains the data generation described in the paper. Comment/uncomment the relevant pieces of code, according to the required setting A/B/C.

"Estimation Procedure" - an R file that contains the pairwise point-estimation, as well as the three bootstrap procedures for variance estimation. The procedure requires that the data variables be named in a certain way, as described in the beginning of the code. Also, the data should be randomly ordered.

Note that the "Prevalents" file contains a function "getZetaTermsNoC" which should be used if censoring is assumed fixed and not random. When censoring is fixed, one should also comment out the parts in the code that estimate the censoring distribution.

In the model for death-after-disease, the disease age can be used as covariate, which is how the estimation procedure is currently written. One could in theory choose to insert other functional forms of the disease age, however this option is not yet implemented.
