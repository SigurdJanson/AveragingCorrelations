# AveragingCorrelations

The correlation between two samples is a biased estimate. This repository provides functions to explore the nature of that bias and correct it.

Many years ago Fisher established that a correlation r is biased. As a consequence he discouraged that we average several correlations. However, sometimes practitioners has absoluetely no chance to avoid it. They would gather a large sample of data if they could but all they can do is getting many small samples. Then we need to average correlations. For that purpose several correction methods have been proposed.

First of all, it is highly recommended to use the sample size when averaging correlations. Do not simply add them up but use their underlying sample size as weight. To get even better results, use one of these correction methods:

* The well-known Fisher Z transformation (Fisher, 1921).
* Hotelling (...).
* Minimum variance by (Olkin & Pratt, 1958).
* Optimised minimum variance using the precise variable k (Olkin & Pratt, 1958). This is a variation of the minimum variance. Olkin & Pratt provided a formula with a parameter $k = 3$. This actually an approximation of an approximation. The actual k they calculated was (-7 + 9*sqrt(2))/2).
* The precise minimum variance estimate also by Olkin & Pratt (1958). This is a complex version. It is a slow function but it is the exact derivation not an approximation.

