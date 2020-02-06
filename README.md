# Averaging Correlations

The correlation between two samples is a biased estimate. This repository provides functions to explore the nature of that bias and correct it.

Many years ago Fisher established that a correlation r is biased. As a consequence he discouraged that we average several correlations. However, sometimes practitioners has absoluetely no chance to avoid it. They would gather a large sample of data if they could but all they can do is getting many small samples. Then we need to average correlations. For that purpose several correction methods have been proposed.

First of all, it is highly recommended to use the sample size when averaging correlations. Do not simply add them up but use their underlying sample size as weight. To get even better results, use one of these correction methods:

* The well-known [Fisher Z transformation](https://en.wikipedia.org/w/index.php?title=Fisher_transformation&oldid=905355965) (Fisher, 1921).
* Hotelling (1953).
* Minimum variance by (Olkin & Pratt, 1958).
* Optimised minimum variance using the precise variable k (Olkin & Pratt, 1958). This is a variation of the minimum variance. Olkin & Pratt provided a formula with a parameter k = 3. This actually an approximation of an approximation. The actual k they calculated was (-7 + 9*sqrt(2))/2).
* The precise minimum variance estimate also by Olkin & Pratt (1958). This is a complex version. It is a slow function but it is an exact formula not an approximation.


Here are a few basic insights concerning the issue.

1. Basically, you can average correlations. Of course, it only makes sense when they come from the same population. Otherwise, you'd be averaging apples with oranges. Of course, it would be preferable to get a single estimate of r from one large sample, but sometimes that is not possible and you have to obtain correlation estimates from several samples or from repeated measurement on the same sample.
2. A correlation coefficient is not unbiased. That means, when you average several correlations it will not converge to the true correlation. It will underestimate the true correlation. This effect is greatest for mid-range correlations around .05.
3. There are several way to correct for that. See [Alexander (1990)](https://link.springer.com/content/pdf/10.3758/BF03334037.pdf) for a brief overview.
4. The Fisher z is the most well known correction. There are lots of others. Some of them are more precise than Fisher z. If you're new to this, Fisher z is still a good start.
5. In any way, you should not simply average the correlations but weight them using the sample size of each correlation.
6. Sample size is also important given the fact that the bias of the correlation coefficient highly depends on the underlying sample sizes (Corey, Dunlap & Burke, 1998).
7. With sample sizes larger than 50 the bias is already pretty small and only shows on the third digit. That is negligible for many applications, which implies that you could ignore any corrections if that is the case in your field (Corey, Dunlap & Burke, 1998).
8. However, still use weighted averages.
9. Keep also in mind that a correction that improves the averaged correlation r may also affect the standard deviation. So, if you want to want to apply a mathematical operation that uses the standard deviation you have to ask this forum again.
10. all the above is only true under the assumption of bivariate normality. Violating normality may change the picture. Most papers do not consider such violations.



# References

Corey, David & Dunlap, William & Burke, Michael. (1998). "Averaging Correlations: Expected Values and Bias in Combined Pearson rs and Fisher's z Transformations". Journal of General Psychology, 125, p. 245-261. 10.1080/00221309809595548. 

Fisher, R. A. (1921). "On the 'probable error' of a coefficient of correlation deduced from a small sample". Metron. 1: 3–32. [(PDF)](https://digital.library.adelaide.edu.au/dspace/bitstream/2440/15169/1/14.pdf)

Hotelling, H. (1953). "New light on the correlation coefficient and its transformations". Journal of the Royal Statistical Society, 15, p. 193-225.

Olkin, I. & Pratt, J.W. (1958). "Unbiased Estimation of CertainCorrelation Coefficients". The Annals of Mathematical Statistics,  29 (1), p. 201-211
