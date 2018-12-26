/*
* Created on Fri Sep 21 10:29:50 2018
* @author: qwzhou
* @email: qiangwei.zhou2013@gmail.com
*/
#include <string>
#include <string.h>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <iterator>
#include <gsl/gsl_cdf.h>
#include <tr1/cmath>
#include <gsl/gsl_sf_gamma.h>
#include <stdio.h>
#include <stdlib.h>

//===========================================================
//fishers_exact(int expReads1, int expReads2, int obsReads1, int obsReads2)
//Calculates significance of read counts between two samples
/*Parameters:
expReads1 - Reads supporting allele 1 (expected)
expReads2 - Reads supporting allele 2 (expected)
obsReads1 - Reads supporting allele 1 (observed)
obsReads2 - Reads supporting allele 2 (observed)
Returns:
p-value P-value from Fisher's Exact Test
*/

static inline double log_sum_log(const double p, const double q) {
  if (p == 0) {return q;}
  else if (q == 0) {return p;}
  const double larger = (p > q) ? p : q;
  const double smaller = (p > q) ? q : p;
  return larger + log(1.0 + exp(smaller - larger));
};

// p(k) =  C(n1, k) C(n2, t - k) / C(n1 + n2, t)
static double log_hyper_g(const size_t k, const size_t n1, const size_t n2, const size_t t) {
  return (gsl_sf_lnfact(n1) - gsl_sf_lnfact(k) - gsl_sf_lnfact(n1 - k) +
          gsl_sf_lnfact(n2) - gsl_sf_lnfact(t - k) - gsl_sf_lnfact(n2 - (t - k)) -
          (gsl_sf_lnfact(n1 + n2) - gsl_sf_lnfact(t) - gsl_sf_lnfact(n1 + n2 - t)));
};

//fishers_exact(meth_factor, coverage_factor, meth_rest, coverage_rest)
//static 
double fishers_exact(size_t a, size_t b, size_t c, size_t d) {
  const size_t m = a + c; // sum of first column
  const size_t n = b + d; // sum of second column
  const size_t k = a + b; // sum of first row
  const double observed = log_hyper_g(a, m, n, k);
  double p = 0.0;
  for (size_t i = (n > k ? 0ul : k - n); i <= std::min(k, m); ++i) {
    const double curr = log_hyper_g(i, m, n, k);
    if (curr <= observed)
      p = log_sum_log(p, curr);
  }
  return exp(p);
};



