#ifndef LONGCALLD_MATH_UTILS_H
#define LONGCALLD_MATH_UTILS_H

int log_likelilood(double prob);
void initialize_lgamma_cache(call_var_opt_t *opt);
double log_bayes_factor(int k, int n, int alpha, int beta, double error_rate, const call_var_opt_t *opt);
double fisher_exact_test(int a, int b, int c, int d, const call_var_opt_t *opt);
double log_betabinom_pmf(int k, int n, int alpha, int beta, const call_var_opt_t *opt);
int median_int(int *arr, int n);
int min_int(int *arr, int n);


#endif