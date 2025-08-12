#include <math.h>
#include <float.h>
#include <stdio.h>
#include "call_var_main.h"

void initialize_lgamma_cache(call_var_opt_t *opt) {
    opt->min_lgamma_i = 0, opt->max_lgamma_i = LONGCALLD_LGAMMA_MAX_I;
    for (int i = opt->min_lgamma_i; i <= opt->max_lgamma_i; i++) {
        opt->lgamma_cache[i] = lgamma(i);
    }
}

double fast_lgamma(int x, const call_var_opt_t *opt) {
    if (x >= opt->min_lgamma_i && x <= opt->max_lgamma_i) {
        return opt->lgamma_cache[x];
    }
    return lgamma(x);
} 

int log_likelilood(double prob) {
    if (prob < 0.00001) return -1000000;
    else if (prob > 0.99999) return 1000000;
    return log(prob);
}

// compare_int
int compare_int(const void *a, const void *b) {
    return (*(int *)a - *(int *)b);
}

int median_int(int *arr, int n) {
    if (n <= 0) return 0; // handle empty array case
    // sort arr
    qsort(arr, n, sizeof(int), (int (*)(const void *, const void *))compare_int);
    int med = n / 2;
    if (n % 2 == 0) {
        // return (arr[med - 1] + arr[med]) / 2; // average of two middle elements
        return arr[med - 1]; // return the lower middle element for even length
    } else {
        return arr[med]; // middle element
    }
}

int min_int(int *arr, int n) {
    if (n <= 0) return 0; // handle empty array case
    int min = arr[0];
    for (int i = 1; i < n; ++i) {
        if (arr[i] < min) min = arr[i];
    }
    return min;
}

// XXX pre-compute lgamma(alpha/beta/alpha+beta) for efficiency

// Log-Beta function: log(B(α, β))
double log_beta(int alpha, int beta, const call_var_opt_t *opt) {
    static int last_alpha = 0, last_beta = 0, last_result = 0;
    if (alpha == last_alpha && beta == last_beta) {
        return last_result;
    }
    return fast_lgamma(alpha, opt) + fast_lgamma(beta, opt) - fast_lgamma(alpha + beta, opt);
}

// Log-Binomial PMF: log(P(k | n, θ))
double log_binom_pmf(int k, int n, double theta, const call_var_opt_t *opt) {
    if (k < 0 || k > n) return -INFINITY;
    if (theta == 0.0) return (k == 0) ? 0.0 : -INFINITY;
    if (theta == 1.0) return (k == n) ? 0.0 : -INFINITY;
    static int last_n = -1, last_k = -1;
    static double log_comb = 0;
    
    if (n != last_n || k != last_k) {
        log_comb = fast_lgamma(n + 1, opt) - fast_lgamma(k + 1, opt) - fast_lgamma(n - k + 1, opt);
        last_n = n;
        last_k = k;
    }
    
    return log_comb + k * log(theta) + (n - k) * log1p(-theta);
}

// Log-Beta-Binomial PMF: log(P(k | n, α, β))
double log_betabinom_pmf(int k, int n, int alpha, int beta, const call_var_opt_t *opt) {
    // fprintf(stderr, "log_betabinom_pmf: k=%d, n=%d, alpha=%d, beta=%d\n", k, n, alpha, beta);
    return fast_lgamma(n + 1, opt) - fast_lgamma(k + 1, opt) - fast_lgamma(n - k + 1, opt)
           + log_beta(k + alpha, n - k + beta, opt) - log_beta(alpha, beta, opt);
}

// Bayes Factor Calculation: log(BF) = log(P(k | n, α, β)) - log(P(k | n, θ))
double log_bayes_factor(int k, int n, int alpha, int beta, double error_rate, const call_var_opt_t *opt) {
    double log_H1 = log_betabinom_pmf(k, n, alpha, beta, opt);
    double log_H0 = log_binom_pmf(k, n, error_rate, opt);
    printf("log_H1: %.5f, log_H0: %.5f\n", log_H1, log_H0); // Debugging line to check values
    printf("H1: %.5f, H0: %.10f\n", exp(log_H1), exp(log_H0)); // Debugging line to check exponentials
    return log_H1 - log_H0;  // log(BF)
    return log_beta(k + alpha, n - k + beta, opt) // log(B(k+α, n-k+β))
           - log_beta(alpha, beta, opt) // log(B(α, β))
              - k * log(error_rate) - (n - k) * log1p(-error_rate); // log(P(k | n, θ))
}

// Helper function to compute hypergeometric probability in log space
double log_hypergeometric(int a, int b, int c, int d, const call_var_opt_t *opt) {
    // Calculate sums once
    const int n1 = a + b;
    const int n2 = c + d;
    const int m1 = a + c;
    const int m2 = b + d;
    const int N = n1 + n2;
    
    // Use symmetry to minimize calculations
    if (n1 > n2) return log_hypergeometric(c, d, a, b, opt);
    if (m1 > m2) return log_hypergeometric(b, a, d, c, opt);
    
    // Calculate using cached lgamma values
    return fast_lgamma(n1 + 1, opt) + fast_lgamma(n2 + 1, opt) + fast_lgamma(m1 + 1, opt) + fast_lgamma(m2 + 1, opt) -
           (fast_lgamma(a + 1, opt) + fast_lgamma(b + 1, opt) + fast_lgamma(c + 1, opt) + fast_lgamma(d + 1, opt) + fast_lgamma(N + 1, opt));
}

// Fisher's exact test (two-tailed)
double fisher_exact_test(int a, int b, int c, int d, const call_var_opt_t *opt) {
    // Compute the probability of the observed table
    double log_p_observed = log_hypergeometric(a, b, c, d, opt);
    double p_observed = exp(log_p_observed);
    
    // For two-tailed test, we sum all probabilities <= p_observed
    double total_p = 0.0;
    
    // Find the minimum and maximum possible values for a
    int min_a = (0 > (a + c) - (a + b + c + d)) ? 0 : (a + c) - (b + d);
    int max_a = (a + b) < (a + c) ? (a + b) : (a + c);
    // Calculate the mode of the hypergeometric distribution
    int mode_a = (int)((a + b) * (a + c) / (double)(a + b + c + d));
    
    // Evaluate tables from mode outward
    for (int delta = 0; delta <= max_a - min_a; delta++) {
        // Check tables above mode
        int current_a = mode_a + delta;
        if (current_a <= max_a) {
            int current_b = (a + b) - current_a;
            int current_c = (a + c) - current_a;
            int current_d = (b + d) - current_b;
            
            if (current_b >= 0 && current_c >= 0 && current_d >= 0) {
                double log_p = log_hypergeometric(current_a, current_b, current_c, current_d, opt);
                double p = exp(log_p);
                if (p <= p_observed + DBL_EPSILON) {
                    total_p += p;
                }
            }
        }
        
        // Check tables below mode (skip mode itself on second iteration)
        if (delta > 0) {
            current_a = mode_a - delta;
            if (current_a >= min_a) {
                int current_b = (a + b) - current_a;
                int current_c = (a + c) - current_a;
                int current_d = (b + d) - current_b;
                
                if (current_b >= 0 && current_c >= 0 && current_d >= 0) {
                    double log_p = log_hypergeometric(current_a, current_b, current_c, current_d, opt);
                    double p = exp(log_p);
                    if (p <= p_observed + DBL_EPSILON) {
                        total_p += p;
                    }
                }
            }
        }
    }
    return total_p;
}

// Fisher's exact test (one-tailed - less)
double fisher_exact_test_left(int a, int b, int c, int d, const call_var_opt_t *opt) {
    double log_p_observed = log_hypergeometric(a, b, c, d, opt);
    double total_p = 0.0;
    
    int min_a = (0 > (a + c) - (a + b + c + d)) ? 0 : (a + c) - (b + d);
    
    for (int current_a = min_a; current_a <= a; current_a++) {
        int current_b = (a + b) - current_a;
        int current_c = (a + c) - current_a;
        int current_d = (b + d) - current_b;
        
        if (current_b >= 0 && current_c >= 0 && current_d >= 0) {
            double log_p = log_hypergeometric(current_a, current_b, current_c, current_d, opt);
            total_p += exp(log_p);
        }
    }
    return total_p;
}

// Fisher's exact test (one-tailed - greater)
double fisher_exact_test_right(int a, int b, int c, int d, const call_var_opt_t *opt) {
    double log_p_observed = log_hypergeometric(a, b, c, d, opt);
    double total_p = 0.0;
    
    int max_a = (a + b) < (a + c) ? (a + b) : (a + c);
    
    for (int current_a = a; current_a <= max_a; current_a++) {
        int current_b = (a + b) - current_a;
        int current_c = (a + c) - current_a;
        int current_d = (b + d) - current_b;
        
        if (current_b >= 0 && current_c >= 0 && current_d >= 0) {
            double log_p = log_hypergeometric(current_a, current_b, current_c, current_d, opt);
            total_p += exp(log_p);
        }
    }
    return total_p;
}

#ifdef _MATH_UTILS_TEST
int main() {
    int k = 8;          // Alternate allele count
    int n = 50;        // Total reads
    int alpha = 2; // Beta prior (somatic)
    int beta = 10; // Beta prior (background)
    double error_rate = 0.01; // Sequencing error

    double logBF = log_bayes_factor(k, n, alpha, beta, error_rate);
    double BF = exp(logBF);

    printf("Log-Bayes Factor: %.2f\n", logBF);
    printf("Bayes Factor: %.2e\n", BF);
    printf("Decision: %s\n", BF > 100 ? "Somatic" : "Error/Germline");

    // test fisher exact test
    int a = 5, b = 10, c = 15, d = 2;
    printf("Fisher's Exact Test (Two-tailed): %.5f\n", fisher_exact_test(a, b, c, d));

    return 0;
}
#endif