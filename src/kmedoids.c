#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

extern int LONGCALLD_VERBOSE;
#define MAX_ITERATIONS 100  // Maximum number of iterations for swap phase

// Function to calculate Manhattan distance between two points
int manhattan_distance(const int *a, const int *b, int dim) {
    int dist = 0;
    for (int i = 0; i < dim; i++) {
        dist += abs(a[i] - b[i]);
    }
    return dist;
}

// Assign points to nearest medoid and calculate total cost
int assign_points(int *assignment, int *medoids, int **data, int n, int dim, int k, int **manhattan_distance_array, int **total_costs, int force_assign) {
    if (!force_assign) {
        if (total_costs[medoids[0]][medoids[1]] != -1) {
            return total_costs[medoids[0]][medoids[1]];
        }
    }
    int total_cost = 0;
    for (int i = 0; i < n; i++) {
        if (data[i][0] == -1) {
            assignment[i] = -1;
            continue;
        }
        int min_dist = INT_MAX;
        int closest_medoid = -1;
        for (int j = 0; j < k; j++) {
            if (manhattan_distance_array[i][medoids[j]] == -1)
                manhattan_distance_array[i][medoids[j]] = manhattan_distance(data[i], data[medoids[j]], dim);
            int dist = manhattan_distance_array[i][medoids[j]];
            if (dist < min_dist) {
                min_dist = dist;
                closest_medoid = medoids[j];
            }
        }
        assignment[i] = closest_medoid;
        total_cost += min_dist;
    }
    total_costs[medoids[0]][medoids[1]] = total_cost;
    total_costs[medoids[1]][medoids[0]] = total_cost;
    return total_cost;
}

int swap_2medoids(int **data, int *medoids, int *best_assignment, int n, int dim, int k, int max_iterations) {
    int changed, iterations = 0;
    int final_best_cost = INT_MAX;
    int *assignment = (int *)malloc(n * sizeof(int));
    int **manhattan_distance_array = (int **)malloc(n * sizeof(int *));
    int **total_cost = (int **)malloc(n * sizeof(int *));
    for (int i = 0; i < n; i++) {
        manhattan_distance_array[i] = (int *)malloc(n * sizeof(int));
        total_cost[i] = (int *)malloc(n * sizeof(int));
        for (int j = 0; j < n; j++) {
            if (i == j) manhattan_distance_array[i][j] = 0;
            else manhattan_distance_array[i][j] = -1;  // manhattan_distance(data[i], data[j], dim);
            total_cost[i][j] = -1;
        }
    }
    do {
        changed = 0;
        iterations++;
        for (int i = 0; i < k; i++) {
            int best_cost = assign_points(assignment, medoids, data, n, dim, k, manhattan_distance_array, total_cost, 0);
            int best_medoid = medoids[i];

            // check if there is a better medoid other than i
            for (int j = 0; j < n; j++) {
                if (data[j][0] == -1) continue;
                for (int l = 0; l < k; l++) {
                    if (medoids[l] == j) continue;
                }
                int temp = medoids[i];
                medoids[i] = j;
                int new_cost = assign_points(assignment, medoids, data, n, dim, k, manhattan_distance_array, total_cost, 0);

                if (new_cost < best_cost) {
                    best_cost = new_cost;
                    best_medoid = j;
                    changed = 1;
                }
                medoids[i] = temp;
            }
            if (final_best_cost > best_cost) final_best_cost = best_cost;
            medoids[i] = best_medoid;
        }
    } while (changed && iterations < max_iterations);


    assign_points(best_assignment, medoids, data, n, dim, k, manhattan_distance_array, total_cost, 1);

    for (int i = 0; i < n; i++) {
        free(manhattan_distance_array[i]);
        free(total_cost[i]);
    } 
    free(manhattan_distance_array); free(total_cost); free(assignment);
    return final_best_cost;
}

// haps: 1/2: hap, 0: unasigned
int init_2medoids(int **xid_profile, int *medoids, int n_medoids, int *haps, int n_reads, int use_hap_info) {
    if (use_hap_info) {
        for (int i = 0; i < n_medoids; ++i) {
            for (int j = 0; j < n_reads; ++j) {
                if (xid_profile[j][0] == -1) continue;
                if (haps[j] == i+1) {
                    medoids[i] = j;
                    break;
                }
            }
        }
    } else { // select initial medoids randomly (first k reads)
        for (int i = 0; i < n_medoids; i++) medoids[i] = i;
    }
    return 0;
}

// if use_hap_info: select initial medoids based on haplotype information
// otherwise: select initial medoids randomly (first k reads)
// up to 2 medoids, if 2 medoids are too similar, only one medoid will be selected
// XXX handle cases where one of the medoids is outlier, i.e., the cluster size of one medoid is too small
// remove outlier, re-assign medoids, may need to re-assign medoids multiple times
// stop: 1) cluster size are unbanlanced but two medoids are similar to each other, 2) cluster size are balanced
int *xid_profile_2medoids(int **xid_profile, int n_dim, int n_medoids, int *haps, int n_reads, int use_hap_info) {
    int *medoids = (int*)malloc(n_medoids * sizeof(int));
    int *best_assignment = (int *)malloc(n_reads * sizeof(int));
    init_2medoids(xid_profile, medoids, n_medoids, haps, n_reads, use_hap_info);
    swap_2medoids(xid_profile, medoids, best_assignment, n_reads, n_dim, n_medoids, MAX_ITERATIONS);
    int *cluster_size = (int *)calloc(n_medoids, sizeof(int));

    if (LONGCALLD_VERBOSE >= 2) {
        for (int i = 0; i < n_reads; i++) {
            fprintf(stderr, "%d ", best_assignment[i]);
            for (int j = 0; j < n_medoids; ++j) {
                if (medoids[j] == best_assignment[i]) {
                    cluster_size[j]++;
                    break;
                }
            }
        } fprintf(stderr, "\nCluster ");
        for (int i = 0; i < n_medoids; i++) {
            fprintf(stderr, "%d: %d\t", i, cluster_size[i]);
        } fprintf(stderr, "\n");
    }

    free(best_assignment);
    return medoids;
}