// gcc -o kmedoids kmedoids.c -lm; cp kmedoids ~/bin/
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

#define DIM 20          // Dimensionality of each data point
#define MAX_ITERATIONS 100  // Maximum number of iterations for swap phase
int N = 0;               // Number of data points

// Function to calculate Manhattan distance between two points
int manhattan_distance(const int *a, const int *b, int dim) {
    int dist = 0;
    for (int i = 0; i < dim; i++) {
        dist += abs(a[i] - b[i]);
    }
    return dist;
}

// Function to initialize medoids
void initialize_medoids(int *medoids, int n, int k) {
    for (int i = 0; i < k; i++) {
        medoids[i] = i;  // Use the first K points as initial medoids
    }
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
        fprintf(stderr, "Point %d: %d -> %d\n", i, assignment[i], closest_medoid);
        assignment[i] = closest_medoid;
        total_cost += min_dist;
    }
    total_costs[medoids[0]][medoids[1]] = total_cost;
    total_costs[medoids[1]][medoids[0]] = total_cost;
    return total_cost;
}

// Assign points to nearest medoid and calculate total cost
// int assign_points(int *assignment, int *medoids, int **data, int n, int dim, int k) {
//     int total_cost = 0;
//     for (int i = 0; i < n; i++) {
//         int min_dist = INT_MAX;
//         int closest_medoid = -1;
//         for (int j = 0; j < k; j++) {
//             int dist = manhattan_distance(data[i], data[medoids[j]], dim);
//             if (dist < min_dist) {
//                 min_dist = dist;
//                 closest_medoid = medoids[j];
//             }
//         }
//         // if (assignment[i] != closest_medoid) {
//         printf("Point %d: %d -> %d\n", i, assignment[i], closest_medoid);
//         // }
//         assignment[i] = closest_medoid;
//         total_cost += min_dist;
//     }
//     return total_cost;
// }

// Swap Phase with max iterations to optimize medoids
// int swap_phase(int *medoids, int *assignment, int **data, int n, int dim, int k, int max_iterations) {
//     int changed, iterations = 0;
//     int final_best_cost = INT_MAX;
//     do {
//         changed = 0;
//         iterations++;
//         for (int i = 0; i < k; i++) {
//             int best_cost = assign_points(assignment, medoids, data, n, dim, k);
//             int best_medoid = medoids[i];

//             for (int j = 0; j < n; j++) {
//                 if (assignment[j] == medoids[i]) continue;

//                 int temp = medoids[i];
//                 medoids[i] = j;
//                 int new_cost = assign_points(assignment, medoids, data, n, dim, k);

//                 if (new_cost < best_cost) {
//                     best_cost = new_cost;
//                     best_medoid = j;
//                     printf("best_assign: %d\n", best_cost);
//                     for (int l = 0; l < n; l++) {
//                         printf("%d ", assignment[l]);
//                     } printf("\n");
//                     changed = 1;
//                 }
//                 medoids[i] = temp;
//             }
//             if (final_best_cost > best_cost) final_best_cost = best_cost;
//             medoids[i] = best_medoid;
//         }
//     } while (changed && iterations < max_iterations);

//     printf("Medoid:\n");
//     for (int j = 0; j < k; j++) {
//         printf("%d: %d\t", j, medoids[j]);
//         for (int i = 0; i < dim; i++) {
//             printf("%5d ", data[medoids[j]][i]);
//         } printf("\n");
//     }

//     printf("%d iter\t", iterations);
//     return final_best_cost;
// }

// Swap Phase with max iterations to optimize medoids
int swap_2medoids(int *medoids, int *assignment, int *best_assignment, int **data, int n, int dim, int k, int max_iterations) {
    int changed, iterations = 0;
    int final_best_cost = INT_MAX;
    int **manhattan_distance_array = (int **)malloc(n * sizeof(int *));
    int **total_cost = (int **)malloc(n * sizeof(int *));
    // int *cur_medoids = (int *)malloc(k * sizeof(int));
    // int *assignment = (int *)malloc(n * sizeof(int));
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
            // for (int j = 0; j < k; j++) cur_medoids[j] = medoids[j];
            // printf("cur_medoids: %d %d\n", cur_medoids[0], cur_medoids[1]);
            int best_cost = assign_points(assignment, medoids, data, n, dim, k, manhattan_distance_array, total_cost, 0);
            fprintf(stderr, "iter%d\tclu%d: medoids: %d,  new_assign: %d\n", iterations, i, medoids[i], best_cost);
            int best_medoid = medoids[i];

            // check if there is a better medoid other than i
            for (int j = 0; j < n; j++) {
                int is_medoid = 0;
                for (int l = 0; l < k; l++) {
                    if (medoids[l] == j) {
                        is_medoid = 1;
                        break;
                    }
                }
                if (is_medoid) continue;
                int temp = medoids[i];
                medoids[i] = j;
                int new_cost = assign_points(assignment, medoids, data, n, dim, k, manhattan_distance_array, total_cost, 0);
                fprintf(stderr, "iter%d\tclu%d: medoids: %d(%d),  new_assign: %d\n", iterations, i, j, medoids[k-i-1], new_cost);

                if (new_cost < best_cost) {
                    best_cost = new_cost;
                    best_medoid = j;
                    printf("iter%d\tclu%d: best_assign: %d\n", iterations, i, best_cost);
                    for (int l = 0; l < n; l++) {
                        printf("%d ", assignment[l]);
                    } printf("\n");
                    changed = 1;
                }
                medoids[i] = temp;
            }
            if (final_best_cost > best_cost) final_best_cost = best_cost;
            medoids[i] = best_medoid;
        }
    } while (changed && iterations < max_iterations);

    printf("Medoid:\n");
    for (int j = 0; j < k; j++) {
        printf("%d: %d\t", j, medoids[j]);
        for (int i = 0; i < dim; i++) {
            printf("%5d ", data[medoids[j]][i]);
        } printf("\n");
    }

    printf("%d iter\n", iterations);

    assign_points(best_assignment, medoids, data, n, dim, k, manhattan_distance_array, total_cost, 1);
    for (int i = 0; i < n; i++) {
        free(manhattan_distance_array[i]);
        free(total_cost[i]);
    } 
    free(manhattan_distance_array); free(total_cost);
    // free(assignment);
    return final_best_cost;
}

// Function to load data from a file
int **load_data(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Failed to open file");
        exit(EXIT_FAILURE);
    }

    int ch;
    while (!feof(file)) {
        ch = fgetc(file);
        if (ch == '\n') N++;
    }
    rewind(file);

    int **data = (int **)malloc(N * sizeof(int *));
    for (int i = 0; i < N; i++) {
        data[i] = (int *)malloc(DIM * sizeof(int));
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < DIM; j++) {
            if (fscanf(file, "%d", &data[i][j]) != 1) {
                fprintf(stderr, "Error reading data at row %d, column %d\n", i + 1, j + 1);
                exit(EXIT_FAILURE);
            }
        }
    }
    fclose(file);
    return data;
}

// Function to calculate silhouette score
double silhouette_score(int *assignment, int **data, int n, int dim, int k) {
    double total_score = 0.0;

    // Loop over each point to calculate its silhouette score
    for (int i = 0; i < n; i++) {
        int current_cluster = assignment[i];

        // Calculate a(i): average intra-cluster distance
        double intra_cluster_dist = 0.0;
        int intra_count = 0;
        for (int j = 0; j < n; j++) {
            if (i != j && assignment[j] == current_cluster) {
                intra_cluster_dist += manhattan_distance(data[i], data[j], dim);
                intra_count++;
            }
        }
        double a_i = (intra_count > 0) ? intra_cluster_dist / intra_count : 0.0;

        // Calculate b(i): minimum average distance to points in a different cluster
        double min_inter_cluster_dist = INT_MAX;
        for (int other_cluster = 0; other_cluster < k; other_cluster++) {
            if (other_cluster == current_cluster) continue;

            double inter_cluster_dist = 0.0;
            int inter_count = 0;
            for (int j = 0; j < n; j++) {
                if (assignment[j] == other_cluster) {
                    inter_cluster_dist += manhattan_distance(data[i], data[j], dim);
                    inter_count++;
                }
            }
            if (inter_count > 0) {
                double avg_inter_dist = inter_cluster_dist / inter_count;
                if (avg_inter_dist < min_inter_cluster_dist) {
                    min_inter_cluster_dist = avg_inter_dist;
                }
            }
        }
        double b_i = min_inter_cluster_dist;

        // Silhouette score for point i
        double s_i = (b_i - a_i) / fmax(a_i, b_i);
        total_score += s_i;
    }

    // Return the average silhouette score for all points
    return total_score / n;
}

int main(int argc, char *argv[]) {
    int **data = load_data(argv[1]);
    int *assignment = calloc(N, sizeof(int));
    int *best_assignment = calloc(N, sizeof(int));

    int optimal_k = 1;
    int min_cost = INT_MAX;
    double max_silhouette_score = -1.0;

    for (int k = 2; k <= 2; k++) {
        printf("\n\nk=%d\n", k);
        int medoids[k];
        initialize_medoids(medoids, N, k);

        // assign_points(assignment, medoids, data, N, DIM, k);
        // int score = swap_phase(medoids, assignment, data, N, DIM, k, MAX_ITERATIONS);
        int score = swap_2medoids(medoids, assignment, best_assignment, data, N, DIM, k, MAX_ITERATIONS);
        for (int j = 0; j < N; j++) {
            printf("%d ", best_assignment[j]);
        } printf("\n");

        double sscore = silhouette_score(best_assignment, data, N, DIM, k);
        printf("Score for k=%d: %d\n", k, score);
        printf("Silhouette score for k=%d: %.10f\n", k, sscore);

        if (sscore > max_silhouette_score) {
            max_silhouette_score = sscore;
            optimal_k = k;
        }
    }

    printf("Optimal number of clusters (k): %d with silhouette score: %.10f\n", optimal_k, max_silhouette_score);


    for (int i = 0; i < N; i++) {
        free(data[i]);
    } free(data); free(assignment); free(best_assignment);
    return 0;
}