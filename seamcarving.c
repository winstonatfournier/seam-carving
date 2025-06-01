#include "seamcarving.h"
#include "c_img.c"
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

double my_sqrt(double x) {
    if (x == 0) return 0;
    double guess = x / 2.0;
    double prev_guess;
    do {
        prev_guess = guess;
        guess = (guess + x / guess) / 2.0;
    } while (prev_guess != guess);
    return guess;
}

void calc_energy(struct rgb_img *im, struct rgb_img **grad) {
    *grad = (struct rgb_img *)malloc(sizeof(struct rgb_img));
    if (*grad == NULL){
        printf("*grad Memory allocation failed\n");
        exit(1);
    }
    ((*grad)->raster) = (uint8_t *)malloc(sizeof(uint8_t)*((im->width)*(im->height)));
    if (((*grad)->raster) == NULL){
        printf("*grad->raster Memory allocation failed\n");
        exit(1);
    }
    int *energies_x = (int *)malloc(sizeof(int)*(3*(im->width)*(im->height)));
    if (energies_x == NULL){
        printf("energies_x Memory allocation failed\n");
        exit(1);
    }
    int *energies_y = (int *)malloc(sizeof(int)*(3*(im->width)*(im->height)));
    if (energies_y == NULL){
        printf("energies_y Memory allocation failed\n");
        exit(1);
    }
    int *energies_z = (int *)calloc(((im->width)*(im->height)), sizeof(int));
    if (energies_z == NULL){
        printf("energies_z Memory allocation failed\n");
        exit(1);
    }
    for (int i = 0; i < im->height; i++){ // iterates through the rows
        for (int j = 0; j < im->width; j++){ // iterates through the columns
            for (int k = 0; k < 3; k++){ // iterates through RGB values for a given pixel
                if (j == 0){
                    energies_x[((3*im->width*i)+(3*j)+k)] = abs((im->raster)[((3*im->width*i)+(3*(j+1))+k)] - (im->raster)[((3*im->width*i)+(3*((im->width)-1))+k)]);
                } else if (j == ((im->width) - 1)){
                    energies_x[((3*im->width*i)+(3*j)+k)] = abs((im->raster)[((3*im->width*i)+(3*(0))+k)] - (im->raster)[((3*im->width*i)+(3*(j-1))+k)]);
                } else {
                    energies_x[((3*im->width*i)+(3*j)+k)] = abs((im->raster)[((3*im->width*i)+(3*(j+1))+k)] - (im->raster)[((3*im->width*i)+(3*(j-1))+k)]);
                }
                if (i == 0){
                    energies_y[((3*im->width*i)+(3*j)+k)] = abs((im->raster)[((3*im->width*((im->height)-1))+(3*j)+k)] - (im->raster)[((3*im->width*(i+1))+(3*j)+k)]);
                } else if (i == ((im->height) - 1)){
                    energies_y[((3*im->width*i)+(3*j)+k)] = abs((im->raster)[((3*im->width*(i-1))+(3*j)+k)] - (im->raster)[((3*im->width*(0))+(3*j)+k)]);
                } else {
                    energies_y[((3*im->width*i)+(3*j)+k)] = abs((im->raster)[((3*im->width*(i-1))+(3*j)+k)] - (im->raster)[((3*im->width*(i+1))+(3*j)+k)]);
                }
                //energies_z[((im->width*i)+j)] += pow((energies_x[((3*im->width*i)+(3*j)+k)]), 2);
                energies_z[((im->width*i)+j)] += ((energies_x[((3*im->width*i)+(3*j)+k)])*(energies_x[((3*im->width*i)+(3*j)+k)]));
                //energies_z[((im->width*i)+j)] += pow((energies_y[((3*im->width*i)+(3*j)+k)]), 2);
                energies_z[((im->width*i)+j)] += ((energies_y[((3*im->width*i)+(3*j)+k)])*(energies_y[((3*im->width*i)+(3*j)+k)]));
            }
        }
    }
    for (int l = 0; l < ((im->width)*(im->height)); l++){
        energies_z[l] = my_sqrt((energies_z[l]));
        energies_z[l] = (energies_z[l])/10;
        ((*grad)->raster)[l] = (uint8_t)energies_z[l];
    }
    (*grad)->height = im->height;
    (*grad)->width = im->width;

    free(energies_x);
    free(energies_y);
    free(energies_z);
    //free((*grad)->raster);
    //free(*grad);
}

void dynamic_seam(struct rgb_img *grad, double **best_arr){
    (*best_arr) = (double *)malloc((grad->width * grad->height) * sizeof(double));
    double ****comp_arr = (double ****)malloc(sizeof(double ***)); // ready comp_arr which'll compare costs across widths
    (*comp_arr) = (double ***)malloc(grad->width * sizeof(double **));
    double *cost_arr = (double *)malloc((grad->width * grad->height) * sizeof(double)); // ready cost_arr which'll record costs locally for a given seam
    int **path_counts = (int **)malloc(sizeof(int *)); // tracks number of seams drawn per bottom row 'width element'
    (*path_counts) = (int *)malloc(grad->width * sizeof(int));

    for (int j = 0; j < grad->width; j++){ // iterate for each bottom row element
        for (int f = 0; f < (grad->width * grad->height); f++){ // initialize cost_arr
            cost_arr[f] = 0;
        }
        cost_arr[((grad->height - 1) * grad->width) + j] = ((grad->raster)[((grad->height - 1) * grad->width) + j]); // put a bottom row element that'll begin the seam
        int path_y[] = {(grad->height - 1)}; // initialize paths
        int path_x[] = {j};
        int path_sz = 1;
        (*path_counts)[j] = 0;
        dynamic_seam_helper(grad, comp_arr, path_counts, cost_arr, (grad->height - 1), j, j, path_y, path_x, path_sz); // call helper
    } // NOTE: this is the structure for comp_arr as you'll see below; pointer -> bottom row elements -> potential paths for each bottom row element -> actual travel costs array for a given path
    for (int a = 0; a < (grad->width * grad->height); a++){ // iterate for every element in the travel cost array
        int min = 1000000; // very big number, should hopefully be bigger than any travel cost we might find
        for (int j = 0; j < grad->width; j++){ // iterate across every bottom row element
            for (int p = 0; p < ((*path_counts)[j]); p++){ // iterate for every path
                if (((*comp_arr)[j][p][a]) < min && ((*comp_arr)[j][p][a]) > 0){
                    min = ((*comp_arr)[j][p][a]); // 'min' records the smallest non-zero travel cost value for that given element to then assign the 'best_arr' matrix
                }
            }
        }
        (*best_arr)[a] = min;
    }
}

void dynamic_seam_helper(struct rgb_img *grad, double ****comp_arr, int **path_counts, double *cost_arr, int i, int j, int j_init, int *path_y, int *path_x, int path_sz){
    if (i < 1){ // base case
        (*path_counts)[j_init] += 1; // add another seam to comp_arr
        (*comp_arr)[j_init] = (double **)realloc((*comp_arr)[j_init], ((*path_counts)[j_init]) * sizeof(double *));
        (*comp_arr)[j_init][((*path_counts)[j_init]) - 1] = (double *)malloc((grad->width * grad->height) * sizeof(double));
        for (int a = 0; a < (grad->width * grad->height); a++){
            (*comp_arr)[j_init][((*path_counts)[j_init]) - 1][a] = (cost_arr)[a]; // copy data from cost_arr to comp_arr
        }
    } else {
        int *path_y_cpy = (int *)malloc(((path_sz) + 1) * sizeof(int)); // create path_cpy to create a local variable tracking the seam's path
        int *path_x_cpy = (int *)malloc(((path_sz) + 1) * sizeof(int));
        for (int c = 0; c < (path_sz); c++){ // copy data
            path_y_cpy[c] = path_y[c];
            path_x_cpy[c] = path_x[c];
        }
        int param1; // establish bounds for recursion based on position; accounting for edge cases
        int param2;
        if (j == 0){
            param1 = 0;
            param2 = 2;
        } else if (j == (grad->width - 1)){
            param1 = -1;
            param2 = 1;
        } else {
            param1 = -1;
            param2 = 2;
        }
        for (int b = param1; b < param2; b++){ // this loop will repeat at most three times, representing three more possible next steps in tracing the seam
            (path_y_cpy)[((grad->height) - i)] = (i - 1); // add next point to path variables
            (path_x_cpy)[((grad->height) - i)] = (j + b);
            uint8_t cost = ((grad->raster)[((i - 1) * grad->width) + (j + b)]); // note cost of travel in the given direction
            double *cost_arr_cpy = (double *)malloc((grad->width * grad->height) * sizeof(double)); // create a local copy of cost_arr
            for (int d = 0; d < (grad->width * grad->height); d++){ // copy data from cost_arr
                (cost_arr_cpy)[d] = (cost_arr)[d];
            }
            for (int k = 0; k < ((grad->height - i) + 1); k++){ // add cost of travel to every point along the path
                (cost_arr_cpy)[(((path_y_cpy)[k]) * grad->width) + ((path_x_cpy)[k])] += cost;
            }
            for (int d = 0; d < (grad->width * grad->height); d++){ // copy data from cost_arr
            }
            dynamic_seam_helper(grad, comp_arr, path_counts, cost_arr_cpy, (i - 1), (j + b), j_init, path_y_cpy, path_x_cpy, (path_sz + 1)); // invoke recursion
        }
    }
}

void recover_path(double *best, int height, int width, int **path){
    (*path) = (int *)malloc(height * sizeof(int));
    int min = 1000000;
    int min_pos = -1;
    for (int b = 0; b < (width); b++){
        if (best[((height - 1) * width) + b] < min){
            min = best[((height - 1) * width) + b];
            min_pos = b;
        }
    }
    (*path)[(height - 1)] = min_pos;
    recover_path_helper(best, height, width, (height - 1), min_pos, path);
}

void recover_path_helper(double *best, int height, int width, int i, int j, int **path){
    if (i >= 1){
        int j_new = 0;
        if (j == 0){
            if (best[((i - 1) * width) + (j + 0)] <= best[((i - 1) * width) + (j + 1)]){
                (*path)[i - 1] = j + 0;
                j_new += 0;
            } else if (best[((i - 1) * width) + (j + 1)] <= best[((i - 1) * width) + (j + 0)]){
                (*path)[i - 1] = j + 1;
                j_new += 1;
            }
        } else if (j == (width - 1)){
            if (best[((i - 1) * width) + (j - 1)] <= best[((i - 1) * width) + (j + 0)]){
                (*path)[i - 1] = j - 1;
                j_new -= 1;
            } else if (best[((i - 1) * width) + (j + 0)] <= best[((i - 1) * width) + (j - 1)]){
                (*path)[i - 1] = j + 0;
                j_new += 0;
            }
        } else {
            if (best[((i - 1) * width) + (j - 1)] <= best[((i - 1) * width) + (j + 0)] && best[((i - 1) * width) + (j - 1)] <= best[((i - 1) * width) + (j + 1)]){
                (*path)[i - 1] = j - 1;
                j_new -= 1;
            } else if (best[((i - 1) * width) + (j + 0)] <= best[((i - 1) * width) + (j - 1)] && best[((i - 1) * width) + (j + 0)] <= best[((i - 1) * width) + (j + 1)]){
                (*path)[i - 1] = j + 0;
                j_new += 0;
            } else if (best[((i - 1) * width) + (j + 1)] <= best[((i - 1) * width) + (j - 1)] && best[((i - 1) * width) + (j + 1)] <= best[((i - 1) * width) + (j + 0)]){
                (*path)[i - 1] = j + 1;
                j_new += 1;
            }
        }
        recover_path_helper(best, height, width, (i - 1), j + j_new, path);
    }
}

void remove_seam(struct rgb_img *src, struct rgb_img **dest, int *path) {
    size_t height = src->height;
    size_t width = src->width - 1;
    
    create_img(dest, height, width);

    for (int y = 0; y < height; y++) {
        int seam_x = path[y];
        int dest_x = 0;

        for (int x = 0; x < src->width; x++) {
            if (x == seam_x) {
                continue;
            }

            (*dest)->raster[3 * (y * width + dest_x) + 0] = src->raster[3 * (y * src->width + x) + 0]; // Red
            (*dest)->raster[3 * (y * width + dest_x) + 1] = src->raster[3 * (y * src->width + x) + 1]; // Green
            (*dest)->raster[3 * (y * width + dest_x) + 2] = src->raster[3 * (y * src->width + x) + 2]; // Blue

            dest_x++;
        }
    }
}