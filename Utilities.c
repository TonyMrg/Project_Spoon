#include "global_parameters.h"

void row_to_col_major(double* src, double* dst, int rows, int cols) {
    int i, j;
#pragma omp parallel for collapse(2) private(i, j)
    for (i = 0; i < rows; ++i) {
        for (j = 0; j < cols; ++j) {
            dst[j * rows + i] = src[i * cols + j];
        }
    }
}

void save_energy(double* energy, int n) {
    FILE* f = fopen("C://Users//kosta//Downloads//energy.csv", "w");
    if (!f) return;
    for (int i = 0; i < n; i++) {
        fprintf(f, "%d,%f\n", i, energy[i]);
    }
    fclose(f);
}

void write_results_txt(double* results, int steps, double dt) {
    int results_size = steps * nom * 2;
    char txt_file_path[200];
    FILE* fp = NULL;

    while (fp == NULL) {
        printf("Input .txt file path: ");
        fgets(txt_file_path, sizeof(txt_file_path), stdin);
        txt_file_path[strcspn(txt_file_path, "\n")] = '\0';
        size_t len = strlen(txt_file_path);
        if (len > 1 && txt_file_path[0] == '"' && txt_file_path[len - 1] == '"') {
            memmove(txt_file_path, txt_file_path + 1, len - 2);
            txt_file_path[len - 2] = '\0';
        }
        
        errno_t err = fopen_s(&fp, txt_file_path, "w");
        if (err != 0 || fp == NULL) {
            printf("Couldn't open .txt file. Try again.\n");
        }
    }

    fprintf(fp, "%d,%d,%lf\n", nom, steps, dt);
    for (int i = 0; i < results_size; i++) {
        fprintf(fp, "%lf\n", results[i]);
    }
    fclose(fp);
}

double* read_results_txt(int* nom, int* steps, double* dt) {    
    char txt_file_path[200];
    FILE* fp = NULL;
    while (fp == NULL) {
        printf("Input .txt file path: ");
        fgets(txt_file_path, sizeof(txt_file_path), stdin);
        txt_file_path[strcspn(txt_file_path, "\n")] = '\0';
        size_t len = strlen(txt_file_path);
        if (len > 1 && txt_file_path[0] == '"' && txt_file_path[len - 1] == '"') {
            memmove(txt_file_path, txt_file_path + 1, len - 2);
            txt_file_path[len - 2] = '\0';
        }

        errno_t err = fopen_s(&fp, txt_file_path, "r");
        if (err != 0 || fp == NULL) {
            printf("Couldn't open .txt file. Try again.\n");
        }
    }

    if (fscanf_s(fp, "%d,%d,%lf", nom, steps, dt) != 3) {
        fprintf(stderr, "Error reading header\n");
        fclose(fp);
        return EXIT_FAILURE;
    }
    int total = (*steps) * (*nom);
    double* results = malloc(total * sizeof(double)); if (!results) {fclose(fp);return EXIT_FAILURE;}
    int index = 0;
    while (index < total && fscanf_s(fp, "%lf", &results[index]) == 1) {
        index++;
    }
    
    fclose(fp);
    return results;
}
