#include <stdio.h>
#include <math.h>
#include <gsl/gsl_fft_complex.h>
#include <stdlib.h>
#include <complex.h>

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

// Function to generate a box function
double* boxFunction(double *x, int size) {
    double* y = malloc(size * sizeof(double));
    if (y == NULL) {
        printf("Memory allocation failed\n");
        exit(1);
    }
    for (int i = 0; i < size; ++i) {
        if (x[i] >= -1.0 && x[i] <= 1.0) {
            y[i] = sqrt(M_PI/2.0);
        } else {
            y[i] = 0.0;
        }
    }
    return y;
}

// Function to roll an array of doubles by half of its length and return the rolled array
double* roll_array(double arr[], int size) {
    // Allocate memory for the rolled array
    double *rolled_array = (double*) malloc(sizeof(double) * size);

    // Find the mid point
    int mid = size / 2;

    // Roll the array
    for (int i = 0; i < size; i++) {
        rolled_array[i] = arr[(i + mid) % size];
    }

    return rolled_array;
}

// Function to perform element-wise complex multiplication and exponentiation
double complex* calculate_FT(double *karr, double *xarr, double *dft_arr, int size) {
    // Allocate memory for the result array
    double complex *result = (double complex*) malloc(sizeof(double complex) * size);

    // Perform element-wise multiplication, multiplication by i, and complex exponentiation
    for (int i = 0; i < size; i++) {
        double complex temp = -I * karr[i] * xarr[0];
        result[i] = (xarr[1] - xarr[0]) * sqrt(size / (2 * (M_PI))) * cpow(M_E, temp) * ((dft_arr[2 * i]) + (I * dft_arr[2 * i + 1]));
    }

    return result;
}

// Function to calculate 2*pi*x/(n*dx^2) for each element in the array
double* calculate_k(double *x, int n) {
    // Allocate memory for the result array
    double *result = (double*) malloc(sizeof(double) * n);

    // Calculate expression for each element
    for (int i = 0; i < n; i++) {
        result[i] = 2 * M_PI * (i - n / 2) / (n * (x[1] - x[0]));
    }

    // Return the array containing the calculated values
    return roll_array(result, n);
}

int main() {
    int i;
    double x;
    int n = 256;
    double x_values[n];
    double y_values[n];
    double data[2 * n];
    double x_min = -50;
    double x_max = 50;
    double FT[n];
    double delta = (x_max - x_min) / (n - 1);

    FILE *outputFile;
    outputFile = fopen("gsl_sinc_fft.txt", "w");
    if (outputFile == NULL) {
        printf("Error opening output file.\n");
        return 1;
    }

    for (i = 0; i < n; i++) {
        x = x_min + i * delta;
        x_values[i] = x;
        y_values[i] = sin(x)/x;
        REAL(data, i) = sin(x)/x;
        IMAG(data, i) = 0.0;
    }

    gsl_fft_complex_radix2_forward(data, 1, n);

    for (i = 0; i < n; i++) {
        data[2 * i] = data[2 * i] / sqrt(n);
        data[2 * i + 1] = data[2 * i + 1] / sqrt(n);
    }

    // Calculate k values
    double *k_values = calculate_k(x_values, n);
    // Calculate FT values
    complex *FT_values = calculate_FT(k_values, x_values, data, n);

    // Calculate analytic values
    double *f_analytic = boxFunction(k_values, n);

    // Write k values and absolute values of FT to file
    for (int i = 0; i < n; i++) {
        fprintf(outputFile, "%.5lf %.5lf\n", k_values[i], cabs(FT_values[i]));
    }

    // Close the file
    fclose(outputFile);

    return 0;
}
