#include <stdio.h>
#include <math.h>
#include <fftw3.h>
#include <stdlib.h>
#include <complex.h>

#define NUM_POINTS 256

double sinc(double x) {
    if (x == 0.0) {
        return 1.0;
    } else {
        return sin(x) / x;
    }
}

// Function to plot two curves from given data using Gnuplot
void plotTwoCurves(double *x1, double *y1, double *x2, double *y2, int size) {
    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");

    if (gnuplotPipe) {
        // Send commands to Gnuplot to plot two curves
        fprintf(gnuplotPipe, "plot '-' with lines title 'Numerical Result', '-' with lines title 'Analytic Plot'\n");
        
        // Send data for the first curve
        for (int i = 0; i < size; ++i) {
            fprintf(gnuplotPipe, "%lf %lf\n", x1[i], y1[i]);
        }
        fprintf(gnuplotPipe, "e\n");
        
        // Send data for the second curve
        for (int i = 0; i < size; ++i) {
            fprintf(gnuplotPipe, "%lf %lf\n", x2[i], y2[i]);
        }
        fprintf(gnuplotPipe, "e\n");

        fflush(gnuplotPipe);
        printf("Press enter to exit...\n");
        getchar(); // Wait for user to press enter

        // Close the pipe
        pclose(gnuplotPipe);
    } else {
        printf("Error opening Gnuplot pipe.\n");
    }
}

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
double complex* calculate_FT(double *karr, double *xarr, fftw_complex *dft_arr, int size) {
    // Allocate memory for the result array
    double complex *result = (double complex*) malloc(sizeof(double complex) * size);

    // Perform element-wise multiplication, multiplication by i, and complex exponentiation
    for (int i = 0; i < size; i++) {
        double complex temp = -I * karr[i] * xarr[0];
        result[i] = (xarr[1] - xarr[0]) * sqrt(size / (2 * (M_PI))) * cpow(M_E, temp) * ((dft_arr[i][0]) + (I * dft_arr[i][1])) ;
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

// Function to calculate DFT using FFTW
fftw_complex* calculate_dft(double *input, int size) {
    // Allocate memory for input and output arrays
    fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size);
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size);

    // Create FFTW plan
    fftw_plan plan = fftw_plan_dft_1d(size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // Initialize input array
    for (int i = 0; i < size; i++) {
        in[i][0] = input[i];  // Real part
        in[i][1] = 0.0;       // Imaginary part
    }

    // Execute FFTW plan
    fftw_execute(plan);

    // Free memory and destroy plan
    fftw_destroy_plan(plan);
    fftw_free(in);

    // Return the DFT coefficients
    return out;
}

int main() {
    // Define the range of x values
    double x_min = -50.0;
    double x_max = 50.0;
    double step = (x_max - x_min) / (NUM_POINTS - 1);

    // Arrays to store x and y values
    double x_values[NUM_POINTS];
    double y_values[NUM_POINTS];
    double FT[NUM_POINTS];

    // Calculate the sinc values and store them in arrays
    for (int i = 0; i < NUM_POINTS; i++) {
        double x = x_min + i * step;
        double y = sinc(x);
        x_values[i] = x;
        y_values[i] = y;
    }

    // Calculate DFT
    fftw_complex *dft = calculate_dft(y_values, NUM_POINTS);

    for (int i = 0; i < NUM_POINTS; i++) {
        dft[i][0] = dft[i][0] / sqrt(NUM_POINTS);
        dft[i][1] = dft[i][1] / sqrt(NUM_POINTS);
    }

    // Calculate k values
    double *k_values = calculate_k(x_values, NUM_POINTS);

    // Calculate FT values
    complex *FT_values = calculate_FT(k_values, x_values, dft, NUM_POINTS);

    // Calculate analytic values
    double *f_analytic = boxFunction(k_values, NUM_POINTS);

    // Extracting k values and absolute values of FT to a file
    FILE *dataFile = fopen("data.txt", "w");
    if (dataFile == NULL) {
        printf("Error opening file.\n");
        return 1;
    }

    // Write k values and absolute values of FT to file
    for (int i = 0; i < NUM_POINTS; i++) {
        fprintf(dataFile, "%.5lf %.5lf\n", k_values[i], cabs(FT_values[i]));
    }

    // Close the file
    fclose(dataFile);

    // Free memory
    fftw_free(dft);

    return 0;
}

