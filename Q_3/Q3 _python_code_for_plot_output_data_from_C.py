import numpy as np
import matplotlib.pyplot as plt

#loading data from c code
data=np.loadtxt('gsl_sinc_fft.txt')
k_value = data[:,0]
Numerical_fourier_transform = data[:,1]



# Define the given function 
def f(x):
    return np.sin(x) / x

# Define the range and number of points for x
xmin = -50
xmax = 50
n = 256

# Calculate the spacing between points
dx = (xmax - xmin) / (n - 1)

# Generate the array of x values
xx = np.linspace(xmin, xmax, n)

# Calculate the function values at the sample points
sample_value = f(xx)

# Perform numerical Fourier transformation
fft_sample = np.fft.fftshift(np.fft.fft(np.fft.fftshift(np.array(sample_value)), norm="ortho"))

# Generate the array of k values
kk = 2 * np.pi * np.fft.fftshift(np.fft.fftfreq(n, dx))

# Calculate the integration factor
int_factor = dx * np.sqrt(n / (2 * np.pi) * np.exp(- 1j * kk * xmin))

# Apply the integration factor to the Fourier transformed data
fft_sample = fft_sample * int_factor


# Plotting the results
plt.figure(figsize=(10, 5))
plt.plot(kk, np.abs(fft_sample), label='Fourier transform using Python code')
plt.plot(k_value, Numerical_fourier_transform, linestyle = "--", label = "Fourier transform using C code")
plt.xlabel('k')
plt.ylabel('|F(k)|')
plt.title('Numerical Fourier Transform Python vs C code (using GSL)')
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
