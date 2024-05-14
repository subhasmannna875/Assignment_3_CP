import numpy as np
import matplotlib.pyplot as plt

# loding data from c code
data_analytic=np.loadtxt('k_analytic_values.txt')
data_numerical=np.loadtxt('k_ft_values.txt')
k_value = data_analytic[:,0]
Analytic_fourier_transform = data_analytic[:,1]
Numerical_fourier_transform = data_numerical[:,1]


# Plotting the results
plt.figure(figsize=(10, 5))
plt.plot(k_value, Numerical_fourier_transform,  label = "Numerical fourier transform")
plt.plot(k_value, Analytic_fourier_transform, linestyle = "--", label='Analytical fourier transform')
plt.xlabel('k')
plt.ylabel('|F(k)|')
plt.title('Fourier Transform of Gaussian function using C code')
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()

