import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

# Define x, y, and xnew to resample at.
x = np.linspace(0, 500, num=10, endpoint=True)
y = np.cos(-x**2/9.0)
xnew = np.linspace(0, 500, num=200, endpoint=True)

# Define interpolators.
f_linear = interp1d(x, y)
f_cubic = interp1d(x, y, kind='cubic')

# Plot.
plt.plot(x, y, 'o', label='data')
plt.plot(xnew, f_linear(xnew), '-', label='linear',color = 'red')
plt.plot(xnew, f_cubic(xnew), '-', label='cubic',color = 'green')
plt.legend(loc='best')
plt.show()