from scipy.interpolate import make_interp_spline
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as dates
from datetime import datetime
from scipy.interpolate import interp1d



import pandas as pd
colums = ['1','2','3','4','5','6','7','8','9']
df  = pd.read_excel('A:\PhD\Machin learning\Book1.xlsx', names = colums)
x = df['1']
y1 = df['2']
y2 = df['3']
y3 = df['4']
y4 = df['5']
y5 = df['6']
y6 = df['7']
y7 = df['8']
y8 = df['9']

new_x = []
new_y1 = []
new_y2 = []
new_y3 = []
new_y4 = []
new_y5 = []
new_y6 = []
new_y7 = []
new_y8 = []

for i in range(0,len(x),4):
    #new_x.append(x[i])
    #new_x.append((x[i+1]+x[i-1])/2)
    new_y1.append(y1[i])
    new_y2.append(y2[i])
    new_y1.append(y3[i])
    new_y1.append(y4[i])
    new_y1.append(y5[i])
    new_y1.append(y6[i])
    new_y1.append(y7[i])
    new_y1.append(y8[i])

    #new_y.append((y[i+1]+y[i-1])/2)
xnew = np.linspace(min(x), max(x), num=500, endpoint=True)

# Define interpolators.
f_linear = interp1d(new_x, new_y1)
f_cubic = interp1d(new_x, new_y1, kind='cubic')
f_cubic2 = interp1d(new_x, new_y2, kind='cubic')
f_cubic3 = interp1d(new_x, new_y3, kind='cubic')
f_cubic4 = interp1d(new_x, new_y4, kind='cubic')
f_cubic5 = interp1d(new_x, new_y5, kind='cubic')
f_cubic6 = interp1d(new_x, new_y6, kind='cubic')
f_cubic7 = interp1d(new_x, new_y7, kind='cubic')
f_cubic8 = interp1d(new_x, new_y8, kind='cubic')
f_cubic9 = interp1d(new_x, new_y9, kind='cubic')




# Plot.
plt.plot(x, y1, '-o', label='data')
plt.plot(x, y2, '-o', label='data')
#plt.plot(xnew, f_linear(xnew), '-', label='linear',color = 'red')
plt.plot(xnew, f_cubic(xnew), '-', label='cubic',color = 'green')
plt.plot(xnew, f_cubic2(xnew), '-', label='cubic',color = 'red')
plt.plot(xnew, f_cubic3(xnew), '-', label='cubic',color = 'black')
plt.plot(xnew, f_cubic4(xnew), '-', label='cubic',color = 'yellow')
plt.plot(xnew, f_cubic5(xnew), '-', label='cubic',color = 'orange')
plt.plot(xnew, f_cubic6(xnew), '-', label='cubic',color = 'brown')
plt.plot(xnew, f_cubic7(xnew), '-', label='cubic',color = 'pink')
plt.plot(xnew, f_cubic8(xnew), '-', label='cubic',color = 'purple')
plt.legend(loc='best')
plt.show()





import numpy as np
import numpy as np
from scipy.interpolate import make_interp_spline
import matplotlib.pyplot as plt

''''# Dataset

y = y1

X_Y_Spline = make_interp_spline(x, y)

# Returns evenly spaced numbers
# over a specified interval.
X_ = np.linspace(x.min(), x.max(), 100)
Y_ = X_Y_Spline(X_)

# Plotting the Graph
plt.plot(x,y,color = 'black')
plt.plot(X_, Y_)
plt.title("Plot Smooth Curve Using the scipy.interpolate.make_interp_spline() Class")
plt.xlabel("X")
plt.ylabel("Y")
plt.show()'''