'''
**
 *
 * Copyright (c) 2025, Georgios Panou
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 *
 * Authors: Pappa A. and Panou G. <geopanou@survey.ntua.gr>
 *
'''

import matplotlib.pyplot as plt
import numpy as np
import os
import shutil


" The linear mathimatical model we based on: y = ax + b "


resultsDir = 'LeastSquares' # Folder for saving the resulting plots of the least square fitting

os.path.exists(resultsDir) and shutil.rmtree(resultsDir) # If directory exists delete it
os.makedirs(resultsDir, exist_ok=True)

'''__________________________________[IMPORT DATA]____________________________'''

FileName = "CSV/trend.csv"

data = np.loadtxt(FileName, delimiter = ';')
x = data[:, 0]  
y = data[:, 1]


n = len(x)              # number of data
m = 2                   # number of the unknown parameters
r = n - m               # degrees of freedom

'''___________________________[CALCULATIONS]_______________________________________'''

A = np.zeros((n, m), dtype = np.float64)        # initialization of the design matrix Α

A[:, 0] = np.array([x])         # calculation of first column of Α
A[:, 1] = np.array(np.ones(n))  # calculation of second column of Α

At = A.T            # transpose matrix A

l = y               # vector of data

N = np.matmul(At, A)        # calculation of Ν matrix

Ninv = np.linalg.inv(N)             # inversion of matrix Ν

u = np.matmul(At, l)            # calculation of u vector

X = np.matmul(Ninv, u)          # solution matrix

v = np.matmul(A, X) - l              # calculation of the vector of residuals

s0 = np.sqrt(np.matmul(v.T, v) / r)         # calculation of the a-posteriori variance factor 

Vx = s0 * s0 * Ninv         # Calculation of the a-posteriori variance-covariance matrix of the vector X ̂_1

s = X[0]
c = X[1]

# calculation of the errors of the parameters
s_s = np.sqrt(Vx[0][0])
s_c = np.sqrt(Vx[1][1])

" Calculation of the best-fit line"
x_opt = x
y_opt = X[0] * x_opt + X[1]

'''___________________________[PRINT RESULTS]_________________________________'''

print(f'\n s = {s:.7f}  +/- {s_s:.10f} ms/day')
print(f' c = {c:.4f}  +/- {s_c:.6f} ms')
print(f"\n\tσ0 = +/- {s0:.4f} ms")


'''_________________________________[PLOTS]________________________________'''

plt.rcParams['figure.figsize'] = (10,8) 
plt.plot(x, y, c = "b", linewidth = 0.5, label = 'SSA') 
plt.plot(x_opt, y_opt, c = 'Red', ls = '-', lw = '1', label = 'best-fit line')
plt.grid(c = 'Gray', ls = '--', lw = 0.4)
plt.title('Linear fitting of trend')
plt.legend(loc='upper right')
plt.xlabel('Time(day)')
plt.ylabel('ΔLOD(ms)')
plt.savefig(f"{resultsDir}/trend.svg")
plt.show()