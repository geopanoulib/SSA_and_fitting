'''
**
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

import numpy as np
import matplotlib.pyplot as plt

resultsDir = 'LeastSquares' # Folder for saving the resulting plots of the least square fitting

'''______________________[SET PARAMETERS]______________________________'''
sCSVFileName = "DLOD_1962_2023.csv"
sDelimiter = ";"


'''____________________[INITIAL VALUES]_______________________________________'''

aInitialValues = np.genfromtxt(f"{resultsDir}/results_harmonic.csv" ,delimiter = sDelimiter)    #Read the data from a csv
 # csv format: Amptitude;Period;Initial phase
 
k = len(aInitialValues) #  number of periodicities
m = 3 * k + 2  # number of the unknown parameters

A0 = aInitialValues[:k, 0]            # Amptitude 
T0 = aInitialValues[:k, 1]             # Period
f0 = aInitialValues[:k, 2]            # Initial phase

s0 = -0.0001194
c0 = 3.0477

'''__________________________________[IMPORT DATA]____________________________'''

data = np.genfromtxt(sCSVFileName, delimiter = sDelimiter)

x = data[:, 0] 
y = data[:, 4]
year = data[:, 1].astype(int)

n = len(x)

'''___________________________[CALCULATIONS]_______________________________________'''

dx = np.ones(m)
threshold = 10**(-4)
it = 0
A = np.ones((n, m))
A[:, m - 2] = x
dL = np.zeros(n)
N = np.zeros((m, m))
u = np.zeros(m)


while (sum(abs(dx)) > threshold and it < 10):
    #print('\n\titeration : {:}'.format(it))
    w = 2 * np.pi / T0
    dL = y - s0 * x - c0
    for i in range(0, k):
    	A[:, i] = np.cos(w[i] * x + f0[i])
    	A[:, i + k] = 2 * np.pi * x * A0[i] * np.sin(w[i] * x + f0[i]) / (T0[i]**2)
    	A[:, i + 2 * k] = -A0[i] * np.sin(w[i] * x + f0[i])
    	dL -= A0[i] * np.cos(w[i] * x + f0[i])
    N = np.dot(np.transpose(A), A)
    u = np.dot(np.transpose(A), dL)
    dx = np.linalg.inv(N)@u
    for i in range(0, k):
        A0[i] += dx[i]
        T0[i] += dx[k + i]
        f0[i] += dx[2 * k + i]
    s0 += dx[m - 2]
    c0 += dx[m - 1]
    it += 1
sumdli2 = np.dot(np.transpose(dL), dL)    
s0_apo = np.sqrt(sumdli2 / (n - m))

Vx = np.linalg.inv(N) * s0_apo**2
sA = np.zeros(k)
sT = np.zeros(k)
sf = np.zeros(k)


for i in range(0, k):
    sA[i] = np.sqrt(Vx[i][i])
    sT[i] = np.sqrt(Vx[k + i][k + i])
    sf[i] = np.sqrt(Vx[2 * k + i][2 * k + i])


s_s = np.sqrt(Vx[m - 2][m - 2])
sc = np.sqrt(Vx[m - 1][m - 1])


'''___________________________[PRINT RESULTS]_________________________________'''

print('n = %d points\t m = %d parameters' % (n, m))
print('iterations =', it)
for i in range(0, k):
    print('A{} = {:.5f} ms\t+/- {:.6f}'.format(i, A0[i], sA[i]))
    print('T{} = {:.2f} days\t+/- {:.3f}'.format(i, T0[i], sT[i]))
    print('f{} = {:.5f} rad\t+/- {:.6f}'.format(i, f0[i], sf[i]))
    print(50 * "-")
print('s = {:.5f} ms/days\t+/- {:.5f} ms/days\n\tc = {:.5f} ms\t+/- {:.5f} ms'.format(s0, s_s, c0, sc))
print('s0_aposteriori = +/- {:.5f} ms\n'.format(s0_apo))

x_opt = np.arange(0,max(x),0.2)

y_opt = np.zeros(len(x_opt))

'''_________________________________[PLOTS]________________________________'''

for i in range(0,len(x_opt)):
    for j in range(0,k):
        y_opt[i] += A0[j]*np.cos(2*np.pi*x_opt[i]/T0[j]+f0[j])
    y_opt[i] += s0*x_opt[i]+c0


plt.rcParams['figure.figsize'] = (10,8) 
plt.xticks(x[::365*5], year[::365*5], rotation=45)
plt.plot(x, y, c = "b", linewidth = 0.5, label = 'SSA') 
plt.plot(x_opt, y_opt, c = 'Red', ls = '-', lw = '1', label = 'best-fit harmonic')
plt.grid(c = 'Gray', ls = '--', lw = 0.4)
plt.title('Ηarmonic fitting of Total Series')
plt.legend()
plt.xlabel('Time(day)')
plt.ylabel('ΔLOD(ms)')
plt.savefig(f"{resultsDir}/TotalSeries.svg")
#plt.show()
plt.close()
