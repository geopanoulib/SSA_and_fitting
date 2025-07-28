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
import pandas as pd


" The mathematical model of the cosine wave we based on: y = A * cos((2 * pi * ti) / T + f) + c "

resultsDir = 'LeastSquares' # Folder for saving the resulting plots of the least square fitting


'''______________________[SET PARAMETERS]____________________________'''

sPeriodic = "periodic1"
sFilePath = "CSV/"

#dictionary with initial values for the period (in some difficult cases you need to set also the ampitude or/and the initial phase)
T0_dict = {'periodic1' : 6750, 'periodic2' : 365, 'periodic3' : 13.66, 'periodic4' : 182.5,
           'periodic5' : 4015, 'periodic6' : 27.5534, 'periodic7' : 13.63, 'periodic8' : 880,
           'periodic9' : 9.1333, 'periodic10' : 2000/1.5, 'periodic11' : 791.6303, 'periodic12' : 121.9517346120, 
           'periodic13' : 609.5531769783, 'periodic14' : 31.8129793394, 'periodic15' : 212.7552901072,
           'periodic16' : 260.8658264003, 'periodic17' : 91.3007758958, 'periodic18' : 191.3170500930, 
           'periodic19' : 9.1206057533,  'periodic20' : 69.6606879610
}


A0_dict = {'periodic1' : 'fix', 'periodic2' : 0.38, 'periodic3' : 'fix', 'periodic4' : 'fix',
          'periodic5' : 0.23, 'periodic6' : 'fix', 'periodic7' : 'fix', 'periodic8' : 'fix',
          'periodic9' : 'fix',  'periodic10' : 'fix', 'periodic11' :'fix', 'periodic12' : 0.03,
          'periodic13' : 0.0309109476, 'periodic14' : 0.0329963246, 'periodic15' :0.0307932557,
          'periodic16' : 0.0169685712, 'periodic17' : 0.0179248216, 'periodic18' : 0.0111759304, 
          'periodic19' : 0.0258861665, 'periodic20' : 0.0183129189, 'periodic21' : 0.0057514996, 
          'periodic22' : 0.0185415377,'periodic23' : 0.0104402646
}
f0_dict = {'periodic1' : 'fix', 'periodic2' :'fix', 'periodic3' : 'fix', 'periodic4' : 'fix',
          'periodic5' : 'fix', 'periodic6' :'fix', 'periodic7' : 1, 'periodic8' : 1,
          'periodic9' : 'fix', 'periodic10' : 1,'periodic11' : 'fix','periodic12' : 1,
          'periodic13': 4.5156383300 , 'periodic14' : -0.1352860593, 'periodic15' : 5.0758826571,
          'periodic16' : 12.9581578543, 'periodic17' : -0.0766351059, 'periodic18' : 1.0384022812,
          'periodic19' : 18.2340573946, 'periodic20' : -0.3238826432
}

'''__________________________________[IMPORT DATA]____________________________'''

data = np.genfromtxt('%s%s.csv'%(sFilePath, sPeriodic), delimiter = ';') # import data from file into an array
td = data[:, 0] # days
yd = data[:, 1] # milliseconds
xmin = 0
xmax = len(yd)
step = 1
t = td[xmin:xmax:step]
y = yd[xmin:xmax:step]

n = len(t)  # number of data
m = 4       # number of the unknown parameters
r = n - m   # degrees of freedom


'''_____________[INITIAL VALUES OF UNKOWN PARAMETERS]_____________________________'''

if (A0_dict[sPeriodic] == 'fix'):
    A0 = (np.max(y) - np.min(y)) / 2            # Amptitude
else:
    A0 = A0_dict[sPeriodic]

T0 = T0_dict[sPeriodic]                     # Period

if (f0_dict[sPeriodic] == 'fix'):
   f0 = np.arccos((y[0] - A0 - np.min(y)) / np.abs(A0))      # Initial phase
else:
    f0 = f0_dict[sPeriodic]

c0 = A0 + np.min(y)                         # Offest along the y-axis

MaximumIterations = 20     
threshold = 10**(-4) 


'''__________________________[CALCULATIONS]____________________________________'''

iterations = 0
s0_i_1 = 1
s0_i = 2
J0 = np.zeros((n, m), dtype = np.float64)      # Initialization of the design matrix J0
while ( (np.abs(s0_i - s0_i_1) > threshold) and (iterations < MaximumIterations) ):
	s0_i_1 = s0_i
	
    # calculation of each column of J0
	J0[:, 0] = np.cos(2 * np.pi * t / T0 + f0) 
	J0[:, 1] = 2 * np.pi * A0 * t * np.sin(2 * np.pi * t / T0 + f0) / (T0**2)
	J0[:, 2] = -A0 * np.sin(2 * np.pi * t / T0 + f0)
	J0[:, 3] = 1.0
    
	J0t = J0.T     # transpose J0 matrix 
	
	F = A0 * np.cos(2 * np.pi * t / T0 + f0) + c0
	
	dl = y - F
	
	N = np.matmul(J0t, J0)         # Ν matrix
	Ninv = np.linalg.inv(N)        # Inversion of Ν matrix
	
	u = np.matmul(J0t, dl)         # vector u
	
	dx = np.matmul(Ninv, u)        # corrections vector dx
	A0 += dx[0]
	T0 += dx[1]
	f0 += dx[2]
	c0 += dx[3]
	
	v = np.matmul(J0, dx) - dl                # residuals vector
	s0_i = np.sqrt(np.matmul(v.T, v) / r)     # the a-posteriori variance factor
	
	iterations += 1


Vx = s0_i * s0_i * Ninv         # The a-posteriori variance-covariance matrix of the inverse of solution vector X ̂_1

# The errors of the parameters
sA = np.sqrt(Vx[0][0])
sT = np.sqrt(Vx[1][1])
sf = np.sqrt(Vx[2][2])
sc = np.sqrt(Vx[3][3])

# Calculation of the best fit harmonic 
t_opt = np.linspace(min(t), max(t), num = 20000)
y_opt = A0 * np.cos(2 * np.pi * t_opt / T0 + f0) + c0
	
'''_____________________________[PRINT RESULTS]___________________________________'''

print('\n n = %d points' % n)
print('\n Iterations = %d' % iterations)
print('\n A = %.10f  +/- %.10f ms' %(abs(A0),sA))
print(' T = %.10f  +/- %.10f days' %(T0,sT))
print('\n phi = %.10f  +/- %.10f rad' %(f0,sf))
print(' c = %.10f   +/- %.10f ms' %(c0,sc))
print('\n σ0 = +/- %.10f ms' % s0_i)

'''_____________________________[PRINT RESULTS]___________________________________'''
# write in csv file the results

# Combine values into a DataFrame
df = pd.DataFrame([[abs(A0), T0, f0]])

# Append to the CSV file
df.to_csv(F"{resultsDir}/results_harmonic.csv", mode='a', header=False, index=False, sep=';')

'''_________________________________[PLOTS]___________________________________'''

plt.rcParams['figure.figsize'] = (10,8) 
plt.plot(t, y, c = "b", linewidth = 0.5, label = 'SSA') 
plt.plot(t_opt, y_opt, c = 'Red', ls = '-', lw = '1', label = 'best-fit harmonic')
plt.grid(c = 'Gray', ls = '--', lw = 0.4)
plt.title('Ηarmonic fitting of %s' %sPeriodic)
plt.legend()
plt.xlabel('Time(day)')
plt.ylabel('ΔLOD(ms)')
plt.savefig(f"{resultsDir}/%s.svg" % sPeriodic)
plt.show()
plt.close()
