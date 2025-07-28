'''
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
'''


import numpy as np
import matplotlib.pyplot as plt
import time
import os
import shutil

from functions import ImportingData
from functions import CalculationOfElementaryMatrices
from functions import ConvertTrajectoryMatrixToTimeSeries

'''______________________[SET PARAMETERS]______________________________'''
sCSVFileName = "DLOD_1962_2023.csv"
sDelimiter = ";"

L = 11140  # Set the window length (L) of the trajectory matrix (2<=L<=N/2)

Ncomp = 100   # Number of components we will examine after the decomposition of trajectory matrix

'''___________________________________________________________________'''

# Create directories to save the plots
plotDir = "Plots" 
W_corrDir = "%s/w_corr"%plotDir
compDir = "%s/Elem_comps"%plotDir

os.path.exists(plotDir) and shutil.rmtree(plotDir) # If directory exists delete it
os.makedirs(plotDir, exist_ok=True)

os.path.exists(W_corrDir) and shutil.rmtree(W_corrDir)
os.makedirs(W_corrDir, exist_ok=True)

os.path.exists(compDir) and shutil.rmtree(compDir)
os.makedirs(compDir, exist_ok=True)

'''______________________[IMPORT DATA]______________________________'''

num, year, dlod = ImportingData(sCSVFileName, sDelimiter) #Read the data from a csv
N = len(num)

if (L>N/2 or L<2):   
    raise ValueError(f'Warning! Give a number from 2 to {N//2}')
K = N - L + 1  # Dimensions of trajectory matrix


'''_______________________[PLOT SETTINGS]______________________________________'''

plt.rcParams['figure.figsize'] = (10,8)    
plt.rcParams['font.size'] = 14             
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['image.cmap'] = 'plasma'

from cycler import cycler
cols = cycler(color=plt.get_cmap('tab20').colors)


'''____________________[PLOT DATA]____________________________________________'''


plt.plot(dlod, color="blue", lw=1.5)
plt.xticks(num[::365*5],year[::365*5], rotation=45)

plt.xlabel('Time(Year)')
plt.ylabel('ΔLOD(ms)')
plt.title('Excess of the length of day')
plt.savefig(f"{plotDir}/dlod_{N}.svg")
plt.close()


'''_________________[CREATE TRAJECTORY MATRIX]____________________________'''

# Create the trajectory matrix (L x K)
columns = [] 
for i in range(0,K):
    columns.append(dlod[i:i+L])
    
X = np.column_stack(columns)  

# Check that the trajectory matix has the right dimensions
c, r = np.shape(X)
if (c != L or r != K):
        raise ValueError("Dimensions of the trajectory matrix are incorrect")


'''__________________[PLOT TRAJECTORY MATRIX]____________________________'''                                        

plt.rcParams['figure.figsize'] = (10,8)    
ax = plt.matshow(X)
plt.xlabel("$L$-Lagged Vectors")
plt.ylabel("$K$-Lagged Vectors")
plt.colorbar(ax.colorbar, fraction=0.045)
ax.colorbar.set_label("DLOD(ms)")
plt.title("The Trajectory Matrix for our Time Series")
plt.savefig(f"{plotDir}/Trajectory.svg")
plt.close()

d = np.linalg.matrix_rank(X) 
print("\n\tThe intrinsic dimensionality of the trajectory space is: ",d)


'''__________________________[SVD]__________________________________________'''

#Singular Value Decomposition (SVD)
start = time.time()

U, Sigma, VT = np.linalg.svd(X)  # U (L x L), Sigma (L), V (K X K)
V = np.transpose(VT)             # Numpy calculates V transpose but we want V

end = time.time()

t = end - start    # Duration of SVD calculation in seconds

print("\n\tSVD: %.2f hours, %.2f minutes, %.4f seconds" %(t//3600, (t%3600)//60, t%60))


'''_______________[CUMULATIVE & RELATIVE CONTRIBUTION PLOTS]___________________'''

fig, ax = plt.subplots(1, 2, figsize=(14,5))

sigma_sumsq = (Sigma**2).sum()

# Plot of Cumulative Contribution
cumulative = (Sigma**2).cumsum() / sigma_sumsq * 100
ax[0].plot(cumulative, lw=2.5)
ax[0].set_xlim(0,Ncomp)        
ax[0].set_title("Cumulative Contribution of $\mathbf{X}_i$ to Trajectory Matrix")
ax[0].set_xlabel("$i$")
ax[0].set_ylabel("Contribution (%)")

# Plot of Relative Contribution
relative = Sigma**2 / sigma_sumsq * 100
ax[1].plot(relative, lw=2.5)
ax[1].set_xlim(0,Ncomp)
ax[1].set_title("Relative Contribution of $\mathbf{X}_i$ to Trajectory Matrix")
ax[1].set_xlabel("$i$")
ax[1].set_ylabel("Contribution (%)")
plt.savefig(f"{plotDir}/Cumulative_Relative_Con_{Ncomp}.svg")
plt.close()


'''____________________[DIAGONAL AVERAGING]_______________________________'''

# Diagonal averaging step for the creation of the elementary time series components
elem_matrices = CalculationOfElementaryMatrices(U,Sigma,V,Ncomp)
F_elem = []         # Empty array for the elementary time series
start = time.time()

for i,X_elem in enumerate(elem_matrices):
    ts = ConvertTrajectoryMatrixToTimeSeries(X_elem)     # Convert to time series with function X_to_TS
    F_elem.append(ts)        # Save all the elementary time series, as columns, to the F_elem array

    
end = time.time()

t = end - start   # In seconds

print("\n\tElementary: %.2f hours, %.2f minutes, %.4f seconds" %(t//3600, (t%3600)//60, t%60))


# Save elementary components in a binary file
with open("elementary.npy", 'wb') as f:
        np.save(f, F_elem)


'''___________________[W-CORRELATION MATRIX]________________________________'''

# Calculation of weights w. Depending on how many times an element appears 
#in the trajectory matrix Χ, which is relevant to its position
w = np.empty(N,dtype=int)
for i in range(N):
    if i<L:
        w[i] = i+1
    elif i<K:
         w[i] = L
    else:
         w[i] = N-i



start = time.time()

F_wnorms = np.array([w.dot(F_elem[i]**2) for i in range(Ncomp)]) # Inner product of w and F_elem[i]
F_wnorms = F_wnorms**(-0.5)      # Calculate the reversed square root, because we will need it next

# Calculation of W-correlation matrix
# Iterate all pairs of i's and j's (i != j). Noting that Wij = Wji.
Wcorr = np.identity(Ncomp+1) # Define matrix Wcorr with dimensions (Ncomp+1) x (Ncomp+1)
for i in range(Ncomp):
    for j in range(i+1,Ncomp):
        Wcorr[i,j] = abs(w.dot(F_elem[i]*F_elem[j]) * F_wnorms[i] * F_wnorms[j])
        Wcorr[j,i] = Wcorr[i,j] 
        
end = time.time()

t = end - start # In seconds

print("\n\tW-correlation matrix: %.2f hours, %.2f minutes, %.4f seconds" %(t//3600, (t%3600)//60, t%60))

'''____________________________[PLOT W-CORR MATRIX]______________________________________'''

#Plot W-correlation matrix by 10 values
for i in range(0, Ncomp+1-10, 10):
    n1 = i
    n2 = i + 10
    
    ax = plt.imshow(Wcorr)
    plt.colorbar(ax.colorbar, fraction=0.045)
    ax.colorbar.set_label("$W_{ij}$")
    plt.xlabel(r"$\tilde{F}_i$")
    plt.ylabel(r"$\tilde{F}_j$")
    plt.xlim(n1-0.5, n2+0.5)
    plt.ylim(n2+0.5, n1-0.5)
    plt.clim(0,1)
    plt.title(f"W-Correlation for Components {n1}-{n2}")
    plt.savefig(f'{W_corrDir}/W-corr({n1}_{n2}).svg')
    plt.close()


'''____________________________[PLOT ELEMENTARY COMPONENTS]______________________________________'''
plt.rcParams['axes.prop_cycle']

for i in range(0,Ncomp+1-4,4):
    n1 = i
    n2 = min(i+4, d) # In case of noiseless time series, where d < n2.
    for i in range(n1,n2+1):
        plt.plot(F_elem[i], lw=2)         # Create plot with elementary components
    
    #plt.plot(dlod, alpha=1, lw=1)         # In the same plot the original time series
    plt.xticks(num[::365*5],year[::365*5], rotation=45)  

    plt.xlabel("Time(Year)") 
    plt.ylabel("ΔLOD(ms)")  
    legend = [r"$\tilde{F}_{%s}$" %i for i in range(n1,n2+1)] + ["$F$"] 
    plt.title(f"Components {n1} to {n2} of our Time Series")      
    plt.legend(legend, loc=(1.01,0.4))
    plt.savefig(f"{compDir}/Components({n1}_{n2}).svg")
    plt.close()


        

        

