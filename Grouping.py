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
import csv
import os
import shutil

from functions import ImportingData

'''______________________[SET PARAMETERS]______________________________'''
sCSVFileName = "DLOD_1962_2023.csv"
sDelimiter = ";"
'''___________________________________________________________________'''

csvDir = 'CSV' # Folder for saving the resulting csv 
os.path.exists(csvDir) and shutil.rmtree(csvDir) # If directory exists delete it
os.makedirs(csvDir, exist_ok=True)

plotsDir = 'periodicities' # Folder for saving the plots of periodicities
os.path.exists(plotsDir) and shutil.rmtree(plotsDir)
os.makedirs(plotsDir, exist_ok=True)

'''_____________________[PLOT SETTINGS]__________________________________'''

plt.rcParams['figure.figsize'] = (10,8)    
plt.rcParams['font.size'] = 14             
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['image.cmap'] = 'plasma'

from cycler import cycler
cols = cycler(color=plt.get_cmap('tab20').colors)


'''______________________[IMPORT DATA]_________________________________'''

with open('elementary.npy', 'rb') as f:
        F_elem = np.load(f)
        
num, year, dlod = ImportingData(sCSVFileName, sDelimiter) #Read the data from a csv
N = len(num)

'''__________________________[GROUPING]____________________________________'''

trend = F_elem[0]
periodic1 = F_elem[1] + F_elem[2]
periodic2 = F_elem[3] + F_elem[4]
periodic3 = F_elem[5] + F_elem[6]
periodic4 = F_elem[7] + F_elem[8]
periodic5 = F_elem[10] + F_elem[11] 
periodic6 = F_elem[12] + F_elem[13]
periodic7 = F_elem[14] + F_elem[15]
periodic8 = F_elem[18] + F_elem[19]
periodic9= F_elem[21] + F_elem[22]
periodic10 = F_elem[26] + F_elem[27]
periodic11 = F_elem[28] + F_elem[29] 
periodic12 = F_elem[32] + F_elem[33]
periodic13 = F_elem[34] + F_elem[35]
periodic14 = F_elem[36] + F_elem[37]
periodic15 = F_elem[43] + F_elem[44]
periodic16 = F_elem[49] + F_elem[50]
periodic17 = F_elem[54] + F_elem[55]
periodic18 = F_elem[56] + F_elem[57]
periodic19 = F_elem[60] + F_elem[61]
periodic20 = F_elem[62] + F_elem[63]
#periodic21 = F_elem[64] + F_elem[65]
#periodic22 = F_elem[66] + F_elem[67]
#periodic23 = F_elem[70] + F_elem[71]


components = [periodic1, periodic2, periodic3, periodic4, periodic5, periodic6,              
              periodic7, periodic8, periodic9, periodic10,periodic11,periodic12,
              periodic13, periodic14, periodic15, periodic16, periodic17, periodic18,
              periodic19, periodic20]

uNumberOfPeriodics = len (components)

noise = dlod - trend
for periodic in components:
    noise -= periodic


'''____________________________[SAVE IN CSV]____________________________________'''

# Save values of each periodicity  in a different csv for leasts-square fitting 
trend_list = []
for n,i in enumerate(trend):
        trend_list.append([n,i])
        
with open(F"{csvDir}/trend.csv", 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile,delimiter=";")
        csvwriter.writerows(trend_list)
csvfile.close() 

for j in range(uNumberOfPeriodics):   # Set numbers of periodicities
    data=[]
    for n,i in enumerate(components[j]):
        data.append([n,i])                # Put enumeration to the data 
    
    with open(f"{csvDir}/periodic{j+1}.csv", 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile,delimiter=";")
        csvwriter.writerows(data)
    csvfile.close()



'''____________________________[PLOTS]________________________________'''

# Plot grouped components together 
plt.plot(dlod, lw=1)
plt.plot(trend)
for i in range(uNumberOfPeriodics):
    plt.plot(components[i])
plt.plot(noise, alpha=0.5)


plt.xticks(num[::365*5],year[::365*5], rotation=45) 
plt.xlabel("Time(Year)")
plt.ylabel("ΔLOD(ms)")
groups = ["trend"]+[ f"periodic {x+1}" for x in range(uNumberOfPeriodics)]+["noise"]


legend = ["original"] + [group for group in groups]
plt.legend(legend, loc=(1.01,0.0))
plt.title("Grouped Time Series Components")
plt.savefig(f"{plotsDir}/Grouped.svg")
plt.show()
plt.close()



# Plot each grouped component seperatelly
#plt.rcParams['figure.figsize'] = (10,3)
plt.plot(trend)
plt.xticks(num[::365*5],year[::365*5], rotation=45)
plt.xlabel("Time(Year)")
plt.ylabel("ΔLOD(ms)")
plt.title("Trend",fontsize=16)
plt.tight_layout()
plt.savefig(f'{plotsDir}/trend.svg')
#plt.show()
plt.close()


for i in range(uNumberOfPeriodics):
    plt.plot(components[i])
    plt.xticks(num[::365*5],year[::365*5], rotation=45)
    plt.xlabel("Time(Year)")
    plt.ylabel("ΔLOD(ms)")
    plt.title(f"Periodic {i+1}", fontsize=16)
    plt.tight_layout()
    plt.savefig(f'{plotsDir}/periodic{i+1}.svg')
    #plt.show()
    plt.close()

plt.plot(noise)
plt.xticks(num[::365*5],year[::365*5], rotation=45)
plt.xlabel("Time(Year)")
plt.ylabel("ΔLOD(ms)")
plt.title("Noise",fontsize=16)
plt.tight_layout()
plt.savefig(f'{plotsDir}/noise.svg')
#plt.show()
plt.close()

