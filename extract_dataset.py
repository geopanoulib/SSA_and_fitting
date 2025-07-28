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

import requests
import pandas as pd

'''______________________[SET PARAMETERS]______________________________'''

sUrl = "https://datacenter.iers.org/data/latestVersion/EOP_20_C04_one_file_1962-now.txt"

'''___________________________________________________________________'''

response = requests.get(sUrl)

sText = response.text

aLines = sText.splitlines()

aDLODdata = []
i = 1
for sLine in aLines:
    aData = sLine.split(" ")
    aData = list(filter(None, aData)) # Filter out empty values in the array
    
    if not aData or aData[0] == '#': # Skip comments and empty lines
        continue
    
    if  aData[0] == '2023' and  aData[1] == '1' and aData[2] == '2': # Stop loop we only want data till 2023-01-01
        break
        
    
    aDLODdata.append([i, aData[0], aData[1], aData[2],float(aData[12])*10**3,float(aData[20])*10**3])
    i +=1
    
df = pd.DataFrame(aDLODdata)

df.to_csv("DLOD_1962_2023.csv", mode='a', header=False, index=False, sep=';')


