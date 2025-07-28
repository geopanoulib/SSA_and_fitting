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

"""
NAME: ImportingData
PURPOSE: Reads data from a csv file and saves in three arrays the 
         enumeration, the year and the dlod values.
PARAMS: - FileName:  name of the csv file 
        - delimeterSymbol: the delimiter of the csv file
RETURN: three arrays (num, year, dlod)
"""
def ImportingData(FileName, delimeterSymbol):
    data = np.genfromtxt(FileName ,delimiter=delimeterSymbol)    #Read the data from a csv
    # csv format: num;year;month;day;dlod;s_dlod
    
    N = len(data) # Total number of data
    
    dlod = []            # DLOD values 
    num = []             # Enumeration (1 to n)
    year = []            # Dates
    
    for idx,row in enumerate(data[:N]):    
        num.append(int(row[0]))
        year.append(int(row[1]))
        dlod.append(float(row[4]))    
        
    return num, year, dlod
    
"""
NAME: CalculationOfElementaryMatrices
PURPOSE: Calculates the elementary matrices X_i of trajectory matrix Χ.
PARAMS: - U: U matrix (L x L)
        - S:  Sigma vector (L)
        - V:  V matrix (K X K)
        - n_comp: numper of components
RETURN: X_elem 
"""   
def CalculationOfElementaryMatrices(U,S,V,n_comp):
    for i in range(0,n_comp + 1):
        outer_product = np.outer(U[:,i],V[:,i]) # Outer product of i-th column of U with i-th column of V 
        X_elem = S[i]*outer_product
        yield X_elem


"""
NAME: ConvertTrajectoryMatrixToTimeSeries
PURPOSE: Calculates the averages of the anti-diagonals of
        the elementary matrix X_i and returns a time series.
PARAMS: - X_i: elementary matrix
RETURN: TimeSeries
""" 
def ConvertTrajectoryMatrixToTimeSeries(X_i):
    X_rev = X_i[::-1]         # Reverse the column ordering of X_i because we need the anti-diagonals, 
                              #but python do calculations with diagonals
    antiDiagonalMeans = []  # Empty list to save the averages of the anti-diagonals.

    # A for loop for a number of indexes
    #that correspond to the position of the diagonals (the positions start from -k+1 to l-1, for matrix k x l)
    for i in range(-X_i.shape[0]+1, X_i.shape[1]):

        antiDiagonal = X_rev.diagonal(i) # Find the X_rev's diagonal in position i

        mean = antiDiagonal.mean()       # Calculate the average of the antidiagonals values

        antiDiagonalMeans.append(mean)  
    
    TimeSeries = np.array(antiDiagonalMeans)  # Convert the list to a numpy array and return it
    return TimeSeries
