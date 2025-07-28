# Singular Spectrum Analysis & Least-Squares fitting

## Code Information

This code contains six executable programs (three for the application of the SSA method and three for the Least-Squares fitting) and one additional file with the functions used from the SSA scripts. The name of the files are:

<u>SSA</u>

- extract_dataset.py

- Investigation.py

- Grouping.py

- functions.py
  
  

<u>Least-Squares fitting</u>

- Linear.py

- Harmonic.py

- Total_harmonic.py


The data where extracted from a given url and thus created a csv file composed by 6 columns seperated by semicolon. (e.g. enumeration;year;month;day;Δlod;weight).
The SSA programs takes as input this csv file where the method will be applied. Here is an example of the files format:
m
```csv
1;1962;1;1;1.723;1.4
2;1962;1;2;1.669;1.4
3;1962;1;3;1.582;1.4
4;1962;1;4;1.496;1.4
5;1962;1;5;1.416;1.4
```

The least-squares fitting programs take as input files the csv that obtained from the SSA program, one per periodicity for the Harmonic script and one for the trend for the Linear script.

## Set parameters

Before executing the script extract_dataset.py, you should set the url where the data will be extracted. 

Before executing the script Investigation.py, you should first set the name (CSVFileName) of your csv file of your data and the window lenght (L) of the trajectory matrix. Also, you could adjust the number of components that will be examined (n_comp). Note that it is better to run the script multiple times with different combinations of the parameters, to see which work best depending the problem and the nature of your data.

In the Grouping script you should set the groups of components that will form the resulting periodicities, based on the investigation you conducted with the findings of the Investigation.py script.

For the Harmonic script you should set the initial values of period (and perhaps amplitude and phase) for each periodicity and then execute it for each one by specifing each time its name (sPeriodic). Note that maybe you should run the script multiple times in order to find the best initial values for the harmonic fitting to coincide. Furthermore, before the least-squares execution do not forget to create a folder "LeastSquares" where the results will be saved.

All files (scripts and csv with data) should be in the same folder.

## Executing Scripts

To execute each script you should run in the linux command line the next command (if you are using python 3):

```bash
python3 Investigation.py
```

The extract_dataset script should be executed first then the Investigation and then the Grouping script. 
Then follows the Linear and Harmonic (for each periodicity) scripts and lastly the Total_harmonic.

