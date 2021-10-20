import os, sys
import openpyxl
import array as arr 
from pathlib import Path
import pandas as pd
import csv 
from matplotlib import pyplot as plt 


## Define CSV file 
df = pd.read_csv("/home/aya/Desktop/Pheno_Work/analysis/data_plots/Xsection_Test_20211020.csv")
x = df['mbreaking'].values
y = df['new_xsc(pb)'].values
xarr = arr.array('d',[-3.8, -3.2])
yarr = arr.array('d',y)

# Xsection Vs M12
xsc_Vsm12 = df.plot(x = "mbreaking", y = "new_xsc(pb)", 
                    marker = ".", color = "blue",
                    label =  "$\sigma (pb)$ \n at tree-level", 
                    linestyle="")
#threshold = 3

plt.xlabel(r"$M_{12}^{3}$")    # r refers to raw string, placed before string, used when writting math expressions
plt.ylabel(r"$\sigma(gg \rightarrow H \rightarrow hh \rightarrow bb4Mu) x 10^{-8} (pb)$")
plt.grid(True)
plt.legend(loc='upper left')
#plt.axhline(threshold, color = "0", alpha = 1, lw=2, linestyle = "--")
#plt.fill_between(-y,-0.1,4, where=y > threshold, color='tab:gray', alpha=0.3)

# if you need to plot certain range of a coloumn values, for eg, 1st 10 values do: 
# x = df['coloumn1_name'].iloc[0:10].values

# Xsection Vs Total Chi Square
xsc_Vschi = df.plot(x = "csqr_total", y = "new_xsc(pb)", 
                    marker = ".", color = "red", alpha = 0.5,
                    label =  "$\sigma (pb)$ \n at tree-level",
                    linestyle="",  fillstyle=None, )

plt.xlabel(r"$\chi^{2}$")
plt.ylabel(r"$\sigma(gg \rightarrow H \rightarrow hh \rightarrow bb4Mu) x 10^{-8} (pb)$")
plt.grid(True)
plt.legend(loc="upper left")

# H2 mass Vs M12 
h2mass_Vsm12 = df.plot(x = "mbreaking", y = "mH(GeV)", 
                       marker = ".", color = "mediumseagreen", 
                       linestyle="")

plt.xlabel(r"$M_{12}^{3}$")
plt.ylabel("H_mass(GeV)")
plt.grid(True)

# Xsection Vs H2 mass 
xsc_Vsh2mass = df.plot(x = "mH(GeV)", y = "new_xsc(pb)", 
                       marker = ".", color = "magenta", 
                       label =  "$\sigma (pb)$ \n at tree-level", 
                       linestyle="")

plt.xlabel("H_mass(GeV)")
plt.ylabel(r"$\sigma(gg \rightarrow H \rightarrow hh \rightarrow bb4Mu) x 10^{-8} (pb)$")
plt.grid(True)
plt.legend(loc='upper right')
