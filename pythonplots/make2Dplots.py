import os
import openpyxl
import sys
from pathlib import Path
import pandas as pd
import csv 
from matplotlib import pyplot as plt 


# Define CSV file 
df = pd.read_csv("/home/aya/Desktop/Pheno_Work/analysis/data_plots/Xsection_Test_20211019.csv")
#df = pd.read_csv("/home/aya/Downloads/Xsection_Test_20210719 - Sheet1(3).csv")
# plot data from file 
#m12_xsc = dfl.plot("Mvalue", "X-section(pb)", marker = "o", color = "red", linestyle = "")

# Xsection Vs M12
xsc_Vsm12 = df.plot(x = "mbreaking", y = "new_xsc(pb)", marker = "o", color = "blue", label =  "$\sigma (pb)$ \n at tree-level", linestyle="")
plt.grid(True)
plt.xlabel(r"$M_{12}^{3}$")
plt.ylabel(r"$\sigma(pb)(gg \rightarrow H \rightarrow hh \rightarrow bb4Mu) x 10^{-8}$")
plt.legend(loc='upper left')

# if you need to plot certain range of a coloumn values, for eg, 1st 10 values do: 
# x = df['coloumn1_name'].iloc[0:10].values

# Xsection Vs Total Chi Square
xsc_Vschi = df.plot(x = "csqr_total", y = "new_xsc(pb)", marker = "o", color = "red", label =  "$\sigma (pb)$ \n at tree-level", linestyle="")
plt.grid(True)
plt.xlabel(r"$\chi^{2}$")
plt.ylabel(r"$\sigma(pb)(gg \rightarrow H \rightarrow hh \rightarrow bb4Mu) x 10^{-8}$")
plt.legend(loc="upper left")

# H2 mass Vs M12 
h2mass_Vsm12 = df.plot(x = "mbreaking", y = "mH(GeV)", marker = "o", color = "green", linestyle="")
plt.grid(True)
plt.xlabel(r"$M_{12}^{3}$")
plt.ylabel("H_mass(GeV)")

# Xsection Vs H2 mass 
xsc_Vsh2mass = df.plot(x = "mH(GeV)", y = "new_xsc(pb)", marker = "o", color = "magenta", label =  "$\sigma (pb)$ \n at tree-level", linestyle="")
plt.grid(True)
plt.xlabel("H_mass(GeV)")
plt.ylabel(r"$\sigma(pb)(gg \rightarrow H \rightarrow hh \rightarrow bb4Mu) x 10^{-8}$")
plt.legend(loc='upper right')
