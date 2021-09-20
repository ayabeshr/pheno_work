import os
import openpyxl
import sys
from pathlib import Path
import pandas as pd
import csv 
from matplotlib import pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D

# Define CSV file 
df = pd.read_csv("/home/aya/Desktop/Pheno_Work/analysis/Xsection_Test_20210817_2.csv")

# Define coloumns in csv file in terms of x, y and z for plotting
x = df['csqr_total']
y = df['mbreaking']
z = df['xsection(pb)']

# plot data from file 

# Define a figure for plotting
fig = plt.figure(figsize=(7,7))
my3dplot = fig.add_subplot(111, projection = '3d')

my3dplot.scatter(x, y, z, marker = "o", color = "blue") 
plt.xlabel(r"$\chi^2$", color="blue")
plt.ylabel(r"$M12 x 10^2$", color="blue")
my3dplot.set_zlabel(r"$\sigma(pb)(gg \rightarrow H \rightarrow hh \rightarrow bb4Mu) x 10^{-11}$", color="red")
my3dplot.view_init(27)

# way to add a legend
#my3dplot.legend(loc='upper left')    # the contents of legend are define in label 
