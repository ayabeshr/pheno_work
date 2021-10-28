import os, sys
from pathlib import Path
import pandas as pd
import csv 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# Define CSV file 
df = pd.read_csv("/home/aya/Pheno_Work/analysis/data_plots/Xsection_Test_20211020.csv")

# Define coloumns in csv file in terms of x, y and z for plotting
x = df['csqr_total']
y = df['mbreaking']
z = df['new_xsc(pb)']

# define figure for plotting 
fig = plt.figure(num=None, figsize=(6,6), dpi=150, facecolor='white', edgecolor='white', frameon=False)
my3dplot = fig.add_subplot(111, projection = '3d')

# plot configuration
my3dplot.scatter(x, y, z, marker = "o", color = "blue") 
plt.xlabel(r"$\chi^2$", color="orange")
plt.ylabel(r"$M_{12}$"+" x"+r"$10^{3}$", color="green")
my3dplot.set_zlabel(r"$\sigma(gg\rightarrow$"+"H"+r"$\rightarrow$"+"hh"+r"$\rightarrow$"+"bb4Mu) x"+r"$10^{-8}$"+" (pb)", color="red")
my3dplot.view_init(27)
plt.savefig("/home/aya/Pheno_Work/analysis/data_plots/newplots/3dplot.png")
