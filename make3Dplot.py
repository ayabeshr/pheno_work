import os
import openpyxl
import sys
from pathlib import Path
import pandas as pd
import csv 
from matplotlib import pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D


# Define CSV file 
df = pd.read_csv("/home/aya/Videos/Book2.csv")

# Define coloumns in csv file in terms of x, y and z for plotting
x = df['Chi']
y = df['param']
z = df['X-section(pb)']

# plot data from file 

# Define a figure for plotting
fig = plt.figure(figsize=(7,7))
my3dplot = fig.add_subplot(111, projection = '3d')

my3dplot.scatter(x, y, z, marker = "o", color = "blue", label = "2HDM-type1\nPhase2 Simulation\n14TeV")     # or you can use my3dplot.scatter(x,y,z) which is equivalent and blue color is default

plt.xlabel("chi.e+02", color = "red")          # or you can equivalently use my3dplot.set_xlabel('chi.e+02', color = "red")
plt.ylabel("m12.e+02", color = "magenta")      # or you can equivalently use my3dplot.set_ylabel('m12.e+02', color = "magneta")
my3dplot.set_zlabel('xsection(pb).e-07', color = "green")

# way to add a legend
my3dplot.legend(loc='upper left')    # the contents of legend are define in label 
