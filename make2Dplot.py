import os
import openpyxl
import sys
from pathlib import Path
import pandas as pd
import csv 
from matplotlib import pyplot as plt 


# Define CSV file 
df = pd.read_csv("/home/aya/Videos/Book2.csv")

# plot data from file 
# Plot-1 
xsc_Vsm12 = df.plot(x = "param", y = "X-section(pb)", marker = "o", color = "blue", linestyle="")
plt.grid(True)
plt.xlabel("m12.e+02")
plt.ylabel("xsection(pb).e-07")

#Plot-2
xsc_Vschi = df.plot(x = "Chi", y = "X-section(pb)", marker = "o", color = "red", linestyle="")
plt.grid(True)
plt.xlabel("chi.e+02")
plt.ylabel("xsection(pb).e-07")
