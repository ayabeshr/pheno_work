############################################
#   script for running MG5 at once for     #
#   several Benchmark points with varying  # 
#   parameter cards                        #                                       
############################################

import os 
import subprocess 

#subprocess.run("ls")
#subprocess.run(["./bin/mg5_aMC", "test2.txt"])

# define mandatory directories 
mg_dir = "/home/aya/programs/MG5_aMC_v2_7_3"
scan_dir = "/home/aya/HEP/toolbox-2.0.3/SSP-1.2.5/Output/THDMScan8M12lam2"

# creat a text file for saving process generation commands 
tf = open("/home/aya/programs/MG5_aMC_v2_7_3/gghhbb4Mu.txt", "w")
tf2 = mg_dir + "/gghhbb4Mu.txt"

# check if the text file has been created successfully 
if(os.path.exists(tf2)):
	print("File " + tf2 + " exists!! \n")
	#subprocess.run("ls")
else:
	print(" \n")
	print("File " + tf2 + " is not found! \n")
	open(mg_dir + "gghhbb4Mu.txt", 'w')
	subprocess.run("ls")
	
tf.write("import model 2HDM_20210524 -modelname \n")
tf.write("define j = g u1 u2 d1 d2 d3 u1bar u2bar d1bar d2bar d3bar \n")

# Iterate over spc files in a given range
for i in range(27, 101):
#for i in range(101, 1000):
#for i in range(1000, 2000):
#for i in range(2000, 4000):
#for i in range(4000, 6000):
#for i in range(6000, 8000):
#for i in range(8000, 9583):				
	sf = scan_dir + "/SPheno.spc.THDM__" + str(i) 
	if(os.path.isfile(sf)):
		tf.write("generate g g > h2 > h1 h1, h1 > d3 d3bar, (h1 > z z, z > e2bar e2, z > e2bar e2) \n")
		tf.write("output gg_hh_bb4Mu_20210731/gg_hh_bb4Mu_" + str(i) + "\n")
		l = ["launch \n", "0 \n", "set nevents 10 \n", "set ebeam1 7000 \n", "set ebeam2 7000 \n"]
		tf.writelines(l)
		tf.write(sf + " \n")
		tf.write("0 \n")  
		
tf.close()
print("Writing Txt file ENDS!")
#subprocess.run(["vim", "gghhbb4Mu.txt"])

subprocess.run(["./bin/mg5_aMC", "gghhbb4Mu.txt"])
print("Running MG ENDS!") 
	
