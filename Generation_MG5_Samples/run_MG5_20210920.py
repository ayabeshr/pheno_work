###################################################
#   script for running MG5 at once for several    #
#   Benchmark points with varying parameter cards # 
###################################################             

#---------------------------------------------------
# To run script in case of MG with versions < 2.8 
# (uses python2) open a terminal and do:   
#
#     cd /path/to/MG_directory                                                               
#     python /path/to/run_MG5.py                                                             
#----------------------------------------------------                                                                                           
# To run script in case of MG with versions >= 2.8.0 
# (uses python3) open a terminal and do: 
#                                                                                             
#     cd /path/to/MG_directory                                                               
#     python3 /path/to/run_MG5.py                                                            
#-----------------------------------------------------

import os 
#import subprocess 

#subprocess.run("ls")
#subprocess.run(["./bin/mg5_aMC", "test2.txt"])

# define mandatory directories 
mg_dir = "/home/Aya/MG5_aMC_v2_7_3"
scan_dir = "/home/Aya/Pheno_Work/THDMScan5M12lam2"

# creat a text file for saving process generation commands 
tf = open("/home/Aya/MG5_aMC_v2_7_3/gghhbb4Mu_new.txt", "w")
tf2 = mg_dir + "/gghhbb4Mu_new.txt"

# check if the text file has been created successfully 
if(os.path.exists(tf2)):
	print("File " + tf2 + " exists!! \n")
	#subprocess.run("ls")
else:
	print(" \n")
	print("File " + tf2 + " is not found! \n")
	open(mg_dir + "gghhbb4Mu_new.txt", 'w')
	subprocess.run("ls")
	
tf.write("import model 2HDM_20210524 --modelname \n")
#tf.write("define j = g u1 u2 d1 d2 d3 u1bar u2bar d1bar d2bar d3bar \n")

# Iterate over spc files in a given range
for i in range(0, 2001):
#for i in range(101, 501):
#for i in range(101, 1000):
#for i in range(1000, 2000):
#for i in range(2000, 4000):
#for i in range(4000, 6000):
#for i in range(6000, 8000):
#for i in range(8000, 9583):				
	sf = scan_dir + "/SPheno.spc.THDM__" + str(i) 
	if(os.path.isfile(sf)):
		tf.write("generate g g > h2 > h1 h1, h1 > d3 d3bar, (h1 > z > e2bar e2 e2bar e2) \n")  # syntax to apply on-shell and off-shell z
		tf.write("output /run/media/Aya/LinuxSpace/signal_20210921_xsc_test/gg_hh_bb4Mu_" + str(i) + "\n")
		l = ["launch \n", "0 \n", "set nevents 10 \n", "set ebeam1 7000 \n", "set ebeam2 7000 \n", "set maxjetflavor 5 \n"]
		tf.writelines(l)
		tf.write(sf + " \n")   # defines path to a valid param_card 
		tf.write("0 \n")  
		
tf.close()
print("Writing Txt file ENDS!")
#subprocess.run(["vim", "gghhbb4Mu.txt"])

# subprocess used in python3
#subprocess.run("alias python=python2")
#subprocess.run(["python", "./bin/mg5_aMC", "gghhbb4Mu.txt"])

# os.system() is used in python2
os.system("python2 ./bin/mg5_aMC gghhbb4Mu_new.txt") 
print("Running MG ENDS!")
