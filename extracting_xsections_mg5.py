#########################################
#  script for extracting x-sections from
#  several generated processes by MG
#########################################


import os 

prcs_mother_dir = "/home/Aya/MG5_aMC_v2_7_3/gg_hh_bb4Mu_20210801"

for i in range(121, 301):
	prcs_daughter_dir = prcs_mother_dir + "/gg_hh_bb4Mu_" + str(i) 
	if(os.path.exists(prcs_daughter_dir)):
		xsc_file = prcs_daughter_dir + "/Events/run_01/run_01_tag_1_banner.txt"
		#with open(xsc_file, "r") as xf:
		xf = open(xsc_file, "r")
		lines = xf.readlines() 
		for line in lines:
			if line[0:20] == "#  Integrated weight":
				print("Xsc for Point nb : [" + str(i) + "] is : " + line) 
				
		xf.close()
		
print("...Extracting Xsections Done...")
print("   ")
	
