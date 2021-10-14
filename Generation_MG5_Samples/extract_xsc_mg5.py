########################################
#  script for extracting x-sections from
#  several generated processes by MG
#########################################

import os, sys 


prcs_mother_dir = "/HEP_DATA/aya/signal_20210921_xsc_test"
xsc = [] 

# loop overall prcs daughter directories
for i in range(400, 410):
	
	prcs_daughter_dir = prcs_mother_dir + "/gg_hh_bb4Mu_" + str(i) 
	
	if(os.path.exists(prcs_daughter_dir)):
		
		xsc_file = prcs_daughter_dir + "/Events/run_01/run_01_tag_1_banner.txt"
		
		with open(xsc_file) as xf:
			lines = xf.readlines()
			for line in lines:
				if line[0:20] == "#  Integrated weight":
					xsection = line[35:55]
					xsc.append(xsection)

# save print output to a txt file
sys.stdout = open("all_xsc.txt", "w")

leng = len(xsc)		
i = 0 

while i < leng:
	print (xsc[i])
	i += 1
	
sys.stdout.close()
	
#print("...Extracting Xsections Done...")
#print("   ")
	
