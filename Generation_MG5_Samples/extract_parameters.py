#####################################
#   script for extracting info for  #
#   a given parameter from output   #
#   spectrum file from SSP_1.2.5    # 
#####################################


import os


spc_dir = "/home/Aya/THDM/THDMScan8M12lam2_copy" 
for i in range(54, 201):
	spc_file = spc_dir + "/SPheno.spc.THDM__" + str(i)
	if(os.path.exists(spc_file)):
		sf = open(spc_file, "r")
		lines = sf.readlines() 
		for line in lines:
			#f = open("/home/Aya/THDM/m12_param.txt", "w")
			#if line[27:35] == "M12input":
				#print(line[9:16])
			#if line[16:25] == "HS Mass 1":
				#print("mh1 for Point [" + str(i) + "] = " + line[4:11])
			if line[16:25] == "HS Mass 2":
				print("mh2 for Point [" + str(i) + "] = " + line[4:11])
		
		sf.close()
			#else:
				#print("ERROR in M12 input for Point [" + str(i) + "]")
