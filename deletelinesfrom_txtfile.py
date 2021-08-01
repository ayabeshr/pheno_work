import os

scan_dir = "/home/aya/HEP/toolbox-2.0.3/SSP-1.2.5/Output/THDMScan8M12lam2"

for i in range(501, 1001):
    sf = scan_dir + "/SPheno.spc.THDM__" + str(i)
    if (os.path.exists(sf)):
        sf2 = open(sf, "r")
        lines = sf2.readlines()
        sf2.close()
       #del lines[847:913]
       #del lines[842:847]
        """del lines[841]
        new_sf = open(sf, "w+")
        for line in lines:
            new_sf.write(line)
        new_sf.close() """
        sf3 = open(sf, "r")
        lines3 = sf3.readlines()
        if lines3[840] == "DECAY        16     0.00000000E+00   # Fv_3":
            print("spc file " + str(i) + " ends with DECAY 16")
        else:
            last_line = lines3[-1]
            print("ERROR!  spc file " + str(i) + " has no DECAY 16 end but instead : " + last_line)
        #print("spc file " + str(i) + " has been modified successfully")

