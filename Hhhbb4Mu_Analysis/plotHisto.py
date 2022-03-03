import os
import ROOT as ru
#from cppyy.interactive import *
#from ROOT import *
from ROOT import TFile, TH1F, TCanvas, TLegend, THStack, TAxis, TH1, TNamed, TTree, TBranch, TLeaf 
from ROOT import kBlue, kGreen, kRed, kCyan, kYellow, kMagenta, kAzure, kOrange, kViolet
import rootlogon
rootlogon.setStyle()
from optparse import OptionParser
import distutils.util as ut 


# CL options can be used with running script

parser = OptionParser()

parser.add_option("-n", "--nbins", dest = "nbins",
                                   help = "set number of bins for plot",
                                   action = "store", type = "int")
                                   
parser.add_option("--xmin", dest = "xmin",
                           help = "set minimum value for x-axis variable",
                           action = "store", type = "float")  
                           
parser.add_option("--xmax", dest = "xmax",
                            help = "set maximum value for x-axis variable",
                            action = "store", type = "float")     
                            
parser.add_option("-o", "--output", dest = "output",
                                    help = "define output root file for saving results",
                                    action = "store", type = "string")    
                                    
parser.add_option("-t", "--tree", dest = "tree",
                                  help = "set tree name in sample input file",
                                  action = "store", type = "string", default = "output_demo")    # default is used in case the user doesn't supply an option   
                                  
parser.add_option("-v", "--variable", dest = "var",
                                      help = "set the variable to be plotted from current tree",
                                      action = "store", type = "string")  
                                  
parser.add_option("--logY", dest = "islogY",
                            help = "set y-axis of canvas or pad to log scale",
                            action = "store_true", default = False)  
                            
parser.add_option("-x", "--xlabel", dest = "xlb",
                                    help = "set x label of stacked histo",
                                    action = "store", type = "string")                                                                                                                                                                        
                           
parser.add_option("-y", "--ylabel", dest = "ylb",
                                    help = "set y label of stacked histo",
                                    action = "store", type = "string", default = "Events")                            
                           
(options, args) = parser.parse_args()

variable = options.var
nbin = options.nbins
xmini = options.xmin
xmaxi = options.xmax  
outfile = options.output  
intree = options.tree 
isLogY = options.islogY 
#xlb = options.                    


samples = "samples.txt"

Lumi = 3.e+03

output = ru.TFile(outfile , "RECREATE")

if (os.path.exists(samples)):     # check if samples configuration file really exists  
	
	smp = open(samples, "r")
	lines = smp.readlines()[1:]    # skip 1st line of file
	
	
	for line in lines:             # loop over info in samples confg file
		
		smp_confg = line.split("|")
		sample_path = smp_confg[0]          # path/to/input sample file
                xsec = float(smp_confg[1])
                Ngen = float(smp_confg[2])
                Nexp = Lumi * xsec
                scalefactor = Nexp / Ngen           # for reweighting samples
                histname = smp_confg[5]
                lineColor = smp_confg[3]
                lineWidth = int(smp_confg[4])
                
                #histo = ru.TH1F(histname, histname, nbin, xmini, xmaxi)
                histo = ru.TH1F(histname, histname, nbin, -10, xmaxi)
                eval("histo.SetLineColor(%s)"%lineColor)
                histo.SetLineWidth(lineWidth)
                histo.SetDirectory(0)              # to be sure that created histo not automatically disappear


	        print( "...working with sample: " + histname )
                print("Scale_Factor = " + str(scalefactor))

                input_file = ru.TFile(sample_path, "READ")
		tr = input_file.Get(intree)
                leaf = ru.TLeaf()
                leaf = tr.GetLeaf("LooseMuons", "lmu_pt")
                #leaf.SetAddress()

                for entry in tr:
                        var = eval("entry.%s"%variable)
                        #var = entry.leaf
                        for v in var:
                                histo.Fill(v, scalefactor)

                input_file.Close()
                output.cd()
                histo.Write()


	print("histos are written to: " + outfile)		
			
	output.Close()
      
      
	# load histos file
	histo_file = ru.TFile(outfile, "READ")   
	
	# call each histo

        h0 = ru.TH1F()
        h1 = ru.TH1F()
        h2 = ru.TH1F()
        h3 = ru.TH1F()
        h4 = ru.TH1F()
        h5 = ru.TH1F()
        h6 = ru.TH1F()
        h7 = ru.TH1F()
        h8 = ru.TH1F()
      
	h0 = histo_file.Get("ZZ4Mu")
	h1 = histo_file.Get("ZZZbb4Mu")
	h2 = histo_file.Get("SM_HTo4Mu")
	h3 = histo_file.Get("ZZbb2Mu")
        h4 = histo_file.Get("ZWpTobbMuNu")
        h5 = histo_file.Get("TTbar")
        h6 = histo_file.Get("DrellYan")
        h7 = histo_file.Get("SM_HTobb")
        h8 = histo_file.Get("HhhTobb4Mu")
        
	
	# define output file for final plot for all samples
	final_histo = ru.TFile("final_histo.root", "RECREATE")
	
	# define canvas 
	cnv = ru.TCanvas("canv", "canv", 400, 400)

        final_histo.cd()
        cnv.cd()

        h0.Draw("HIST")
        h1.Draw("HIST SAME")
        h2.Draw("HIST SAME")
        h3.Draw("HIST SAME")
        h4.Draw("HIST SAME")
        h5.Draw("HIST SAME")
        h6.Draw("HIST SAME")
        h7.Draw("HIST SAME")
        h8.Draw("HIST SAME")

	if ( isLogY ):              # set y-axis log scale for already plotted hs                  
		cnv.SetLogy()
		

	# define legend 
	leg = ru.TLegend(0.55,0.7,0.85,0.85)     #suitable dimensions for legend at upper left of pad

        leg.AddEntry(h0, "ZZ4Mu"       , "l")
	leg.AddEntry(h1, "ZZZbb4Mu"    , "l")
	leg.AddEntry(h2, "SM_HTo4Mu"   , "l")
	leg.AddEntry(h3, "ZZbb2Mu"     , "l")
        leg.AddEntry(h4, "ZWpTobbMuNu" , "l")
        leg.AddEntry(h5, "TTbar"       , "l")
        leg.AddEntry(h6, "DrellYan"    , "l")
        leg.AddEntry(h7, "SM_HTobb"    , "l")
        leg.AddEntry(h8, "HhhTobb4Mu"  , "l")

        
	#leg.SetHeader("(3000 fb^{-1}, 14TeV)","C")
        leg.SetMargin(0.4)
        leg.SetNColumns(2)
	#leg.SetColumnSeparation(0.1)
        leg.SetFillColor(0)
        leg.SetLineColor(0)
        
	leg.Draw("same")          # plot legend in the same pad of hs           
	
	final_histo.Close()
	histo_file.Close()
	
	if (os.path.exists("final_histo.root")):
		print("finshed plotting, plot saved in: final_histo.root")

