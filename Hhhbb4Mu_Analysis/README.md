Hhhbb4Mu Analysis: 

1- Studying expectation for BSM second heavy neutral scalar Higgs production and decay into two SM neutral scalar Higgs with final state 2b-jets+4Muons at CMS Phase2 upgrade with Lumimosity 3000 1/fb.

2- Signal and Background data files generated using MG5+Pythia8+Delphes.
 
 
BTagging in CMS Phase2 Delphes Simulation Cards:

1- There three Working Points WP which are Loose, Medium and Tight for b-jets, 
   defined in Delphes card, each WP is stored in a certain BitNumber position,
   for e.g; BitNumber=0 (default) means stored in 1st position, BitNumber=1 
   stored in 2nd position, and so on.....
   
2- Each BitNumber stores certain bit values, specifically defined in CMS Phase2
   Card as :
   
   
       WP       BitNumber    Bit Values
   
     Loose          0          1,3,5,7
     Medium         1          2,3,6,7
     Tight          2          4,5,6,7

   
    Bit Value       0    1    2    3    4    5    6    7
   
    Binary form    000  001  010  011  100  101  110  111 
   
   
3- Each jet has a certain bit value stored in a variable called Jet_BTag[] in 
   Jet container in Delphes tree.  
   
