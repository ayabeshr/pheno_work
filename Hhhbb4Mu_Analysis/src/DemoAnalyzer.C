#define DemoAnalyzer_cxx

R__LOAD_LIBRARY(libPhysics)
//R__LOAD_LIBRARY(/home/aya/programs/Delphes-3.4.2/libDelphes.so)
R__ADD_LIBRARY_PATH(/home/aya/programs/Delphes-3.4.2)
R__LOAD_LIBRARY(libDelphes)

#include "/home/aya/Videos/DemoAnalyzer.h"
#include <TStyle.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>  
#include <TMath.h>
#include "THStack.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>


void DemoAnalyzer::Loop()
{
//   In a ROOT session, you can do:
//      root> .L DemoAnalyzer.C
//      root> DemoAnalyzer t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//    This is the loop skeleton where:
//     jentry is the global entry number in the chain
//     ientry is the entry number in the current Tree
//     Note that the argument to GetEntry must be:
//     jentry for TChain::GetEntry
//     ientry for TTree::GetEntry and TBranch::GetEntry
//
//    To read only selected branches, Insert statements like:
//     METHOD1:
//      fChain->SetBranchStatus("*",0);  // disable all branches
//      fChain->SetBranchStatus("branchname",1);  // activate branchname
//     METHOD2: replace line
//      fChain->GetEntry(jentry);       //read all branches
//      by  b_branchname->GetEntry(ientry); //read only this branch
  
   
	
  if (fChain == 0) return;
   
   
   // TH1D
   // various muon combinations bet 1,2,3,4 muons gives various Z masses, choose combination wich is nearer to Z mass 91 GeV to be Za, and the other with mass < 91 GeV take it to be Zb of signal and call them Za, Zb
   
   // size of available objects
   TH1D *h_muons_size_Loose;
   TH1D *h_jets_size;
   TH1D *h_MET_size;
   
   // MET    
   TH1D *h_MET;
   TH1D *h_eta_MET;
   TH1D *h_phi_MET;
   
   // All Muons   
   TH1D *h_pt_allMuons_Loose;
   TH1D *h_eta_allMuons_Loose;
   TH1D *h_phi_allMuons_Loose;
   
   // All Jets
   TH1D *h_pt_allJets;
   TH1D *h_eta_allJets;
   TH1D *h_phi_allJets;
   
   TH1D *h_mZ12;
   TH1D *h_mZ34;
   TH1D *h_mZ13;
   TH1D *h_mZ24;
   TH1D *h_mZ14;
   TH1D *h_mZ23;
      
   TH1D *h_DR_mu1mu2;     // define those histo
   TH1D *h_DR_mu3mu4;
   TH1D *h_DR_mu1mu3;
   TH1D *h_DR_mu2mu4;
   TH1D *h_DR_mu1mu4;
   TH1D *h_DR_mu2mu3;
      
   TH1D *h_mZa_4mu;             // Za and Zb are Z from h1->Za Zb
   TH1D *h_mZb_4mu;
   TH1D *h_pt_Za;
   TH1D *h_pt_Zb;
   TH1D *h_eta_Za;
   TH1D *h_eta_Zb;
   TH1D *h_phi_Za;
   TH1D *h_phi_Zb;
      
   // h1 -> Za Zb
   TH1D *h_mh1_ZaZb;
   TH1D *h_pt_h1_ZaZb;
   TH1D *h_eta_h1_ZaZb;
   TH1D *h_phi_h1_ZaZb;
      
   TH1D *h_mb_jet_1;
   TH1D *h_mb_jet_2;
   TH1D *h_pt_b_jet_1;
   TH1D *h_pt_b_jet_2;
   TH1D *h_eta_b_jet_1;
   TH1D *h_eta_b_jet_2;
   TH1D *h_phi_b_jet_1;
   TH1D *h_phi_b_jet_2;
   TH1D *h_DR_b1b2 ;
   
   // h1 -> b b~
   TH1D *h_mh1_b1b2;
   TH1D *h_pt_h1_b1b2;
   TH1D *h_eta_h1_b1b2;
   TH1D *h_phi_h1_b1b2;
   
   // h2 -> h1 h1 
   TH1D *h_mh2_h1h1;
   TH1D *h_pt_h2_h1h1;
   TH1D *h_eta_h2_h1h1;
   TH1D *h_phi_h2_h1h1;
      
   // Size of all Muons
   h_muons_size_Loose = new TH1D("h_muonsloose_size", "h_muonsloose_size", 10, 0., 10.);
   h_muons_size_Loose->GetXaxis()->SetTitle("Number of Muons");
   h_muons_size_Loose->GetYaxis()->SetTitle("Number of Events"); 
   
   // Size of all Jets
   h_jets_size = new TH1D("h_jets_size", "h_jets_size", 10, 0., 10.);
   h_jets_size->GetXaxis()->SetTitle("Number of Jets");
   h_jets_size->GetYaxis()->SetTitle("Number of Events"); 
   
   // Size of MET 
   h_MET_size = new TH1D("h_MET_size", "h_MET_size", 10, 0., 10.);
   h_MET_size->GetXaxis()->SetTitle("Number of MET");
   h_MET_size->GetYaxis()->SetTitle("Number of Events"); 
   
   // MET
   h_MET = new TH1D("h_MET", "h_MET", 100, 0., 200.);
   h_MET->GetXaxis()->SetTitle("MET[GeV]");
   h_MET->GetYaxis()->SetTitle("Number of Events");
   
   // Eta MET
   h_eta_MET = new TH1D("h_eta_MET", "h_eta_MET", 16, -8., 8.);
   h_eta_MET->GetXaxis()->SetTitle("#eta");
   h_eta_MET->GetYaxis()->SetTitle("Number of Events");
   
   // Phi MET 
   h_phi_MET = new TH1D("h_phi_MET", "h_phi_MET", 20, -10., 10.);
   h_phi_MET->GetXaxis()->SetTitle("#phi");
   h_phi_MET->GetYaxis()->SetTitle("Number of Events");
   
   // pT of all Muons
   h_pt_allMuons_Loose = new TH1D("h_pt_allMuonsLoose", "h_pt_allMuonsLoose", 100, 0., 100.);
   h_pt_allMuons_Loose->GetXaxis()->SetTitle("p_{T} [GeV/c]");
   h_pt_allMuons_Loose->GetYaxis()->SetTitle("Number of Events");
   
   // Eta of all Muons
   h_eta_allMuons_Loose = new TH1D("h_eta_allMuonsLoose", "h_eta_allMuonsLoose", 16, -8., 8.);
   h_eta_allMuons_Loose->GetXaxis()->SetTitle("#eta");
   h_eta_allMuons_Loose->GetYaxis()->SetTitle("Number of Events");
   
   // Phi of all Muons
   h_phi_allMuons_Loose= new TH1D("h_phi_allMuonsLoose", "h_phi_allMuonsLoose", 20, -10., 10.);
   h_phi_allMuons_Loose->GetXaxis()->SetTitle("#phi");
   h_phi_allMuons_Loose->GetYaxis()->SetTitle("Number of Events");

   // pT of all Jets
   h_pt_allJets = new TH1D("h_pt_allJets", "h_pt_allJets", 100, 0., 100.);
   h_pt_allJets->GetXaxis()->SetTitle("p_{T} [GeV/c]");
   h_pt_allJets->GetYaxis()->SetTitle("Number of Events");
   
   // Eta of all Jets
   h_eta_allJets = new TH1D("h_eta_allJets", "h_eta_allJets", 16, -8., 8.);
   h_eta_allJets->GetXaxis()->SetTitle("#eta");
   h_eta_allJets->GetYaxis()->SetTitle("Number of Events");
   
   // Phi of all Jets
   h_phi_allJets= new TH1D("h_phi_allJets", "h_phi_allJets", 20, -10., 10.);
   h_phi_allJets->GetXaxis()->SetTitle("#phi");
   h_phi_allJets->GetYaxis()->SetTitle("Number of Events");

   // Combinations 1234
   
   // Z12
   h_mZ12 = new TH1D("h_mZ12", "h_mZ12", 75, 0., 150.);
   h_mZ12->GetXaxis()->SetTitle("Invariant mass of dimuon [GeV/c^{2}]"); 
   h_mZ12->GetYaxis()->SetTitle("Number of Events");
   
   // Z34
   h_mZ34 = new TH1D("h_mZ34", "h_mZ34", 75, 0., 150.);
   h_mZ34->GetXaxis()->SetTitle("Invariant mass of dimuon [GeV/c^{2}]"); 
   h_mZ34->GetYaxis()->SetTitle("Number of Events");
   
   // Combinations 1324
   
   // Z13
   h_mZ13 = new TH1D("h_mZ13", "h_mZ13", 75, 0., 150.);
   h_mZ13->GetXaxis()->SetTitle("Invariant mass of dimuon [GeV/c^{2}]"); 
   h_mZ13->GetYaxis()->SetTitle("Number of Events");
   
   // Z24
   h_mZ24 = new TH1D("h_mZ24", "h_mZ24", 75, 0., 150.);
   h_mZ24->GetXaxis()->SetTitle("Invariant mass of dimuon [GeV/c^{2}]"); 
   h_mZ24->GetYaxis()->SetTitle("Number of Events");
   
   // Combinations 1423
   
   // Z14
   h_mZ14 = new TH1D("h_mZ14", "h_mZ14", 75, 0., 150.);
   h_mZ14->GetXaxis()->SetTitle("Invariant mass of dimuon [GeV/c^{2}]"); 
   h_mZ14->GetYaxis()->SetTitle("Number of Events");
   
   // Z23
   h_mZ23 = new TH1D("h_mZ23", "h_mZ23", 75, 0., 150.);
   h_mZ23->GetXaxis()->SetTitle("Invariant mass of dimuon [GeV/c^{2}]"); 
   h_mZ23->GetYaxis()->SetTitle("Number of Events");
   
   // Za 
   
   // Za Invariant Mass of dimuons closest to Z mass ~ 91. GeV
   h_mZa_4mu = new TH1D("h_mZa_4mu", "h_mZa_4mu_closest to Z Mass", 120, 0., 120.);
   h_mZa_4mu->GetXaxis()->SetTitle("Invariant mass of dimuon [GeV/c^{2}]"); 
   h_mZa_4mu->GetYaxis()->SetTitle("Number of Events");
   
   h_pt_Za = new TH1D("h_pt_Za", "h_pt_Za", 500, 0., 500.); 
   h_pt_Za->GetXaxis()->SetTitle("p_{T} [GeV/c]");
   h_pt_Za->GetYaxis()->SetTitle("Number of Events");
   
   h_eta_Za = new TH1D("h_eta_Za", "h_eta_Za", 16, -8., 8.); 
   h_eta_Za->GetXaxis()->SetTitle("#eta");
   h_eta_Za->GetYaxis()->SetTitle("Number of Events");
   
   h_phi_Za = new TH1D("h_phi_Za", "h_phi_Za", 20, -10., 10.); 
   h_phi_Za->GetXaxis()->SetTitle("#phi");
   h_phi_Za->GetYaxis()->SetTitle("Number of Events");
     
   // Zb
   
   // Za Invariant Mass of dimuons not close to Z mass ~ 91. GeV 
   h_mZb_4mu = new TH1D("h_mZb_4mu", "h_mZb_4mu_Not Close to Z Mass", 120, 0., 120.);
   h_mZb_4mu->GetXaxis()->SetTitle("Invariant mass of dimuon [GeV/c^{2}]"); 
   h_mZb_4mu->GetYaxis()->SetTitle("Number of Events");
   
   h_pt_Zb = new TH1D("h_pt_Zb", "h_pt_Zb", 500, 0., 500.); 
   h_pt_Zb->GetXaxis()->SetTitle("p_{T} [GeV/c]");
   h_pt_Zb->GetYaxis()->SetTitle("Number of Events");
   
   h_eta_Zb = new TH1D("h_eta_Zb", "h_eta_Zb", 16, -8., 8.); 
   h_eta_Zb->GetXaxis()->SetTitle("#eta");
   h_eta_Zb->GetYaxis()->SetTitle("Number of Events");
   
   h_phi_Zb = new TH1D("h_phi_Zb", "h_phi_Zb", 20, -10., 10.); 
   h_phi_Zb->GetXaxis()->SetTitle("#phi");
   h_phi_Zb->GetYaxis()->SetTitle("Number of Events");
     
   h_mb_jet_1 = new TH1D("h_mb_jet_1", "h_mb_jet_1", 100, 0., 100.);
   h_mb_jet_1->GetXaxis()->SetTitle(" mass of b1 jet [GeV/C^{2}]");
   h_mb_jet_1->GetYaxis()->SetTitle("Number of Events");
      
   h_mb_jet_2 = new TH1D("h_mb_jet_2", "h_mb_jet_2", 100, 0., 100.);
   h_mb_jet_2->GetXaxis()->SetTitle(" mass of b2 jet [GeV/C^{2}]");
   h_mb_jet_2->GetYaxis()->SetTitle("Number of Events");
      
   h_pt_b_jet_1 = new TH1D("h_pt_b_jet_1", "h_pt_b_jet_1", 500, 0., 500.);
   h_pt_b_jet_1->GetXaxis()->SetTitle("p_{T} (GeV/C)");
   h_pt_b_jet_1->GetYaxis()->SetTitle("Number of Events");
   
   h_pt_b_jet_2 = new TH1D("h_pt_b_jet_2", "h_pt_b_jet_2", 500, 0., 500.);
   h_pt_b_jet_2->GetXaxis()->SetTitle("p_{T} (GeV/C)");
   h_pt_b_jet_2->GetYaxis()->SetTitle("Number of Events");
   
   h_eta_b_jet_1 = new TH1D("h_eta_b_jet_1", "h_eta_b_jet_1", 20, -8., 8.);
   h_eta_b_jet_1->GetXaxis()->SetTitle("#eta");
   h_eta_b_jet_1->GetYaxis()->SetTitle("Number of Events");
   
   h_eta_b_jet_2 = new TH1D("h_eta_b_jet_2", "h_eta_b_jet_2", 20, -8., 8.);
   h_eta_b_jet_2->GetXaxis()->SetTitle("#eta");
   h_eta_b_jet_2->GetYaxis()->SetTitle("Number of Events");
   
   h_phi_b_jet_1 = new TH1D("h_phi_b_jet_1", "h_phi_b_jet_1", 20, -10., 10.);
   h_phi_b_jet_1->GetXaxis()->SetTitle("#phi");
   h_phi_b_jet_1->GetYaxis()->SetTitle("Number of Events");
   
   h_phi_b_jet_2 = new TH1D("h_phi_b_jet_2", "h_phi_b_jet_2", 20, -10., 10.);
   h_phi_b_jet_2->GetXaxis()->SetTitle("#phi");
   h_phi_b_jet_2->GetYaxis()->SetTitle("Number of Events");
     
   h_DR_b1b2 = new TH1D("h_DR_b1b2", "DR between b1 b2 jets", 20, 0., 20.);
   h_DR_b1b2->GetXaxis()->SetTitle("#Delta R");
   h_DR_b1b2->GetYaxis()->SetTitle("Number of Events");
   
   // DR seperation bet 4Muons  
   h_DR_mu1mu2 = new TH1D("h_DR_mu1mu2", "DR between muon1 muon2", 20, 0., 20.);
   h_DR_mu1mu2->GetXaxis()->SetTitle("#Delta R_{12}");
   h_DR_mu1mu2->GetYaxis()->SetTitle("Number of Events");
      
   h_DR_mu3mu4 = new TH1D("h_DR_mu3mu4", "DR between muon3 muon4", 20, 0., 10.);
   h_DR_mu3mu4->GetXaxis()->SetTitle("#Delta R_{34}");
   h_DR_mu3mu4->GetYaxis()->SetTitle("Number of Events");
      
   h_DR_mu1mu3 = new TH1D("h_DR_mu1mu3", "DR between muon1 muon3", 20, 0., 10.);
   h_DR_mu1mu3->GetXaxis()->SetTitle("#Delta R_{13}");
   h_DR_mu1mu3->GetYaxis()->SetTitle("Number of Events");
      
   h_DR_mu2mu4 = new TH1D("h_DR_mu2mu4", "DR between muon2 muon4", 20, 0., 10.);
   h_DR_mu2mu4->GetXaxis()->SetTitle("#Delta R_{24}");
   h_DR_mu2mu4->GetYaxis()->SetTitle("Number of Events");
   
   h_DR_mu1mu4 = new TH1D("h_DR_mu1mu4", "DR between muon1 muon4", 20, 0., 10.);
   h_DR_mu1mu4->GetXaxis()->SetTitle("#Delta R_{14}");
   h_DR_mu1mu4->GetYaxis()->SetTitle("Number of Events");
   
   h_DR_mu2mu3 = new TH1D("h_DR_mu2mu3", "DR between muon2 muon3", 20, 0., 10.);
   h_DR_mu2mu3->GetXaxis()->SetTitle("#Delta R_{23}");
   h_DR_mu2mu3->GetYaxis()->SetTitle("Number of Events");
   
   
   // DR seperation bet 2 b jets 
   h_DR_b1b2 = new TH1D("h_DR_b1b2", "h_DR_b1b2", 20, 0., 10.); 
   h_DR_b1b2->GetXaxis()->SetTitle("#Delta R");
   h_DR_b1b2->GetYaxis()->SetTitle("Number of Events");
   
   // h1 -> Za Zb
   h_mh1_ZaZb = new TH1D("h_mh1_ZaZb", "h_mh1_ZaZb", 100, 0., 200.);
   h_mh1_ZaZb->GetXaxis()->SetTitle("Invariant mass of 4 muons (GeV/c^{2})");
   h_mh1_ZaZb->GetYaxis()->SetTitle("Number of Events");
   
   h_pt_h1_ZaZb = new TH1D("h_pt_h1_ZaZb", "h_pt_h1_ZaZb", 500, 0., 500.);
   h_pt_h1_ZaZb->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   h_pt_h1_ZaZb->GetYaxis()->SetTitle("Number of Events");
   
   h_eta_h1_ZaZb = new TH1D("h_eta_h1_ZaZb", "h_eta_h1_ZaZb", 16, -8., 8.);
   h_eta_h1_ZaZb->GetXaxis()->SetTitle("#eta");
   h_eta_h1_ZaZb->GetYaxis()->SetTitle("Number of Events");
      
   h_phi_h1_ZaZb = new TH1D("h_phi_h1_ZaZb", "h_phi_h1_ZaZb", 20, -10., 10.);
   h_phi_h1_ZaZb->GetXaxis()->SetTitle("#phi");
   h_phi_h1_ZaZb->GetYaxis()->SetTitle("Number of Events");
   
   // h1 -> b b~
   h_mh1_b1b2 = new TH1D("h_mh1_b1b2", "h_mh1_b1b2", 100, 0., 200.);
   h_mh1_b1b2->GetXaxis()->SetTitle("Invariant mass of 2 b-jets (GeV/c^{2})");
   h_mh1_b1b2->GetYaxis()->SetTitle("Number of Events");
   
   h_pt_h1_b1b2 = new TH1D("h_pt_h1_b1b2", "h_pt_h1_b1b2", 500, 0., 500.);
   h_pt_h1_b1b2->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   h_pt_h1_b1b2->GetYaxis()->SetTitle("Number of Events");
   
   h_eta_h1_b1b2 = new TH1D("h_eta_h1_b1b2", "h_eta_h1_b1b2", 16, -8., 8.);
   h_eta_h1_b1b2->GetXaxis()->SetTitle("#eta");
   h_eta_h1_b1b2->GetYaxis()->SetTitle("Number of Events"); 
      
   h_phi_h1_b1b2 = new TH1D("h_phi_h1_b1b2", "h_phi_h1_b1b2", 20, -10., 10.);
   h_phi_h1_b1b2->GetXaxis()->SetTitle("#phi");
   h_phi_h1_b1b2->GetYaxis()->SetTitle("Number of Events"); 
   
   // h2 -> h1 h1 
   h_mh2_h1h1 = new TH1D("h_mh2_h1h1", "h_mh2_h1h1", 350, 0., 700.);
   h_mh2_h1h1->GetXaxis()->SetTitle("Invariant mass of 2 b-jets + 4 muons (GeV/C^{2})");
   h_mh2_h1h1->GetYaxis()->SetTitle("Number of Events"); 
   
   h_pt_h2_h1h1 = new TH1D("h_pt_h2_h1h1", "h_pt_h2_h1h1", 500, 0., 500.);
   h_pt_h2_h1h1->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   h_pt_h2_h1h1->GetYaxis()->SetTitle("Number of Events");
   
   h_eta_h2_h1h1 = new TH1D("h_eta_h2_h1h1", "h_eta_h2_h1h1", 16, -8., 8.);
   h_eta_h2_h1h1->GetXaxis()->SetTitle("#eta");
   h_eta_h2_h1h1->GetYaxis()->SetTitle("Number of Events");
   
   h_phi_h2_h1h1 = new TH1D("h_phi_h2_h1h1", "h_phi_h2_h1h1", 16, -8., 8.);
   h_phi_h2_h1h1->GetXaxis()->SetTitle("#phi");
   h_phi_h2_h1h1->GetYaxis()->SetTitle("Number of Events");
   
   
   //================================================================================================//
   //                                   DATASET SIGNAL H->hh->bb4Mu  NoPU                            //
   //================================================================================================//  
   TFile *indata = new TFile("/media/aya/PACKUP/Aya/signal_20210831/BMP1_hh_bb_4Mu.root", "READ");  
   
   //================================================================================================//
   //                                   DATASET Background 14 TeV                                    //
   //================================================================================================// 
   // DY 
   //TFile *input_file = new TFile("/media/aya/LinuxSpace/MyWork_Final_Samples/DY_BG_14TeV_SMfull_GEN/DY_BG_14TeV_SMfull_pythia8_CMSPhaseII-0PU_GEN-SIM.root", "READ");
   
   
   // Results root file
   //TFile *op_file = new TFile("/home/aya/Desktop/Pheno_Work/analysis/Results/output_demo_DY.root", "RECREATE"); 
   TFile *op_file = new TFile("/home/aya/Desktop/Pheno_Work/analysis/Results/output_demo_BMP1_hh_bb_4Mu.root", "RECREATE"); 
   
   //------------------------WEIGHT Calculation---------------------------
  
   float Lumi_data = 3.e+3;    // in 1/fb
   //Lumi_mc = nEvents/xsection(fb);
  
   // signal
   float Lumi_mc = 100000./42.32e-11;  // BMP1
   
   //float Lumi_mc = 1000000./898225.;   // DY
   
    
   float wt = Lumi_data/Lumi_mc;
   
   
   /*===================================================================================*/  
  /*------------------------------Looping over ALL Events-----------------------------*/
  /*==================================================================================*/ 
   
   
   Long64_t nentries = fChain->GetEntriesFast();
   //Long64_t nentries = 10;
   Long64_t nbytes = 0, nb = 0;
   
   for (Long64_t jentry = 0; jentry < nentries; jentry++) {
   
	 cout << "******START EVENT LOOP!******    ,    Event nb = " << jentry << endl; 
	 Long64_t ientry = LoadTree(jentry);
         if (ientry < 0) break;
         nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;
	   
	// Loop over all Muons Loose
         cout << "start loop overall loose muons" << endl;
         
	 for (Int_t i = 0; i < MuonLoose_size; i++){
	       
		 h_muons_size_Loose->Fill(i);
		 h_pt_allMuons_Loose->Fill(MuonLoose_PT[i], wt);
		 h_eta_allMuons_Loose->Fill(MuonLoose_Eta[i], wt);
		 h_phi_allMuons_Loose->Fill(MuonLoose_Phi[i], wt);
	  }
	  cout << "end loop overall loose muons" << endl;
      
         // Loop overall Jets
          cout << "start loop overall jets" << endl;
          for (Int_t i = 0; i < Jet_size; i++){
	        
		h_jets_size->Fill(i);
	        h_pt_allJets->Fill(Jet_PT[i], wt);
	        h_eta_allJets->Fill(Jet_Eta[i], wt);
	        h_phi_allJets->Fill(Jet_Phi[i], wt);
          }
          cout << "end loop overall jets" << endl;
      
	   
          // Loop over MET
          cout << "start loop overall MET" << endl;
          for (Int_t i = 0; i < MissingET_size; i++){
	      
		h_MET_size->Fill(i);
	        h_MET->Fill(MissingET_MET[i], wt);
	        h_eta_MET->Fill(MissingET_Eta[i], wt);
                h_eta_MET->Fill(MissingET_Phi[i], wt);
          }
          cout << "end loop overall MET" << endl;
      
	   
      
		  
		  if ( MuonLoose_size > 3 ){ // I have at least 4 Muons per event 
			  
			  // TLorentzVector declarations 
              TLorentzVector mu1, mu2, mu3, mu4, Z12, Z34, Z13, Z24, Z14, Z23, Za, Zb, b1, b2, h1, h2, H;
      
              // Set Pt, eta, phi mass for 4 muons TLV
              double muon_mass = 0.105658375;  // in GeV
              //double m_mu1, m_mu2, m_mu3, m_mu4;
              double pt_mu1, pt_mu2, pt_mu3, pt_mu4;
              double eta_mu1, eta_mu2, eta_mu3, eta_mu4;
              double phi_mu1, phi_mu2, phi_mu3, phi_mu4;
      
              // Initialize variables
              pt_mu1 = -9999.; pt_mu2 = -9999.; pt_mu3 = -9999.; pt_mu4 = -9999.;
              eta_mu1 = -9999.; eta_mu2 = -9999.; eta_mu3 = -9999.; eta_mu4 = -9999.;
              phi_mu1 = -9999.; phi_mu2 = -9999.; phi_mu3 = -9999.; phi_mu4 = -9999.;
      
              mu1.SetPtEtaPhiM(MuonLoose_PT[0], MuonLoose_Eta[0], MuonLoose_Phi[0], muon_mass); 
              mu2.SetPtEtaPhiM(MuonLoose_PT[1], MuonLoose_Eta[1], MuonLoose_Phi[1], muon_mass);
              mu3.SetPtEtaPhiM(MuonLoose_PT[2], MuonLoose_Eta[2], MuonLoose_Phi[2], muon_mass);
              mu4.SetPtEtaPhiM(MuonLoose_PT[3], MuonLoose_Eta[3], MuonLoose_Phi[3], muon_mass);  
      
              pt_mu1 = mu1.Pt();
              eta_mu1 = mu1.Eta();
              phi_mu1 = mu1.Phi();
      
              pt_mu2 = mu2.Pt();
              eta_mu2 = mu2.Eta();
              phi_mu2 = mu2.Phi();
      
              pt_mu3 = mu3.Pt();
              eta_mu3 = mu3.Eta();
              phi_mu3 = mu3.Phi();
      
              pt_mu4 = mu2.Pt();
              eta_mu4 = mu2.Eta();
              phi_mu4 = mu2.Phi();
      
              // Calculate DR bet each 2 muons for all possible combinations 1234, 1324, 1423
              double DR_mu12, DR_mu34, DR_mu13, DR_mu24, DR_mu14, DR_mu23; 
      
              // Initialize variables 
              DR_mu12 = -9999.; DR_mu34 = -9999.; DR_mu13 = -9999.; DR_mu24 = -9999.; DR_mu14 = -9999.; DR_mu23 = -9999.; 
      
              DR_mu12 = TMath::Sqrt(TMath::Power((eta_mu1 - eta_mu2), 2) + TMath::Power((phi_mu1 - phi_mu2), 2));
              DR_mu34 = TMath::Sqrt(TMath::Power((eta_mu3 - eta_mu4), 2) + TMath::Power((phi_mu3 - phi_mu4), 2));
              DR_mu13 = TMath::Sqrt(TMath::Power((eta_mu1 - eta_mu3), 2) + TMath::Power((phi_mu1 - phi_mu3), 2));
              DR_mu24 = TMath::Sqrt(TMath::Power((eta_mu2 - eta_mu4), 2) + TMath::Power((phi_mu2 - phi_mu4), 2));
              DR_mu14 = TMath::Sqrt(TMath::Power((eta_mu1 - eta_mu4), 2) + TMath::Power((phi_mu1 - phi_mu4), 2));
              DR_mu23 = TMath::Sqrt(TMath::Power((eta_mu2 - eta_mu3), 2) + TMath::Power((phi_mu2 - phi_mu3), 2));
      
              h_DR_mu1mu2->Fill(DR_mu12, wt);
              h_DR_mu3mu4->Fill(DR_mu34, wt);
              h_DR_mu1mu3->Fill(DR_mu13, wt);
              h_DR_mu2mu4->Fill(DR_mu24, wt);
              h_DR_mu1mu4->Fill(DR_mu14, wt);
              h_DR_mu2mu3->Fill(DR_mu23, wt); 
       
       
             /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
              ^                                         ^ 
              ^            Determine Za, Zb             ^
              ^                                         ^
              ^              for h -> Z Z               ^
              ^                                         ^ 
              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
      
      
              // 4 Muons various combination pairs 1234, 1324, 1423
              double mZ = 91.1876; // in GeV 
              double mZ12, mZ34, mZ13, mZ24, mZ14, mZ23, dZc1, dZc2, dZc3;
              double pt_Z12, pt_Z34, pt_Z13, pt_Z24, pt_Z14, pt_Z23;
              double eta_Z12, eta_Z34, eta_Z13, eta_Z24, eta_Z14, eta_Z23;
              double phi_Z12, phi_Z34, phi_Z13, phi_Z24, phi_Z14, phi_Z23;
              double dZ12, dZ23, dZ34, dZ13, dZ14, dZ24;
              double mZa, mZb, ptZa, ptZb, etaZa, etaZb, phiZa, phiZb;
      
              // Initialize variables
              mZ12 = -9999.; mZ34 = -9999.; mZ13 = -9999.; mZ24 = -9999.; mZ14 = -9999.; mZ23 = -9999.; dZc1 = -9999.; dZc2 = -9999.; dZc3 = -9999.;
              pt_Z12 = -9999.; pt_Z34 = -9999.; pt_Z13 = -9999.; pt_Z24 = -9999.; pt_Z14 = -9999.; pt_Z23 = -9999.;
              eta_Z12 = -9999.; eta_Z34 = -9999.; eta_Z13 = -9999.; eta_Z24 = -9999.; eta_Z14 = -9999.; eta_Z23 = -9999.;
              phi_Z12 = -9999.; phi_Z34 = -9999.; phi_Z13 = -9999.; phi_Z24 = -9999.; phi_Z14 = -9999.; phi_Z23 = -9999.;
              dZ12 = -9999.; dZ23 = -9999.; dZ34 = -9999.; dZ13 = -9999.; dZ14 = -9999.; dZ24 = -9999.;
              mZa = -9999.; mZb = -9999.; ptZa = -9999.; ptZb = -9999.; etaZa = -9999.; etaZb = -9999.; phiZa = -9999.; phiZb = -9999.;
      
      
              // Sure that muon pairs are of opposite signs 
              if ( MuonLoose_Charge[0] + MuonLoose_Charge[1] + MuonLoose_Charge[2] + MuonLoose_Charge[3] == 0){

		   // Start 4 muon combination 1234 
		   if ( MuonLoose_Charge[0] + MuonLoose_Charge[1] == 0){  // mu1, mu2 
			  
                        if ( MuonLoose_Charge[2] + MuonLoose_Charge[3] == 0){  // mu3, mu4 
				  
		             Z12 = mu1 + mu2;
                             mZ12 = Z12.M();
                             pt_Z12 = Z12.Pt();
                             eta_Z12 = Z12.Eta();
                             phi_Z12 = Z12.Phi();
                  
                             Z34 = mu3 + mu4;
                             mZ34 = Z34.M();
                             pt_Z34 = Z34.Pt();
                             eta_Z34 = Z34.Eta();
                             phi_Z34 = Z34.Phi();
			      
			     if (mZ12 > 0.) h_mZ12->Fill(mZ12, wt);
			     if (mZ34 > 0.) h_mZ34->Fill(mZ34, wt);
			     
			 }
	            }
		  
		    dZ12 = fabs(mZ12 - mZ);
		    dZ34 = fabs(mZ34 - mZ);
		    
                    // condition ? result_if_true : result_if_false  -> syntax for using ? conditional operator 
		    //dZc1 = ( dZ12 < dZ34 ) ? dZ12 : dZ34;
		    if ( dZ12 < dZ34 ) {
				      
			 dZc1 = dZ12;
		         cout << "dZ mass for combination 1234 = dZ12 = " << dZc1 << endl;
	            }	
		    
		    else 
		    {
			 dZc1 = dZ34;
		         cout << "dZ mass for combination 1234 = dZ34 = " << dZc1 << endl;
		    }
		    
		    // Start 4 muon combination 1324 
		    if ( MuonLoose_Charge[0] + MuonLoose_Charge[2] == 0){  // mu1, mu3
				
		         if ( MuonLoose_Charge[1] + MuonLoose_Charge[3] == 0){ // mu2, mu4
					
			       Z13 = mu1 + mu3;
			       mZ13 = Z13.M();
			       pt_Z13 = Z13.Pt();
                               eta_Z13 = Z13.Eta();
                               phi_Z13 = Z13.Phi();
					
			       Z24 = mu2 + mu4;
			       mZ24 = Z24.M();
			       pt_Z24 = Z24.Pt();
                               eta_Z24 = Z24.Eta();
                               phi_Z24 = Z24.Phi();
					
			       if (mZ13 > 0.) h_mZ13->Fill(mZ13, wt);
			       if (mZ24 > 0.) h_mZ24->Fill(mZ24, wt);
			   }
		      }
		    
		      dZ13 = fabs(mZ13 - mZ);
		      dZ24 = fabs(mZ24 - mZ);
		
		      //dZc2 = ( dZ13 < dZ24 ) ? dZ13 : dZ24; 
		      if ( dZ13 < dZ24 ) {
				      
			   dZc2 = dZ13;
		           cout << "dZ mass for combination 1324 = dZ13 = " << dZc2 << endl;
		      }	
		       
		      else 
		      {
			    dZc2 = dZ24;
		            cout << "dZ mass for combination 1324 = dZ24 = " << dZc2 << endl;
		      }
		    
		      // Start 4 muon combination 1423 
		      if ( MuonLoose_Charge[0] + MuonLoose_Charge[3] == 0){  // mu1, mu4
				
			   if ( MuonLoose_Charge[1] + MuonLoose_Charge[2] == 0){ // mu2, mu3
					
			        Z14 = mu1 + mu4;
			        mZ14 = Z14.M();
		                pt_Z14 = Z14.Pt();
			        eta_Z14 = Z14.Eta();
			        phi_Z14 = Z14.Phi();
					
			        Z23 = mu2 + mu3;
				mZ23 = Z23.M();
			        pt_Z23 = Z23.Pt();
			        eta_Z23 = Z23.Eta();
		                phi_Z23 = Z23.Phi();
					
			        if (mZ14 > 0.) h_mZ14->Fill(mZ14, wt);
			        if (mZ23 > 0.) h_mZ23->Fill(mZ23, wt);
			    }
		       }
		    
		       dZ14 = fabs(mZ14 - mZ);
		       dZ23 = fabs(mZ23 - mZ);
		
		       //dZc3 = ( dZ14 < dZ23 ) ? dZ14 : dZ23; 
		       if ( dZ14 < dZ23 ) {
				      
			    dZc3 = dZ14;
		            cout << "dZ mass for combination 1423 = dZ14 = " << dZc3 << endl;
		       }	
		         
		       else 
		       {
			    dZc3 = dZ23;
		            cout << "dZ mass for combination 1423 = dZ23 = " << dZc3 << endl;
		       } 
       
       
                       // Choose dZc of the least value bet dZc1, dZc2, dZc3 to be Za, Zb combination
		       if ( dZc1 < dZc2 && dZc1 < dZc3 ){  // dZc1 < dZc2, dZc3
				
		            if ( dZ12 < dZ34 ){
					
			         Za = Z12;
			         mZa = mZ12;
                                 ptZa = pt_Z12;
                                 etaZa = eta_Z12;
                                 phiZa = phi_Z12;
					
			         Zb = Z34;
				 mZb = mZ34;
                                 ptZb = pt_Z34;
                                 etaZb = eta_Z34;
			         phiZb = phi_Z34;
					
			     }
				
			     else
		             {
			          Za = Z34;
			          mZa = mZ34;
                                  ptZa = pt_Z34;
                                  etaZa = eta_Z34;
                                  phiZa = phi_Z34;
                    
				  Zb = Z12;
			          mZb = mZ12;
                                  ptZb = pt_Z12;
                                  etaZb = eta_Z12;
                                  phiZb = phi_Z12;
					
			      } 
			
			 } // end if dZc1   
			
			 else if ( dZc2 < dZc1 && dZc2 < dZc3 ){  // dZc2 < dZc1, dZc3
				
				   if ( dZ13 < dZ24 ){
					
					Za = Z13;
					mZa = mZ13;
                                        ptZa = pt_Z13;
                                        etaZa = eta_Z13;
                                        phiZa = phi_Z13;
					
					Zb = Z24;
					mZb = mZ24;
                                        ptZb = pt_Z24;
                                        etaZb = eta_Z24;
                                        phiZb = phi_Z24;
                    
				    }
				 
				    else
				    {      
					 Za = Z24;
					 mZa = mZ24;
                                         ptZa = pt_Z24;
                                         etaZa = eta_Z24;
                                         phiZa = phi_Z24;
					
					 Zb = Z13;
					 mZb = mZ13;
                                         ptZb = pt_Z13;
                                         etaZb = eta_Z13;
                                         phiZb = phi_Z13;
					
				    }
				
			   } // end else if dZc2
		    
		           else 
			   {  // dZc3 < dZc1, dZc2
					
			        if ( dZ14 < dZ23 ){
					
			             Za = Z14;
		                     mZa = mZ14;
                                     ptZa = pt_Z14;
                                     etaZa = eta_Z14;
                                     phiZa = phi_Z14;
					
				     Zb = Z23;
			             mZb = mZ23;
                                     ptZb = pt_Z23;
                                     etaZb = eta_Z23;
                                     phiZb = phi_Z23;
					
				 }
				
				 else
				 {
			              Za = Z23;
			              mZa = mZ23;
                                      ptZa = pt_Z23;
                                      etaZa = eta_Z23;
                                      phiZa = phi_Z23;
					
				      Zb = Z14;
			              mZb = mZ14;
                                      ptZb = pt_Z14;
                                      etaZb = eta_Z14;
                                      phiZb = phi_Z14;
				   	
		                  }
				
			    } // end else dZc3
			
                            if ( mZa > 40. && mZa < 120.){
				if ( mZb > 12. && mZb < 120. ){
					
				       h_mZa_4mu->Fill(mZa, wt);            
                                       h_mZb_4mu->Fill(mZb, wt);
                                       h_pt_Za->Fill(ptZa, wt);
                                       h_pt_Zb->Fill(ptZb, wt);
                                       h_eta_Za->Fill(etaZa, wt);
                                       h_eta_Zb->Fill(etaZb, wt);
		                       h_phi_Za->Fill(phiZa, wt);
		                       h_phi_Zb->Fill(phiZb, wt); 
		            
		            
		           //============================
		           // Reconstruct h1 from Za, Zb
		           //============================
            
                           double mh1_ZaZb, pt_h1_ZaZb, eta_h1_ZaZb, phi_h1_ZaZb;
                   
                           // Initialize Variables
                           mh1_ZaZb = -9999.; pt_h1_ZaZb = -9999.; eta_h1_ZaZb = -9999.; phi_h1_ZaZb = -9999.;
                   
                           h1 = Za + Zb;
                           mh1_ZaZb = h1.M();
                           pt_h1_ZaZb = h1.Pt();
                           eta_h1_ZaZb = h1.Eta();
                           phi_h1_ZaZb = h1.Phi();
            
                           h_mh1_ZaZb->Fill(mh1_ZaZb, wt);
                           h_pt_h1_ZaZb->Fill(pt_h1_ZaZb, wt);
                           h_eta_h1_ZaZb->Fill(eta_h1_ZaZb, wt);
                           h_phi_h1_ZaZb->Fill(phi_h1_ZaZb, wt);
		            
		            
		            
		         /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                           ^                                         ^ 
                           ^            Determine b1, b2             ^
                           ^                                         ^
                           ^              for h -> b1 b2             ^
                           ^                                         ^ 
                           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
       
                           bool found_bjet = false;
                           int nbjets = 0;            // total nb of b-jets found per event
                   
		           vector<Int_t> bjet_indx;    // saves index of jet tagged as b-jet 
                           bjet_indx.clear();
            
                          // Loop and select B-Tagged jets
                          for ( Int_t i = 0; i < Jet_size; i++){
					   
				// Examine which of jets are BTagged
                                UInt_t jet_bTag = Jet_BTag[i]; 
                       
				//std::cout << "Jet [" << i << "] BTag = " << jet_bTag << endl; 
             
                                if ( jet_bTag == 1) found_bjet = true; 
				 
                                    if ( found_bjet ){
						   
					  nbjets++;
				          bjet_indx.push_back(i);  // save index of jet tagged as b-jet 
				    } 	
		            } // end loop overall jets	
		 	       
		 	    cout << "========================================" << endl;
		 	    cout << " number of b-jets =  " << nbjets          << endl;
		            cout << "========================================" << endl;	
		           
		            if ( nbjets > 1 ){ // I have at least 2 b jets 
					   
				  double  mb1, pt_b1, eta_b1, phi_b1, mb2, pt_b2, eta_b2, phi_b2;
				       
				  // Initialize Variables
				  mb1 = -9999.; pt_b1 = -9999.; eta_b1 = -9999.; phi_b1 = -9999.; mb2 = -9999.; pt_b2 = -9999.; eta_b2 = -9999.; phi_b2 = -9999.;
			
			          Int_t b1_indx = bjet_indx[0];
			          Int_t b2_indx = bjet_indx[1];
			       
			           cout << "nbjets = " << nbjets << " b1-jet indx = " << b1_indx << " b2-jet indx = " << b2_indx << endl;
			           
			           b1.SetPtEtaPhiM(Jet_PT[b1_indx], Jet_Eta[b1_indx], Jet_Phi[b1_indx], Jet_Mass[b1_indx]);
			           b2.SetPtEtaPhiM(Jet_PT[b2_indx], Jet_Eta[b2_indx], Jet_Phi[b2_indx], Jet_Mass[b2_indx]);
						
			           mb1 = b1.M();
			           pt_b1 = b1.Pt();		    
                                   eta_b1 = b1.Eta();          
                                   phi_b1 = b1.Phi();           
                        
                                   mb2 = b2.M();
			           pt_b2 = b2.Pt();		    
                                   eta_b2 = b2.Eta();          
                                   phi_b2 = b2.Phi(); 
            
                                   double DR_b1b2; 
                       
                                   // Initialize Variable
                                   DR_b1b2 = -9999.;
            
                                   DR_b1b2 = TMath::Sqrt(TMath::Power((eta_b1 - eta_b2), 2) + TMath::Power((phi_b1 - phi_b2), 2));  
            
                                   h_mb_jet_1->Fill(mb1, wt);
                                   h_mb_jet_2->Fill(mb2, wt);
                                   h_pt_b_jet_1->Fill(pt_b1, wt);
                                   h_pt_b_jet_2->Fill(pt_b2, wt);
                                   h_eta_b_jet_1->Fill(eta_b1, wt);
                                   h_eta_b_jet_2->Fill(eta_b2, wt);
                                   h_phi_b_jet_1->Fill(phi_b1, wt);
                                   h_phi_b_jet_2->Fill(phi_b2, wt);
                                   h_DR_b1b2->Fill(DR_b1b2, wt); 
                         
                         
                                  //============================
		                  // Reconstruct h2 from b1, b2
		                  //============================
             
                                   double mh2_b1b2, pt_h2_b1b2, eta_h2_b1b2, phi_h2_b1b2;
                       
                                   // Initialize Variables
                                   mh2_b1b2 = -9999.; pt_h2_b1b2 = -9999.; eta_h2_b1b2 = -9999.; phi_h2_b1b2 = -9999.;
                        
                                   h2 = b1 +b2;
                                   mh2_b1b2 = h2.M();
                                   pt_h2_b1b2 = h2.Pt();
                                   eta_h2_b1b2 = h2.Eta();
                                   phi_h2_b1b2 = h2.Phi();
                
                                   h_mh1_b1b2->Fill(mh2_b1b2, wt);
                                   h_pt_h1_b1b2->Fill(pt_h2_b1b2, wt);
                                   h_eta_h1_b1b2->Fill(eta_h2_b1b2, wt);
                                   h_phi_h1_b1b2->Fill(phi_h2_b1b2, wt);
                
                
                      //======================================//
		              //                                       //
		              //      Reconstruct BSM H from SM h      // 
		              //               H -> h h                //
		              //                                       //
		              //=======================================//
                
                       double mH_hh, pt_H_hh, eta_H_hh, phi_H_hh;
               
                       // Initialize Variables
                       mH_hh = -9999.; pt_H_hh = -9999.; eta_H_hh = -9999.; phi_H_hh = -9999.;
                       
                       H = h1 + h2;
                       mH_hh = H.M(); 
                       pt_H_hh = H.Pt();
                       eta_H_hh = H.Eta();
                       phi_H_hh = H.Phi();
               
                       h_mh2_h1h1->Fill(mH_hh, wt);
                       h_pt_h2_h1h1->Fill(pt_H_hh, wt);
                       h_eta_h2_h1h1->Fill(eta_H_hh, wt);
                       h_phi_h2_h1h1->Fill(phi_H_hh, wt);
                
                        
		            } // end if nbjets > 1 
		            
		        } // end if mZb
		    } // end if mZa 
		  
        } // end if on 4 muons charge 
      } // end if MuonLoose_size > 3
    } // end loop over all muons  
			   
				  
 
   } // end loop overall events
   
   std::cout << "Loop Ends!" << endl; 
   std::cout << "Writing Histograms!" << endl;  
   
   h_muons_size_Loose->Write();
   h_pt_allMuons_Loose->Write();
   h_eta_allMuons_Loose->Write();
   h_phi_allMuons_Loose->Write();
   h_jets_size->Write();
   h_pt_allJets->Write();
   h_eta_allJets->Write();
   h_phi_allJets->Write();
   h_MET_size->Write();
   h_MET->Write();
   h_eta_MET->Write();
   h_eta_MET->Write();  
   h_DR_mu1mu2->Write();
   h_DR_mu3mu4->Write();
   h_DR_mu1mu3->Write();
   h_DR_mu2mu4->Write();
   h_DR_mu1mu4->Write();
   h_DR_mu2mu3->Write();
   h_mZ12->Write();
   h_mZ34->Write();
   h_mZ13->Write();
   h_mZ24->Write();
   h_mZ14->Write();
   h_mZ23->Write();
   h_mZa_4mu->Write();            
   h_mZb_4mu->Write();
   h_pt_Za->Write();
   h_pt_Zb->Write();
   h_eta_Za->Write();
   h_eta_Zb->Write();
   h_phi_Za->Write();
   h_phi_Zb->Write();
   h_mh1_ZaZb->Write();
   h_pt_h1_ZaZb->Write();
   h_eta_h1_ZaZb->Write();
   h_phi_h1_ZaZb->Write();
   h_mb_jet_1->Write();
   h_mb_jet_2->Write();
   h_pt_b_jet_1->Write();
   h_pt_b_jet_2->Write();
   h_eta_b_jet_1->Write();
   h_eta_b_jet_2->Write();
   h_phi_b_jet_1->Write();
   h_phi_b_jet_2->Write();
   h_DR_b1b2->Write(); 
   h_mh1_b1b2->Write();
   h_pt_h1_b1b2->Write();
   h_eta_h1_b1b2->Write();
   h_phi_h1_b1b2->Write();
   h_mh2_h1h1->Write();
   h_pt_h2_h1h1->Write();
   h_eta_h2_h1h1->Write();
   h_phi_h2_h1h1->Write();
   
   std::cout << "Writing Histograms ends!" << endl;
   
   op_file->Write();
   std::cout << "write output file! " << endl;
   
   op_file->Close();   
   
   std::cout << "saving..." << endl;    
   std::cout << "ROOT file: " << op_file << " has been created sucessfully!" << endl;
   
   
}
