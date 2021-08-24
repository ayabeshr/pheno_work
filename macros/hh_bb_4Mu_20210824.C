#define hh_bb_4Mu_20210824_cxx
#include "hh_bb_4Mu_20210824.h"
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

void hh_bb_4Mu_20210824::Loop()
{
//   In a ROOT session, you can do:
//      root> .L hh_bb_4Mu_20210824.C
//      root> hh_bb_4Mu_20210824 t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      
      // TH1D
      // various muon combinations bet 1,2,3,4 muons gives various Z masses, choose combination wich is nearer to Z mass 91 GeV to be Za, and the other with mass < 91 GeV take it to be Zb of signal and call them Za, Zb
   
      // size of available objects
      TH1D *h_muons_size;
      TH1D *h_jets_size;
      TH1D *h_MET_size;
   
      // MET    
      TH1D *h_MET;
      TH1D *h_eta_MET;
      TH1D *h_phi_MET;
   
      // All Muons   
      TH1D *h_pt_allMuons;
      TH1D *h_eta_allMuons;
      TH1D *h_phi_allMuons;
   
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
      //TH1D *h_phi_Za;
      //TH1D *h_phi_Zb;
      
      // h1 -> Za Zb
      TH1D *h_mh1_ZaZb;
      TH1D *h_pt_h1_ZaZb;
      TH1D *h_eta_h1_ZaZb;
      
      TH1D *h_mb_jet_1;
      TH1D *h_mb_jet_2;
      TH1D *h_pt_b_jet_1;
      TH1D *h_pt_b_jet_2;
      TH1D *h_eta_b_jet_1;
      TH1D *h_eta_b_jet_2;
      //TH1D *h_phi_b_jet_1;
      //TH1D *h_phi_b_jet_2;
      TH1D *h_DR_b1b2 ;
   
      // h1 -> b b~
      TH1D *h_mh1_b1b2;
      TH1D *h_pt_h1_b1b2;
      TH1D *h_eta_h1_b1b2;
   
      // h2 -> h1 h1 
      TH1D *h_mh2_h1h1;
      TH1D *h_pt_h2_h1h1;
      TH1D *h_eta_h2_h1h1;
      
      
      // Declarying Variables
      // those variables obtained by defining TLorentzVector objects used as indicators to get object pt, eta,...
   
      double mb1, mb2; 
      double pt_b1 , pt_b2; 
      double eta_b1, eta_b2; 
      double phi_b1, phi_b2; 
   
      double m_mu1, m_mu2, m_mu3, m_mu4;
      double pt_mu1, pt_mu2, pt_mu3, pt_mu4;
      double eta_mu1, eta_mu2, eta_mu3, eta_mu4;
      double phi_mu1, phi_mu2, phi_mu3, phi_mu4;
   
      double mZ12, mZ23, mZ34, mZ13, mZ14, mZ24;
      double dZ12, dZ23, dZ34, dZ13, dZ14, dZ24;
      double mZa, mZb, ptZa, ptZb, etaZa, etaZb; 
      double mh1_ZaZb, mh1_b1b2;
      double pt_h1_ZaZb, pt_h1_b1b2;
      double eta_h1_ZaZb, eta_h1_b1b2;
      double phi_h1_ZaZb, phi_h1_b1b2;
   
      double mh2_h1h1, pt_h2_h1h1, eta_h2_h1h1, phi_h2_h1h1;
      
      // Size of all Muons
      h_muons_size = new TH1D("h_muons_size", "h_muons_size", 10, 0., 10.);
      h_muons_size->GetXaxis()->SetTitle("Number of Muons");
      h_muons_size->GetYaxis()->SetTitle("Number of Events"); 
   
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
      h_pt_allMuons = new TH1D("h_pt_allMuons", "h_pt_allMuons", 100, 0., 100.);
      h_pt_allMuons->GetXaxis()->SetTitle("p_{T} [GeV/c]");
      h_pt_allMuons->GetYaxis()->SetTitle("Number of Events");
   
      // Eta of all Muons
      h_eta_allMuons = new TH1D("h_eta_allMuons", "h_eta_allMuons", 16, -8., 8.);
      h_eta_allMuons->GetXaxis()->SetTitle("#eta");
      h_eta_allMuons->GetYaxis()->SetTitle("Number of Events");
   
      // Phi of all Muons
      h_phi_allMuons= new TH1D("h_phi_allMuons", "h_phi_allMuons", 20, -10., 10.);
      h_phi_allMuons->GetXaxis()->SetTitle("#phi");
      h_phi_allMuons->GetYaxis()->SetTitle("Number of Events");

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
   
      h_pt_Za = new TH1D("h_pt_Za", "h_pt_Za", 100, 0., 100.); 
      h_pt_Za->GetXaxis()->SetTitle("p_{T} [GeV/c]");
      h_pt_Za->GetYaxis()->SetTitle("Number of Events");
   
      h_eta_Za = new TH1D("h_eta_Za", "h_eta_Za", 16, -8., 8.); 
      h_eta_Za->GetXaxis()->SetTitle("#eta");
      h_eta_Za->GetYaxis()->SetTitle("Number of Events");
   
    /*  h_phi_Za = new TH1D("h_phi_Za", "h_phi_Za", 20, -10., 10.); 
      h_phi_Za->GetXaxis()->SetTitle("#phi");
      h_phi_Za->GetYaxis()->SetTitle("Number of Events");
    */
      // Zb
   
      // Za Invariant Mass of dimuons not close to Z mass ~ 91. GeV 
      h_mZb_4mu = new TH1D("h_mZb_4mu", "h_mZb_4mu_Not Close to Z Mass", 120, 0., 120.);
      h_mZb_4mu->GetXaxis()->SetTitle("Invariant mass of dimuon [GeV/c^{2}]"); 
      h_mZb_4mu->GetYaxis()->SetTitle("Number of Events");
   
      h_pt_Zb = new TH1D("h_pt_Zb", "h_pt_Zb", 100, 0., 100.); 
      h_pt_Zb->GetXaxis()->SetTitle("p_{T} [GeV/c]");
      h_pt_Zb->GetYaxis()->SetTitle("Number of Events");
   
      h_eta_Zb = new TH1D("h_eta_Zb", "h_eta_Zb", 16, -8., 8.); 
      h_eta_Zb->GetXaxis()->SetTitle("#eta");
      h_eta_Zb->GetYaxis()->SetTitle("Number of Events");
   
     /* h_phi_Zb = new TH1D("h_phi_Zb", "h_phi_Zb", 20, -10., 10.); 
      h_phi_Zb->GetXaxis()->SetTitle("#phi");
      h_phi_Zb->GetYaxis()->SetTitle("Number of Events");
     */
   
      h_pt_b_jet_1 = new TH1D("h_pt_b_jet_1", "h_pt_b_jet_1", 100, 0., 100.);
      h_pt_b_jet_1->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      h_pt_b_jet_1->GetYaxis()->SetTitle("Number of Events");
   
      h_pt_b_jet_2 = new TH1D("h_pt_b_jet_2", "h_pt_b_jet_2", 100, 0., 100.);
      h_pt_b_jet_2->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      h_pt_b_jet_2->GetYaxis()->SetTitle("Number of Events");
   
      h_eta_b_jet_1 = new TH1D("h_eta_b_jet_1", "h_eta_b_jet_1", 20, -8., 8.);
      h_eta_b_jet_1->GetXaxis()->SetTitle("#eta");
      h_eta_b_jet_1->GetYaxis()->SetTitle("Number of Events");
   
      h_eta_b_jet_2 = new TH1D("h_eta_b_jet_2", "h_eta_b_jet_2", 20, -8., 8.);
      h_eta_b_jet_2->GetXaxis()->SetTitle("#eta");
      h_eta_b_jet_2->GetYaxis()->SetTitle("Number of Events");
   
    /*  h_phi_b_jet_1 = new TH1D("h_phi_b_jet_1", "h_phi_b_jet_1", 20, -10., 10.);
      h_phi_b_jet_1->GetXaxis()->SetTitle("#phi");
      h_phi_b_jet_1->GetYaxis()->SetTitle("Number of Events");
   
      h_phi_b_jet_2 = new TH1D("h_phi_b_jet_2", "h_phi_b_jet_2", 20, -10., 10.);
      h_phi_b_jet_2->GetXaxis()->SetTitle("#phi");
      h_phi_b_jet_2->GetYaxis()->SetTitle("Number of Events");
     */
   
      // DR seperation bet 4Muons  
      h_DR_mu1mu2 = new TH1D("h_DR_mu1mu2", "DR between muon1 muon2", 20, 0., 20.);
      h_DR_mu1mu2->GetXaxis()->SetTitle("#Delta R");
      h_DR_mu1mu2->GetYaxis()->SetTitle("Number of Events");
      
      h_DR_mu3mu4 = new TH1D("h_DR_mu3mu4", "DR between muon3 muon4", 20, 0., 10.);
      h_DR_mu3mu4->GetXaxis()->SetTitle("#Delta R");
      h_DR_mu3mu4->GetYaxis()->SetTitle("Number of Events");
      
      h_DR_mu1mu3 = new TH1D("h_DR_mu1mu3", "DR between muon1 muon3", 20, 0., 10.);
      h_DR_mu1mu3->GetXaxis()->SetTitle("#Delta R");
      h_DR_mu1mu3->GetYaxis()->SetTitle("Number of Events");
      
      h_DR_mu2mu4 = new TH1D("h_DR_mu2mu4", "DR between muon2 muon4", 20, 0., 10.);
      h_DR_mu2mu4->GetXaxis()->SetTitle("#Delta R");
      h_DR_mu2mu4->GetYaxis()->SetTitle("Number of Events");
   
      h_DR_mu1mu4 = new TH1D("h_DR_mu1mu4", "DR between muon1 muon4", 20, 0., 10.);
      h_DR_mu1mu4->GetXaxis()->SetTitle("#Delta R");
      h_DR_mu1mu4->GetYaxis()->SetTitle("Number of Events");
   
      h_DR_mu2mu3 = new TH1D("h_DR_mu2mu3", "DR between muon2 muon3", 20, 0., 10.);
      h_DR_mu2mu3->GetXaxis()->SetTitle("#Delta R");
      h_DR_mu2mu3->GetYaxis()->SetTitle("Number of Events");
   
   
      // DR seperation bet 2 b jets 
      h_DR_b1b2 = new TH1D("h_DR_b1b2", "h_DR_b1b2", 20, 0., 10.); 
      h_DR_b1b2->GetXaxis()->SetTitle("#Delta R");
      h_DR_b1b2->GetYaxis()->SetTitle("Number of Events");
   
      // h1 -> Za Zb
      h_mh1_ZaZb = new TH1D("h_mh1_ZaZb", "h_mh1_ZaZb", 100, 0., 200.);
      h_mh1_ZaZb->GetXaxis()->SetTitle("Invariant mass of 4 muons (GeV/c^{2})");
      h_mh1_ZaZb->GetYaxis()->SetTitle("Number of Events");
   
      h_pt_h1_ZaZb = new TH1D("h_pt_h1_ZaZb", "h_pt_h1_ZaZb", 100, 0., 100.);
      h_pt_h1_ZaZb->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      h_pt_h1_ZaZb->GetYaxis()->SetTitle("Number of Events");
   
      h_eta_h1_ZaZb = new TH1D("h_eta_h1_ZaZb", "h_eta_h1_ZaZb", 16, -8., 8.);
      h_eta_h1_ZaZb->GetXaxis()->SetTitle("#eta");
      h_eta_h1_ZaZb->GetYaxis()->SetTitle("Number of Events");
   
      // h1 -> b b~
      h_mh1_b1b2 = new TH1D("h_mh1_b1b2", "h_mh1_b1b2", 100, 0., 200.);
      h_mh1_b1b2->GetXaxis()->SetTitle("Invariant mass of 2 b-jets (GeV/c^{2})");
      h_mh1_b1b2->GetYaxis()->SetTitle("Number of Events");
   
      h_pt_h1_b1b2 = new TH1D("h_pt_h1_b1b2", "h_pt_h1_b1b2", 100, 0., 100.);
      h_pt_h1_b1b2->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      h_pt_h1_b1b2->GetYaxis()->SetTitle("Number of Events");
   
      h_eta_h1_b1b2 = new TH1D("h_eta_h1_b1b2", "h_eta_h1_b1b2", 16, -8., 8.);
      h_eta_h1_b1b2->GetXaxis()->SetTitle("#eta");
      h_eta_h1_b1b2->GetYaxis()->SetTitle("Number of Events"); 
   
      // h2 -> h1 h1 
      h_mh2_h1h1 = new TH1D("h_mh2_h1h1", "h_mh2_h1h1", 350, 0., 700.);
      h_mh2_h1h1->GetXaxis()->SetTitle("Invariant mass of 2 b-jets + 4 muons (GeV/c^{2})");
      h_mh2_h1h1->GetYaxis()->SetTitle("Number of Events"); 
   
      h_pt_h2_h1h1 = new TH1D("h_pt_h2_h1h1", "h_pt_h2_h1h1", 100, 0., 100.);
      h_pt_h2_h1h1->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      h_pt_h2_h1h1->GetYaxis()->SetTitle("Number of Events");
   
      h_eta_h2_h1h1 = new TH1D("h_eta_h2_h1h1", "h_eta_h2_h1h1", 16, -8., 8.);
      h_eta_h2_h1h1->GetXaxis()->SetTitle("#eta");
      h_eta_h2_h1h1->GetYaxis()->SetTitle("Number of Events");
      
      
      // Datasets NoPU 
      TFile *input_file = new TFile("/home/Aya/DY_BG_14TeV_SMfull_GEN/DY_BG_14TeV_SMfull_pythia8_CMSPhaseII-0PU_GEN-SIM.root", "READ");
   
      // Results root file
      TFile *output_file = new TFile("/home/Aya/out_h2h1h1_bb4Mu_PhaseII_0PU.root");
      
      
      //------------------------WEIGHT Calculation---------------------------
  
      float Lumi_data = 3.e+3;    // in 1/fb
      //Lumi_mc = nEvents/xsection(fb);
  
     // Lumi_mc for process: DY at 14TeV 
     float Lumi_mc = 1000000./898225.; 
     float wt = Lumi_data/Lumi_mc;
      
      
      
      
      std::cout << "just a test!" << endl; 
   }
}
