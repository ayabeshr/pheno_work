////////////////////////////////////////////////////////
//      To run macro in a ROOT session do:            //
//     gSystem->Load("/path/to/libDelphes.so");       //
//     .L path/to/DemoAnalyzer.C                      //
//     DemoAnalyzer e;                                //
//     e.Loop();                                      //
////////////////////////////////////////////////////////


#define DemoAnalyzer_cxx

R__LOAD_LIBRARY(libPhysics)
//R__LOAD_LIBRARY(/home/aya/programs/Delphes-3.4.2/libDelphes.so)
R__ADD_LIBRARY_PATH(/home/aya/programs/Delphes-3.4.2)
R__LOAD_LIBRARY(libDelphes)

//#include "/HEP_DATA/aya/DemoAnalyzer.h"
#include "/home/aya/Pheno_Work/analysis/macro/checkChargePairs_1234_1324/DemoAnalyzer.h"
#include "/home/aya/Pheno_Work/analysis/macro/checkChargePairs_1234_1324/newtree_leptonsadded.h"
//#include "/home/aya/Pheno_Work/analysis/macro/cuts_observable.C"
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
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TFormula.h"
#include <TMath.h>
#include "THStack.h"
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>



void DemoAnalyzer::Loop()
{

   if (fChain == 0) return;
   
   
   // Array act as a counter stores how many events passed a certain selection, 
   Int_t NEvents[19]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};    // size = nb of selection cuts, each element referes to certain selection  
   
   
   //================================================================================================//
   //                                   DATASET SIGNAL H->hh->bb4Mu  NoPU                            //
   //================================================================================================//  
   
   // TFile *in_sig = new TFile("/media/aya/PACKUP/Aya/final_analysis/inputs/gghhbb4Mu_final.root", "READ");
   //TFile *in_sig = new TFile("/media/aya/LinuxSpace/Pheno_Work_2/analysis/gghhbb4M_AllDiagrams_20210915.root", "READ");
    TFile *in_sig = new TFile("/home/aya/programs/Delphes-3.4.2/sig_1e04Events.root", "READ");   //xsec order 10-4 pb 
   // TFile *in_sig = new TFile("/home/aya/programs/Delphes-3.4.2/sig3_10kEvents_xsecOrder-5_pb.root", "READ");
   // TFile *in_sig = new TFile("/home/aya/programs/Delphes-3.4.2/sig3_10kEvents_xsecOrder-6_pb.root", "READ");
   
   //================================================================================================//
   //                                   DATASET Background 14 TeV                                    //
   //================================================================================================// 

   // TFile *TT_CMS_Snowmass = new TFile("/media/aya/PACKUP/Aya/final_analysis/inputs/TT_TuneCUETP8M2T4_14TeV-powheg-pythia8_1751_2.root", "READ");
   // TFile* in_ttbar    = new TFile("/media/aya/PACKUP/Aya/final_analysis/inputs/ttbar.root"    , "READ");
   //  TFile *in_ZZ4Mu    = new TFile("/media/aya/PACKUP/Aya/final_analysis/inputs/ZZ4Mu.root"    , "READ");
   //  TFile *in_DY       = new TFile("/media/aya/PACKUP/Aya/final_analysis/inputs/DrellYan.root"       , "READ");   // haven't produced yet DY.root
   // TFile *in_ZZbb2Mu  = new TFile("/media/aya/PACKUP/Aya/final_analysis/inputs/ZZbb2Mu.root"  , "READ");
   //  TFile *in_ZWpbbMuNL= new TFile("/media/aya/PACKUP/Aya/final_analysis/inputs/ZWpbbMuNL.root", "READ");
   // TFile *in_ZZZbb4Mu = new TFile("/media/aya/PACKUP/Aya/final_analysis/inputs/ZZZbb4Mu.root" , "READ");
   // TFile *in_smHTobb    = new TFile("/media/aya/PACKUP/Aya/final_analysis/inputs/sm_ggToHTobb_NLO.root"   , "READ");
   // TFile *in_sm_H4Mu = new TFile("/media/aya/PACKUP/Aya/final_analysis/inputs/H4MU.root", "READ");
   
   //================================================================================================//
   //                                         Output Root files                                      //
   //================================================================================================//
   
   // TFile* out_sig    = new TFile("/media/aya/PACKUP/Aya/final_analysis/outputs/output_demo_ggHhhbb4M_final.root"   , "RECREATE");   // loose b-jets
   // TFile* out_sig    = new TFile("/home/aya/Pheno_Work/analysis/output_demo_ggHhhbb4M_final.root"   , "RECREATE");               // tight b-jets
   // TFile* out_sig    = new TFile("/media/aya/PACKUP/Aya/final_analysis/outputs/output_demo_ggHhhbb4M_2_final.root"   , "RECREATE"); 
   // TFile *out_sig = new TFile("/home/aya/Pheno_Work/analysis/out_sig2_1e04Events.root", "RECREATE");
   // TFile *out_sig = new TFile("/home/aya/Pheno_Work/analysis/out_sig3_1e04Events_unweighted_v2.root", "RECREATE");
    //  TFile *out_sig = new TFile("/home/aya/Pheno_Work/analysis/output_MuCut_TestZmass/out_sig3_1e04Events_unweighted_MuCUTs_TestZmass.root", "RECREATE");
    // TFile *out_sig = new TFile("/home/aya/Pheno_Work/analysis/output_FullSel_cutflow2/out2_sig3_1e04Events_unweighted_FullSel_cutflow2.root", "RECREATE");
    TFile *out_sig = new TFile("/home/aya/Pheno_Work/analysis/signal_20220420.root", "RECREATE");
   //TFile *out_sig = new TFile("/home/aya/Pheno_Work/analysis/out_sig3_10kEvents_xsecOrder-5_pb_unweighted_v2.root", "RECREATE");
   // TFile *out_sig = new TFile("/home/aya/Pheno_Work/analysis/out_sig3_10kEvents_xsecOrder-6_pb_unweighted_v2.root", "RECREATE");

    // TFile *out_sig = new TFile("/home/aya/Pheno_Work/analysis/out_sig3_1e04Events_unweighted_withoutDR_blep.root", "RECREATE");

    // TFile *out_TT_CMS_Snowmass = new TFile("/media/aya/PACKUP/Aya/final_analysis/outputs/out_TT_TuneCUETP8M2T4_14TeV-powheg-pythia8_1751_2.root", "RECREATE");
    // TFile* out_ttbar    = new TFile("/media/aya/PACKUP/Aya/final_analysis/outputs/output_demo_ttbar.root"             , "RECREATE");
   //TFile* out_ttbar    = new TFile("/media/aya/PACKUP/Aya/final_analysis/outputs/output_demo_ttbar_unweighted_v2.root"             , "RECREATE");
     //  TFile* out_ttbar    = new TFile("/home/aya/Pheno_Work/analysis/output_MuCut_TestZmass/output_demo_ttbar_unweighted_MuCUTs_TestZmass.root"             , "RECREATE");
     // TFile *out_ttbar = new TFile("/home/aya/Pheno_Work/analysis/output_FullSel_cutflow2/out_ttbar_unweighted_FullSel_cutflow2.root", "RECREATE");
   // TFile* out_ttbar    = new TFile("/media/aya/PACKUP/Aya/final_analysis/outputs/output_demo_ttbar_unweighted_tight_bjets.root"             , "RECREATE");
   // TFile* out_ZZ4Mu= new TFile("/media/aya/PACKUP/Aya/final_analysis/outputs/output_demo_ZZ4Mu.root"             , "RECREATE");
   // TFile* out_ZZ4Mu = new TFile("/media/aya/PACKUP/Aya/final_analysis/outputs/output_demo_ZZ4Mu_unweighted.root"             , "RECREATE");
   //TFile* out_ZZ4Mu = new TFile("/home/aya/Pheno_Work/analysis/output_MuCut_TestZmass/output_demo_ZZ4Mu_unweighted_MuCUTs_TestZmass.root"             , "RECREATE");
     //   TFile *out_ZZ4Mu = new TFile("/home/aya/Pheno_Work/analysis/output_FullSel_cutflow2/out_ZZ4Mu_unweighted_FullSel_cutflow2.root", "RECREATE");
   // TFile* out_DY         = new TFile("/media/aya/PACKUP/Aya/final_analysis/outputs/output_demo_DY_unweighted.root"                , "RECREATE");
   // TFile* out_DY         = new TFile("/home/aya/Pheno_Work/analysis/output_MuCut_TestZmass/output_demo_DY_unweighted_MuCUTs_TestZmass.root"                , "RECREATE");
   //  TFile *out_DY = new TFile("/home/aya/Pheno_Work/analysis/output_FullSel_cutflow2/out_DY_unweighted_FullSel_cutflow2.root", "RECREATE");
   //TFile* out_ZZbb2Mu    = new TFile("/media/aya/PACKUP/Aya/final_analysis/outputs/output_demo_ZZbb2Mu.root"           , "RECREATE");
    // TFile* out_ZZbb2Mu    = new TFile("/media/aya/PACKUP/Aya/final_analysis/outputs/output_demo_ZZbb2Mu_unweighted.root"           , "RECREATE");
  // TFile* out_ZZbb2Mu = new TFile("/home/aya/Pheno_Work/analysis/output_MuCut_TestZmass/output_demo_ZZbb2Mu_unweighted_MuCUTs_TestZmass.root"             , "RECREATE");
  //  TFile *out_ZZbb2Mu = new TFile("/home/aya/Pheno_Work/analysis/output_FullSel_cutflow2/out_ZZbb2Mu_unweighted_FullSel_cutflow2.root", "RECREATE");
    //TFile* ZZZbb4Mu   = new TFile("/media/aya/PACKUP/Aya/final_analysis/outputs/output_demo_ZZZbb4Mu.root"          , "RECREATE");
    // TFile* out_ZZZbb4Mu = new TFile("/home/aya/Pheno_Work/analysis/output_MuCut_TestZmass/output_demo_ZZZbb4Mu_unweighted_MuCUTs_TestZmass.root"             , "RECREATE");
    // TFile *out_ZZZbb4Mu = new TFile("/home/aya/Pheno_Work/analysis/output_FullSel_cutflow2/out_ZZZbb4Mu_unweighted_FullSel_cutflow2.root", "RECREATE");
    //  TFile* out_ZWpbbMuNL  = new TFile("/home/aya/Pheno_Work/analysis/output_MuCut_TestZmass/output_demo_ZWpbbMuNL_unweighted_MuCUTs_TestZmass.root"         , "RECREATE");
    //TFile *out_ZWpbbMuNL = new TFile("/home/aya/Pheno_Work/analysis/output_FullSel_cutflow2/out_ZWpbbMuNL_unweighted_FullSel_cutflow2.root", "RECREATE");
    /* TFile* ZWp2MuMuNL = new TFile("/media/aya/LinuxSpace/Pheno_Work_2/final_analysis/output_demo_ZWpbbMuNL.root"    , "RECREATE");
   TFile* ZWmbbMuNL  = new TFile("/media/aya/LinuxSpace/Pheno_Work_2/final_analysis/output_demo_ZWmbbMuNL.root"    , "RECREATE");
   TFile* ZWm2MuMuNL = new TFile("/media/aya/LinuxSpace/Pheno_Work_2/final_analysis/output_demo_ZWmbbMuNL.root"    , "RECREATE"); 
   TFile* WW2Mu2NL   = new TFile("/media/aya/LinuxSpace/Pheno_Work_2/final_analysis/output_demo_WW2Mu2NL.root"     , "RECREATE"); */
   //  TFile* out_sm_HTobb      = new TFile("/home/aya/Pheno_Work/analysis/output_MuCut_TestZmass/output_smHbb_unweighted_MuCUTs_TestZmass.root"        , "RECREATE");
   // TFile *out_sm_HTobb = new TFile("/home/aya/Pheno_Work/analysis/output_FullSel_cutflow2/out_sm_HTobb_unweighted_FullSel_cutflow2.root", "RECREATE");
   // TFile* out_sm_H4Mu = new TFile("/media/aya/PACKUP/Aya/final_analysis/outputs/output_sm_H4Mu_unweighted.root","RECREATE");
   // TFile* out_sm_H4Mu = new TFile("/home/aya/Pheno_Work/analysis/output_MuCut_TestZmass/output_sm_H4Mu_unweighted_MuCUTs_TestZmass.root","RECREATE");
   // TFile *out_sm_H4Mu = new TFile("/home/aya/Pheno_Work/analysis/output_FullSel_cutflow2/out_sm_H4Mu_unweighted_FullSel_cutflow2.root", "RECREATE");

   //------------------------WEIGHT Calculation---------------------------
  
   float Lumi_data = 3.e+03;    // in 1/fb
	
   //-----------------------------------//
   //  Lumi_mc = nEvents/xsection(fb);  //
   //-----------------------------------//	

   int scale_factor = 1.e+02;    // to be consistent with theory
   float org_xsc_pb = 1.2416e-04;    //original xsec in pb
   float org_xsc_fb = org_xsc_pb * 1.e+03;
   float sf_xsc_pb = org_xsc_pb * scale_factor;
   float sf_xsc_fb = sf_xsc_pb * 1.e+03;   // scaled xsec in fb 
     
   //float Lumi_mc = 1.e+04 / sf_xsc_fb;
   float Lumi_mc = 1.e+04 / org_xsc_fb;
   //float Lumi_mc = 1.e+05/ 3.554e-05;       // signal: gghhbb4Mu     xsec should be at least of oder 10 to have at least 1 event to be seen or 10^2 to have ~ 10 events
   // float Lumi_mc = 1.e+05/ 1.247e-01;         // original xsec in fb for highest xsec found in model using scan16 
   //float Lumi_mc = 1.e+06/898225.;          // DY
   // float Lumi_mc = 1.e+06/7559.17;            // ttbar
   // float Lumi_mc = 1.e+06/ 12.16;           // ZZ4M  
   //float Lumi_mc = 1.e+06/106.6518;         // ZZbb2Mu
   // float Lumi_mc = 1.e+06/5.519569e-03;    // ZZZbb4Mu
   //float Lumi_mc = 1.e+06/280.1;            // ZWpbbMuNL
   // float Lumi_mc = 1.e+06/10870;              // sm_HTobb 
   // xsec for smHTo4Mu is :  0.3351 fb 
   // float wt = Lumi_data/Lumi_mc;
   //float wt = 1.;    // examine plots with unweighted events
   
   
   //------------------------END WEIGHT CALC----------------------------


   
   //---------------------Define 1D Vectors--------------------

			  
   vector<Int_t> v_muon_idx;            // saves muon index if fullfill Object Selection

   // define vectors for 4 muons of signal
   vector<Float_t> lep1;
   vector<Float_t> lep2;
   vector<Float_t> lep3;
   vector<Float_t> lep4;


   vector<Int_t> bjet_indx ;            // saves tight b-jet index 
   vector<Int_t> bjet_indx_AfterSel; 	// saves index of tight b-jets after applying series of event selections on them
  
   vector<Int_t> v_b1_jet_Charge;
   vector<Int_t> v_b2_jet_Charge;

   vector<Float_t> dr_b1b2;
   vector<Float_t> dr_h1h2;

   int NEvents_13;


   //---------------------Define BOOLEANS---------------------------
   bool axpEvent;
   bool isCase1;
   bool isCase2;
   bool catgr1;    // category1 for mu comb 1234
   bool catgr2;    // category2 for mu comb 1423
   bool catgr3;    // category3 for mu comb 1324

  
   //---------------------START Defining Tree Branches--------------------
   
 
   // b_Gen_Muon           = new_tree->Branch("Gen_Muons"           ,   &genMuons  );
   // b_Gen_Bjets          = new_tree->Branch("Gen_BJets"           ,   &gen_bjet  );
   // b_LooseMuons         = new_tree->Branch("LooseMuons"          ,   &loose_mu  );
   // b_deltaR_muons_loose = new_tree->Branch("DeltaR_LooseMuons"   ,   &drMuonsLoose   );
   b_TightMuons         = new_tree->Branch("TightMuons"          ,   &tight_mu  );
   b_muon1              = new_tree->Branch("Muon1"               ,   &muon1  );
   b_muon2              = new_tree->Branch("Muon2"               ,   &muon2  );
   b_muon3              = new_tree->Branch("Muon3"               ,   &muon3  );
   b_muon4              = new_tree->Branch("Muon4"               ,   &muon4  );
   b_ZcombMass          = new_tree->Branch("ZcombMass"           ,   &mz );
   // b_deltaR_muons_tight = new_tree->Branch("DeltaR_TightMuons"   ,   &drMuonsTight  );
   // b_Jets               = new_tree->Branch("Jets"                ,   &jets      );
   //b_DR_LooseBJets_4Mu  = new_tree->Branch("DR_LooseBJets_4Muons",   &drLooseBjet4Mu      );
   // b_DR_MedBJets_4Mu    = new_tree->Branch("DR_MedBJets_4Muons"  ,   &drMedBjet4Mu      );
   // b_DR_TightBJets_4Mu  = new_tree->Branch("DR_TightBJets_4Muons",   &drTightBjet4Mu      );
   // b_MET                = new_tree->Branch("MET"                 ,   &met       );
   // b_4Muons_beforeCut   = new_tree->Branch("Four_Muons_beforeCut",   &four_muons_beforeCut);
   // b_deltaR_muons       = new_tree->Branch("DeltaR_Muons"        ,   &drMuons   );
   // b_dR_bjet_mu_NoCUT   = new_tree->Branch("DR_bjet_mu_NoCUT"    ,   &dR_bjets_mu_NoCUT); 
   // b_4Muons_FullSel     = new_tree->Branch("4Muons_after_FullSel",   &four_muons_afterFullSel);
   // b_dR_4Muons_FullSel  = new_tree->Branch("DR_4Muons_after_FullSel",&drMuons_FullSel);
   b_Za_signal          = new_tree->Branch("ZaOFSignal"          ,   &za_sig    );
   b_Zb_signal          = new_tree->Branch("ZbOFSignal"          ,   &zb_sig    );
   /* b_b1_signal          = new_tree->Branch("bjet_1_OFSignal"     ,   &b1_sig    );
   b_b2_signal          = new_tree->Branch("bjet_2_OFSignal"     ,   &b2_sig    );
   b_1st_Higgs_Ofsignal = new_tree->Branch("1st_Higgs_OFSignal"  ,   &smhiggs1  );
   b_2nd_Higgs_Ofsignal = new_tree->Branch("2nd_Higgs_OFSignal"  ,   &smhiggs2  );
   b_BSM_Higgs_Ofsignal = new_tree->Branch("Heavy_Higgs_OFSignal",   &heavyHiggs);
   b_dr_b1b2            = new_tree->Branch("DR_2bjets_signal"    ,   &dr_b1b2,   "dr_b1b2/F");
   b_dr_h1h2            = new_tree->Branch("DR_2SM_HIGSS"        ,   &dr_h1h2,   "dr_h1h2/F");
   */
   
   //---------------------END Defining Tree Branches--------------------



  //---------------------START Defining 1D&2D HISTOS--------------------

  // TH1F *no_loose_muons_NOCUT = new TH1F("NbLooseMuons_nocut", "Number LooseMuons_nocut", 50, 0, 50);
   //  no_loose_muons_NOCUT->GetXaxis()->SetTitle("number");
   // no_loose_muons_NOCUT->GetYaxis()->SetTitle("Events");

   /* TH1F *hh_dr = new TH1F("dr_hh", "DR H4Mu, Hbb", 20, 0, 10);
   hh_dr->GetXaxis()->SetTitle("#Delta R");
   hh_dr->GetYaxis()->SetTitle("Events");

   TH1F *b1b2_dr = new TH1F("dr_b1b2", "DR b1,b2", 20, 0, 10);
   b1b2_dr->GetXaxis()->SetTitle("#Delta R");
   b1b2_dr->GetYaxis()->SetTitle("Events");


   TH1F *H_mass = new TH1F("H_mass", "BSM Heavy Higgs Mass", 700, 100, 2000);
   H_mass->GetXaxis()->SetTitle("m_{4mu+jj} [GeV]");
   H_mass->GetYaxis()->SetTitle("Events");
   */

   TH2F *ZZ = new TH2F("Z1Z2_masses", "Z1 & Z2 mass in GeV", 50,0,100, 60,60,140);
   ZZ->GetXaxis()->SetTitle("Z2 mass [GeV]");
   ZZ->GetYaxis()->SetTitle("Z1 mass [GeV]");

   TH2F *hh = new TH2F("hh_masses", "SM hh masses in GeV", 60,60,140, 60,60,140);
   hh->GetXaxis()->SetTitle("m_{4muons} [GeV]");
   hh->GetYaxis()->SetTitle("m_{bb} [GeV]");

   /* TH1F *H_pt = new TH1F("H_pt", "BSM Heavy Higgs p_{T}", 500, 0, 1000);
   H_pt->GetXaxis()->SetTitle("p_{T} [GeV/C]");
   H_pt->GetYaxis()->SetTitle("Events");

   TH1F *H_eta = new TH1F("H_eta", "BSM Heavy Higgs Eta", 20, -10, 10);
   H_eta->GetXaxis()->SetTitle("#eta");
   H_eta->GetYaxis()->SetTitle("Events");

   TH1F *H_phi = new TH1F("H_phi", "BSM Heavy Higgs Phi", 20, -10, 10);
   H_phi->GetXaxis()->SetTitle("#phi");
   H_phi->GetYaxis()->SetTitle("Events");
   */
   
   //---------------------END Defining 1D HISTOS--------------------


   
   
   
   
  /*===================================================================================*/  
  /*------------------------------Looping over ALL Events-----------------------------*/
  /*==================================================================================*/ 
   
   
   Long64_t nentries = fChain->GetEntriesFast();
   // Long64_t nentries = 10;
    cout << "nEvents = " << nentries << endl;

  
   Long64_t nbytes = 0, nb = 0;
   
   for (Long64_t jentry = 0; jentry < nentries; jentry++)  
   {
      cout << "******START EVENT LOOP!******    ,    Event nb = " << jentry << endl; 
	  
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      cout << "nbytes for current entry = " << nbytes << endl;
      // if (Cut(ientry) < 0) continue;

      NEvents_13 = 0;

      //============Initialize BOOLEANS============
      axpEvent = false;
      isCase1  = false;
      isCase2  = false; 
      catgr1   = false;
      catgr2   = false;
      catgr3   = false;
     
      //============Initialize vectors=============

      v_muon_idx.clear();

      lep1.clear();
      lep2.clear();
      lep3.clear();
      lep4.clear();

      muon1.pT.clear();
      muon1.eta.clear();
      muon1.phi.clear();
      muon1.mass.clear();

      muon2.pT.clear();
      muon2.eta.clear();
      muon2.phi.clear();
      muon2.mass.clear();

      muon3.pT.clear();
      muon3.eta.clear();
      muon3.phi.clear();
      muon3.mass.clear();

      muon4.pT.clear();
      muon4.eta.clear();
      muon4.phi.clear();
      muon4.mass.clear();

      mz.mZ_12.clear();
      mz.mZ_34.clear();
      mz.mZ_13.clear();
      mz.mZ_24.clear();
      mz.mZ_14.clear();
      mz.mZ_23.clear();

      
      bjet_indx.clear();
      bjet_indx_AfterSel.clear();

      v_b1_jet_Charge.clear();
      v_b2_jet_Charge.clear();
      
      /* genMuons.gen_muon_pt.clear();
      genMuons.gen_muon_eta.clear();
      genMuons.gen_muon_phi.clear();
      
      gen_bjet.gen_bjet_pt.clear();
      gen_bjet.gen_bjet_eta.clear();
      gen_bjet.gen_bjet_phi.clear(); 

      loose_mu.lmu_pt.clear();
      loose_mu.lmu_eta.clear();
      loose_mu.lmu_phi.clear();
      loose_mu.lmu_size.clear();

      loose_mu.muon1_pt.clear();
      loose_mu.muon2_pt.clear();
      loose_mu.muon3_pt.clear();
      loose_mu.muon4_pt.clear();
      loose_mu.muon1_eta.clear();
      loose_mu.muon2_eta.clear();
      loose_mu.muon3_eta.clear();
      loose_mu.muon4_eta.clear();
      loose_mu.muon1_phi.clear();
      loose_mu.muon2_phi.clear();
      loose_mu.muon3_phi.clear();
      loose_mu.muon4_phi.clear();
      loose_mu.dxy_mu1.clear();
      loose_mu.dxy_mu2.clear();
      loose_mu.dxy_mu3.clear();
      loose_mu.dxy_mu4.clear();

      loose_mu.loose_mu12_mass.clear();
      loose_mu.loose_mu34_mass.clear();
      loose_mu.loose_mu13_mass.clear();
      loose_mu.loose_mu24_mass.clear();
      loose_mu.loose_mu14_mass.clear();
      loose_mu.loose_mu23_mass.clear();

      drMuonsLoose.v_DR_mu1mu2.clear();
      drMuonsLoose.v_DR_mu3mu4.clear();
      drMuonsLoose.v_DR_mu1mu3.clear();
      drMuonsLoose.v_DR_mu2mu4.clear();
      drMuonsLoose.v_DR_mu1mu4.clear();
      drMuonsLoose.v_DR_mu2mu3.clear();  */

      tight_mu.tmu_pt.clear();
      tight_mu.tmu_eta.clear();
      tight_mu.tmu_phi.clear();
      tight_mu.tmu_size.clear();

      /*  tight_mu.muon1_pt.clear();
      tight_mu.muon2_pt.clear();
      tight_mu.muon3_pt.clear();
      tight_mu.muon4_pt.clear();
      tight_mu.muon1_eta.clear();
      tight_mu.muon2_eta.clear();
      tight_mu.muon3_eta.clear();
      tight_mu.muon4_eta.clear();
      tight_mu.muon1_phi.clear();
      tight_mu.muon2_phi.clear();
      tight_mu.muon3_phi.clear();
      tight_mu.muon4_phi.clear();
      tight_mu.dxy_mu1.clear();
      tight_mu.dxy_mu2.clear();
      tight_mu.dxy_mu3.clear();
      tight_mu.dxy_mu4.clear();
      tight_mu.dz_mu1.clear();
      tight_mu.dz_mu2.clear();
      tight_mu.dz_mu3.clear();
      tight_mu.dz_mu4.clear(); 

      tight_mu.tight_mu12_mass.clear();
      tight_mu.tight_mu34_mass.clear();
      tight_mu.tight_mu13_mass.clear();
      tight_mu.tight_mu24_mass.clear();
      tight_mu.tight_mu14_mass.clear();
      tight_mu.tight_mu23_mass.clear();

      drMuonsTight.v_DR_mu1mu2.clear();
      drMuonsTight.v_DR_mu3mu4.clear();
      drMuonsTight.v_DR_mu1mu3.clear();
      drMuonsTight.v_DR_mu2mu4.clear();
      drMuonsTight.v_DR_mu1mu4.clear();
      drMuonsTight.v_DR_mu2mu3.clear(); 

      jets.jet_pt.clear();
      jets.jet_eta.clear();
      jets.jet_phi.clear();
      jets.jet_size.clear();

      jets.loose_bjet_pt.clear();
      jets.loose_bjet_eta.clear();
      jets.loose_bjet_phi.clear();
      jets.loose_bjet_size.clear();

      jets.med_bjet_pt.clear();
      jets.med_bjet_eta.clear();
      jets.med_bjet_phi.clear();
      jets.med_bjet_size.clear();

      jets.tight_bjet_pt.clear();
      jets.tight_bjet_eta.clear();
      jets.tight_bjet_phi.clear();
      jets.tight_bjet_size.clear();
      jets.v_all_loose_bjets_pt.clear();
      jets.v_all_loose_bjets_eta.clear();

      drLooseBjet4Mu.DR_LooseBJet_LooseMu1.clear();
      drLooseBjet4Mu.DR_LooseBJet_LooseMu2.clear();
      drLooseBjet4Mu.DR_LooseBJet_LooseMu3.clear();
      drLooseBjet4Mu.DR_LooseBJet_LooseMu4.clear();
      drLooseBjet4Mu.DR_LooseBJet_TightMu1.clear();
      drLooseBjet4Mu.DR_LooseBJet_TightMu2.clear();
      drLooseBjet4Mu.DR_LooseBJet_TightMu3.clear();
      drLooseBjet4Mu.DR_LooseBJet_TightMu4.clear();

      drMedBjet4Mu.DR_MedBJet_LooseMu1.clear();
      drMedBjet4Mu.DR_MedBJet_LooseMu2.clear();
      drMedBjet4Mu.DR_MedBJet_LooseMu3.clear();
      drMedBjet4Mu.DR_MedBJet_LooseMu4.clear();
      drMedBjet4Mu.DR_MedBJet_TightMu1.clear();
      drMedBjet4Mu.DR_MedBJet_TightMu2.clear();
      drMedBjet4Mu.DR_MedBJet_TightMu3.clear();
      drMedBjet4Mu.DR_MedBJet_TightMu4.clear();

      drTightBjet4Mu.DR_TightBJet_LooseMu1.clear();
      drTightBjet4Mu.DR_TightBJet_LooseMu2.clear();
      drTightBjet4Mu.DR_TightBJet_LooseMu3.clear();
      drTightBjet4Mu.DR_TightBJet_LooseMu4.clear();
      drTightBjet4Mu.DR_TightBJet_TightMu1.clear();
      drTightBjet4Mu.DR_TightBJet_TightMu2.clear();
      drTightBjet4Mu.DR_TightBJet_TightMu3.clear();
      drTightBjet4Mu.DR_TightBJet_TightMu4.clear();

      met.met_MET.clear();
      met.met_eta.clear();
      met.met_phi.clear();
      
      four_muons_beforeCut.muon1_pt.clear();
      four_muons_beforeCut.muon2_pt.clear();
      four_muons_beforeCut.muon3_pt.clear();
      four_muons_beforeCut.muon4_pt.clear();
      four_muons_beforeCut.muon1_eta.clear();
      four_muons_beforeCut.muon2_eta.clear();
      four_muons_beforeCut.muon3_eta.clear();
      four_muons_beforeCut.muon4_eta.clear();
      four_muons_beforeCut.muon1_phi.clear();
      four_muons_beforeCut.muon2_phi.clear();
      four_muons_beforeCut.muon3_phi.clear();
      four_muons_beforeCut.muon4_phi.clear();
      drMuons.v_DR_mu1mu2.clear();
      drMuons.v_DR_mu3mu4.clear();
      drMuons.v_DR_mu1mu3.clear();
      drMuons.v_DR_mu2mu4.clear();
      drMuons.v_DR_mu1mu4.clear();
      drMuons.v_DR_mu2mu3.clear();  
      zmass.Z12_mass.clear();
      zmass.Z34_mass.clear();
      zmass.Z13_mass.clear();
      zmass.Z24_mass.clear();
      zmass.Z14_mass.clear();
      zmass.Z23_mass.clear();
      dR_bjets_mu_NoCUT.v_DR_bmu1.clear();
      dR_bjets_mu_NoCUT.v_DR_bmu2.clear();
      dR_bjets_mu_NoCUT.v_DR_bmu3.clear();
      dR_bjets_mu_NoCUT.v_DR_bmu4.clear();
      four_muons_afterFullSel.muon1_pt.clear();
      four_muons_afterFullSel.muon2_pt.clear();
      four_muons_afterFullSel.muon3_pt.clear();
      four_muons_afterFullSel.muon4_pt.clear();
      four_muons_afterFullSel.muon1_eta.clear();
      four_muons_afterFullSel.muon2_eta.clear();
      four_muons_afterFullSel.muon3_eta.clear();
      four_muons_afterFullSel.muon4_eta.clear();
      four_muons_afterFullSel.muon1_phi.clear();
      four_muons_afterFullSel.muon2_phi.clear();
      four_muons_afterFullSel.muon3_phi.clear();
      four_muons_afterFullSel.muon4_phi.clear();
      drMuons_FullSel.v_DR_mu1mu2.clear();
      drMuons_FullSel.v_DR_mu3mu4.clear();
      drMuons_FullSel.v_DR_mu1mu3.clear();
      drMuons_FullSel.v_DR_mu2mu4.clear();
      drMuons_FullSel.v_DR_mu1mu4.clear();
      drMuons_FullSel.v_DR_mu2mu3.clear();*/

      za_sig.Za_mass.clear();
      za_sig.Za_pt.clear();
      za_sig.Za_eta.clear(); 

      zb_sig.Zb_mass.clear();
      zb_sig.Zb_pt.clear();
      zb_sig.Zb_eta.clear();

      b1_sig.b1_mass.clear();
      b1_sig.b1_pt.clear();
      b1_sig.b1_eta.clear();
      b1_sig.DR_b1_mu1.clear();
      b1_sig.DR_b1_mu2.clear();
      b1_sig.DR_b1_mu3.clear();
      b1_sig.DR_b1_mu4.clear();

      b2_sig.b2_mass.clear();
      b2_sig.b2_pt.clear();
      b2_sig.b2_eta.clear();
      b2_sig.DR_b2_mu1.clear();
      b2_sig.DR_b2_mu2.clear();
      b2_sig.DR_b2_mu3.clear();
      b2_sig.DR_b2_mu4.clear();

      smhiggs1.v_h1_mass.clear();
      smhiggs1.v_h1_pt.clear();
      smhiggs1.v_h1_eta.clear();

      smhiggs2.v_h2_mass.clear();
      smhiggs2.v_h2_pt.clear();
      smhiggs2.v_h2_eta.clear();

      heavyHiggs.v_H_mass.clear();
      heavyHiggs.v_H_pt.clear();
      heavyHiggs.v_H_eta.clear();

      dr_b1b2.clear();
      dr_h1h2.clear(); 
            
      
       //------------------------GEN PARTICLES------------------------------ 	 
		 
	   // Get Thresholds for pT, eta for generated objects 
	 
	  // Looping overall Gen_Particles 	 
      /* for (Int_t i = 0; i < Particle_size; i++){
		 
	        // Get particle pdg id 
	        int p_id = Particle_PID[i];
		 
	        // Check for Muons
	        if ( p_id == 13 ) { // pdg_id = 13 for Muon

		  Float_t gen_mu_pt = Particle_PT[i];
	          Float_t gen_mu_eta = Particle_Eta[i];
		  Float_t gen_mu_phi = Particle_Phi[i];
		       
		  genMuons.gen_muon_pt.push_back(gen_mu_pt);
		  genMuons.gen_muon_eta.push_back(gen_mu_eta);
		  genMuons.gen_muon_phi.push_back(gen_mu_phi);
		       
		      
	        } // end if on id 13 	
	     
	        // Check for b quarks
	        if ( p_id == 5 ) { // pdg_id = 5 for b quark 
			 
		  Float_t gen_b_pt = Particle_PT[i];
		  Float_t gen_b_eta = Particle_Eta[i];
		  Float_t gen_b_phi = Particle_Phi[i]; 
			
		  gen_bjet.gen_bjet_pt.push_back(gen_b_pt);
		  gen_bjet.gen_bjet_eta.push_back(gen_b_eta);
		  gen_bjet.gen_bjet_phi.push_back(gen_b_phi);
			 
			
	         } // end if on id 5
		 
           } // end loop over gen particles    
           
      */
      
      //------------Loop over Loose Muons--------------
      /*    for (Int_t i = 0; i < MuonLoose_size; i++){
	
	Float_t muon_pt = MuonLoose_PT[i];
	Float_t muon_eta = MuonLoose_Eta[i];
	Float_t muon_phi = MuonLoose_Phi[i];
			
	loose_mu.lmu_pt.push_back(muon_pt);
	loose_mu.lmu_eta.push_back(muon_eta);
	loose_mu.lmu_phi.push_back(muon_phi);
	
		
       }

       // Number of Loose Muons / event
       Int_t nb_loose_muons = MuonLoose_size;
       loose_mu.lmu_size.push_back(nb_loose_muons);

       // kinematics of 4 muons no cuts 
       loose_mu.muon1_pt.push_back(  MuonLoose_PT[0]  );
       loose_mu.muon2_pt.push_back(  MuonLoose_PT[1]  );
       loose_mu.muon3_pt.push_back(  MuonLoose_PT[2]  );
       loose_mu.muon4_pt.push_back(  MuonLoose_PT[3]  );
       loose_mu.muon1_eta.push_back( MuonLoose_Eta[0] );
       loose_mu.muon2_eta.push_back( MuonLoose_Eta[1] );
       loose_mu.muon3_eta.push_back( MuonLoose_Eta[2] );
       loose_mu.muon4_eta.push_back( MuonLoose_Eta[3] );
       loose_mu.muon1_phi.push_back( MuonLoose_Phi[0] );
       loose_mu.muon2_phi.push_back( MuonLoose_Phi[1] );
       loose_mu.muon3_phi.push_back( MuonLoose_Phi[2] );
       loose_mu.muon4_phi.push_back( MuonLoose_Phi[3] );

       // impact parameter (distance from interaction point "Primary Vertex" determined by dxy max value) of 4 Muons no cuts,  later in cuts dxy_muon < dxy max value for this muon
       loose_mu.dxy_mu1.push_back( MuonLoose_D0[0] );
       loose_mu.dxy_mu2.push_back( MuonLoose_D0[1] );
       loose_mu.dxy_mu3.push_back( MuonLoose_D0[2] );
       loose_mu.dxy_mu4.push_back( MuonLoose_D0[3] );

       
       // Delta R between 4 Muons
       Float_t DR_Loose_mu12, DR_Loose_mu34, DR_Loose_mu13, DR_Loose_mu24, DR_Loose_mu14, DR_Loose_mu23; 
      
       DR_Loose_mu12 = TMath::Sqrt(TMath::Power((MuonLoose_Eta[0] - MuonLoose_Eta[1]), 2) + TMath::Power((MuonLoose_Phi[0] - MuonLoose_Phi[1]), 2));
       DR_Loose_mu34 = TMath::Sqrt(TMath::Power((MuonLoose_Eta[2] - MuonLoose_Eta[3]), 2) + TMath::Power((MuonLoose_Phi[2] - MuonLoose_Phi[3]), 2));
       DR_Loose_mu13 = TMath::Sqrt(TMath::Power((MuonLoose_Eta[0] - MuonLoose_Eta[2]), 2) + TMath::Power((MuonLoose_Phi[0] - MuonLoose_Phi[2]), 2));
       DR_Loose_mu24 = TMath::Sqrt(TMath::Power((MuonLoose_Eta[1] - MuonLoose_Eta[3]), 2) + TMath::Power((MuonLoose_Phi[1] - MuonLoose_Phi[3]), 2));
       DR_Loose_mu14 = TMath::Sqrt(TMath::Power((MuonLoose_Eta[0] - MuonLoose_Eta[3]), 2) + TMath::Power((MuonLoose_Phi[0] - MuonLoose_Phi[3]), 2));
       DR_Loose_mu23 = TMath::Sqrt(TMath::Power((MuonLoose_Eta[1] - MuonLoose_Eta[2]), 2) + TMath::Power((MuonLoose_Phi[1] - MuonLoose_Phi[2]), 2));

       drMuonsLoose.v_DR_mu1mu2.push_back(DR_Loose_mu12);
       drMuonsLoose.v_DR_mu3mu4.push_back(DR_Loose_mu34);
       drMuonsLoose.v_DR_mu1mu3.push_back(DR_Loose_mu13);
       drMuonsLoose.v_DR_mu2mu4.push_back(DR_Loose_mu24);
       drMuonsLoose.v_DR_mu1mu4.push_back(DR_Loose_mu14);
       drMuonsLoose.v_DR_mu2mu3.push_back(DR_Loose_mu23);


       // mass of opposite charge muon pairs
       TLorentzVector lo_mu1, lo_mu2,lo_mu3, lo_mu4;
       double muon_mass = 0.105658375;  // in GeV

       lo_mu1.SetPtEtaPhiM(MuonLoose_PT[0], MuonLoose_Eta[0], MuonLoose_Phi[0], muon_mass); 
       lo_mu2.SetPtEtaPhiM(MuonLoose_PT[1], MuonLoose_Eta[1], MuonLoose_Phi[1], muon_mass);
       lo_mu3.SetPtEtaPhiM(MuonLoose_PT[2], MuonLoose_Eta[2], MuonLoose_Phi[2], muon_mass);
       lo_mu4.SetPtEtaPhiM(MuonLoose_PT[3], MuonLoose_Eta[3], MuonLoose_Phi[3], muon_mass);  
      
       double m12_l = (lo_mu1 + lo_mu2).M();
       double m34_l = (lo_mu3 + lo_mu4).M();
       double m13_l = (lo_mu1 + lo_mu3).M();
       double m24_l = (lo_mu2 + lo_mu4).M();
       double m14_l = (lo_mu1 + lo_mu4).M();
       double m23_l = (lo_mu2 + lo_mu3).M();

       loose_mu.loose_mu12_mass.push_back(m12_l);
       loose_mu.loose_mu34_mass.push_back(m34_l);
       loose_mu.loose_mu13_mass.push_back(m13_l);
       loose_mu.loose_mu24_mass.push_back(m24_l);
       loose_mu.loose_mu14_mass.push_back(m14_l);
       loose_mu.loose_mu23_mass.push_back(m23_l);




       
      //------------Loop over Tight Muons--------------
      for (Int_t i = 0; i < MuonTight_size; i++){
	
	Float_t muon_pt = MuonTight_PT[i];
	Float_t muon_eta = MuonTight_Eta[i];
	Float_t muon_phi = MuonTight_Phi[i];
			
	tight_mu.tmu_pt.push_back(muon_pt);
	tight_mu.tmu_eta.push_back(muon_eta);
	tight_mu.tmu_phi.push_back(muon_phi);
			
		
       }

      // Number of Tight Muons / event
      Int_t nb_tight_muons = MuonTight_size;
      tight_mu.tmu_size.push_back(nb_tight_muons);

      tight_mu.muon1_pt.push_back(  MuonTight_PT[0]  );
      tight_mu.muon2_pt.push_back(  MuonTight_PT[1]  );
      tight_mu.muon3_pt.push_back(  MuonTight_PT[2]  );
      tight_mu.muon4_pt.push_back(  MuonTight_PT[3]  );
      tight_mu.muon1_eta.push_back( MuonTight_Eta[0] );
      tight_mu.muon2_eta.push_back( MuonTight_Eta[1] );
      tight_mu.muon3_eta.push_back( MuonTight_Eta[2] );
      tight_mu.muon4_eta.push_back( MuonTight_Eta[3] );
      tight_mu.muon1_phi.push_back( MuonTight_Phi[0] );
      tight_mu.muon2_phi.push_back( MuonTight_Phi[1] );
      tight_mu.muon3_phi.push_back( MuonTight_Phi[2] );
      tight_mu.muon4_phi.push_back( MuonTight_Phi[3] );

      // impact parameter in xy plane (distance from interaction point "Primary Vertex" determined by dxy max value) of 4 Muons no cuts,  later in cuts dxy_muon < dxy max value for this muon
      tight_mu.dxy_mu1.push_back( abs(MuonTight_D0[0]) );
      tight_mu.dxy_mu2.push_back( abs(MuonTight_D0[1]) );
      tight_mu.dxy_mu3.push_back( abs(MuonTight_D0[2]) );
      tight_mu.dxy_mu4.push_back( abs(MuonTight_D0[3]) );

      // impact parameter component in z-direction along beam line 
      tight_mu.dz_mu1.push_back( abs(MuonTight_DZ[0]) );
      tight_mu.dz_mu2.push_back( abs(MuonTight_DZ[1]) );
      tight_mu.dz_mu3.push_back( abs(MuonTight_DZ[2]) );
      tight_mu.dz_mu4.push_back( abs(MuonTight_DZ[3]) );
 
      
      // Delta R between 4 Muons
       Float_t DR_Tight_mu12, DR_Tight_mu34, DR_Tight_mu13, DR_Tight_mu24, DR_Tight_mu14, DR_Tight_mu23; 
      
       DR_Tight_mu12 = TMath::Sqrt(TMath::Power((MuonTight_Eta[0] - MuonTight_Eta[1]), 2) + TMath::Power((MuonTight_Phi[0] - MuonTight_Phi[1]), 2));
       DR_Tight_mu34 = TMath::Sqrt(TMath::Power((MuonTight_Eta[2] - MuonTight_Eta[3]), 2) + TMath::Power((MuonTight_Phi[2] - MuonTight_Phi[3]), 2));
       DR_Tight_mu13 = TMath::Sqrt(TMath::Power((MuonTight_Eta[0] - MuonTight_Eta[2]), 2) + TMath::Power((MuonTight_Phi[0] - MuonTight_Phi[2]), 2));
       DR_Tight_mu24 = TMath::Sqrt(TMath::Power((MuonTight_Eta[1] - MuonTight_Eta[3]), 2) + TMath::Power((MuonTight_Phi[1] - MuonTight_Phi[3]), 2));
       DR_Tight_mu14 = TMath::Sqrt(TMath::Power((MuonTight_Eta[0] - MuonTight_Eta[3]), 2) + TMath::Power((MuonTight_Phi[0] - MuonTight_Phi[3]), 2));
       DR_Tight_mu23 = TMath::Sqrt(TMath::Power((MuonTight_Eta[1] - MuonTight_Eta[2]), 2) + TMath::Power((MuonTight_Phi[1] - MuonTight_Phi[2]), 2));

       drMuonsTight.v_DR_mu1mu2.push_back(DR_Tight_mu12);
       drMuonsTight.v_DR_mu3mu4.push_back(DR_Tight_mu34);
       drMuonsTight.v_DR_mu1mu3.push_back(DR_Tight_mu13);
       drMuonsTight.v_DR_mu2mu4.push_back(DR_Tight_mu24);
       drMuonsTight.v_DR_mu1mu4.push_back(DR_Tight_mu14);
       drMuonsTight.v_DR_mu2mu3.push_back(DR_Tight_mu23);

       // mass of opposite charge muon pairs
       TLorentzVector ti_mu1, ti_mu2,ti_mu3, ti_mu4;
      

       ti_mu1.SetPtEtaPhiM(MuonTight_PT[0], MuonTight_Eta[0], MuonTight_Phi[0], muon_mass); 
       ti_mu2.SetPtEtaPhiM(MuonTight_PT[1], MuonTight_Eta[1], MuonTight_Phi[1], muon_mass);
       ti_mu3.SetPtEtaPhiM(MuonTight_PT[2], MuonTight_Eta[2], MuonTight_Phi[2], muon_mass);
       ti_mu4.SetPtEtaPhiM(MuonTight_PT[3], MuonTight_Eta[3], MuonTight_Phi[3], muon_mass);  
      
       double m12_t = (ti_mu1 + ti_mu2).M();
       double m34_t = (ti_mu3 + ti_mu4).M();
       double m13_t = (ti_mu1 + ti_mu3).M();
       double m24_t = (ti_mu2 + ti_mu4).M();
       double m14_t = (ti_mu1 + ti_mu4).M();
       double m23_t = (ti_mu2 + ti_mu3).M();
   
       tight_mu.tight_mu12_mass.push_back(m12_t);
       tight_mu.tight_mu34_mass.push_back(m34_t);
       tight_mu.tight_mu13_mass.push_back(m13_t);
       tight_mu.tight_mu24_mass.push_back(m24_t);
       tight_mu.tight_mu14_mass.push_back(m14_t);
       tight_mu.tight_mu23_mass.push_back(m23_t); 



       
      //----------------Loop over Jets-----------------

      // counters for bjets
      Int_t nb_loose_bjets  = 0;
      Int_t nb_med_bjets = 0;
      Int_t nb_tight_bjets  = 0;

      for ( Int_t i = 0; i < Jet_size; i++ ){
		   
	 Float_t j_pt = Jet_PT[i];
	 Float_t j_eta = Jet_Eta[i];
	 Float_t j_phi = Jet_Phi[i];
	 UInt_t jet_bTag = Jet_BTag[i];
		   
	 jets.jet_pt.push_back(j_pt);
	 jets.jet_eta.push_back(j_eta);
	 jets.jet_phi.push_back(j_phi);


	 // loose bjets
	 Bool_t isLoose = ( jet_bTag & (1 << 0) );

	 if ( isLoose == 1 ){

            jets.loose_bjet_pt.push_back(j_pt);
            jets.loose_bjet_eta.push_back(j_eta);
            jets.loose_bjet_phi.push_back(j_phi);
            nb_loose_bjets++;

	    // calc DR bjet & each of 4 muons

	    double DR_loosebjet_looseMu1 = TMath::Sqrt(TMath::Power((Jet_Eta[i] - MuonLoose_Eta[0]), 2) + TMath::Power((Jet_Phi[i] - MuonLoose_Phi[0]), 2));
            double DR_loosebjet_looseMu2 = TMath::Sqrt(TMath::Power((Jet_Eta[i] - MuonLoose_Eta[1]), 2) + TMath::Power((Jet_Phi[i] - MuonLoose_Phi[1]), 2));
	    double DR_loosebjet_looseMu3 = TMath::Sqrt(TMath::Power((Jet_Eta[i] - MuonLoose_Eta[2]), 2) + TMath::Power((Jet_Phi[i] - MuonLoose_Phi[2]), 2));
            double DR_loosebjet_looseMu4 = TMath::Sqrt(TMath::Power((Jet_Eta[i] - MuonLoose_Eta[3]), 2) + TMath::Power((Jet_Phi[i] - MuonLoose_Phi[3]), 2));

	    double DR_loosebjet_tightMu1 = TMath::Sqrt(TMath::Power((Jet_Eta[i] - MuonTight_Eta[0]), 2) + TMath::Power((Jet_Phi[i] - MuonTight_Phi[0]), 2));
            double DR_loosebjet_tightMu2 = TMath::Sqrt(TMath::Power((Jet_Eta[i] - MuonTight_Eta[1]), 2) + TMath::Power((Jet_Phi[i] - MuonTight_Phi[1]), 2));
	    double DR_loosebjet_tightMu3 = TMath::Sqrt(TMath::Power((Jet_Eta[i] - MuonTight_Eta[2]), 2) + TMath::Power((Jet_Phi[i] - MuonTight_Phi[2]), 2));
            double DR_loosebjet_tightMu4 = TMath::Sqrt(TMath::Power((Jet_Eta[i] - MuonTight_Eta[3]), 2) + TMath::Power((Jet_Phi[i] - MuonTight_Phi[3]), 2));

	    drLooseBjet4Mu.DR_LooseBJet_LooseMu1.push_back(DR_loosebjet_looseMu1);
            drLooseBjet4Mu.DR_LooseBJet_LooseMu2.push_back(DR_loosebjet_looseMu2);
            drLooseBjet4Mu.DR_LooseBJet_LooseMu3.push_back(DR_loosebjet_looseMu3);
            drLooseBjet4Mu.DR_LooseBJet_LooseMu4.push_back(DR_loosebjet_looseMu4);
            drLooseBjet4Mu.DR_LooseBJet_TightMu1.push_back(DR_loosebjet_tightMu1);
            drLooseBjet4Mu.DR_LooseBJet_TightMu2.push_back(DR_loosebjet_tightMu2);
            drLooseBjet4Mu.DR_LooseBJet_TightMu3.push_back(DR_loosebjet_tightMu3);
            drLooseBjet4Mu.DR_LooseBJet_TightMu4.push_back(DR_loosebjet_tightMu4);
	    
	 } 


	 // medium bjets
	 Bool_t isMedium = ( jet_bTag & (1 << 1) );

	 if ( isMedium == 1 ){

            jets.med_bjet_pt.push_back(j_pt);
            jets.med_bjet_eta.push_back(j_eta);
            jets.med_bjet_phi.push_back(j_phi);
            nb_med_bjets++;

	    // calc DR bjet & each of 4 muons

	    double DR_medbjet_looseMu1 = TMath::Sqrt(TMath::Power((Jet_Eta[i] - MuonLoose_Eta[0]), 2) + TMath::Power((Jet_Phi[i] - MuonLoose_Phi[0]), 2));
            double DR_medbjet_looseMu2 = TMath::Sqrt(TMath::Power((Jet_Eta[i] - MuonLoose_Eta[1]), 2) + TMath::Power((Jet_Phi[i] - MuonLoose_Phi[1]), 2));
	    double DR_medbjet_looseMu3 = TMath::Sqrt(TMath::Power((Jet_Eta[i] - MuonLoose_Eta[2]), 2) + TMath::Power((Jet_Phi[i] - MuonLoose_Phi[2]), 2));
            double DR_medbjet_looseMu4 = TMath::Sqrt(TMath::Power((Jet_Eta[i] - MuonLoose_Eta[3]), 2) + TMath::Power((Jet_Phi[i] - MuonLoose_Phi[3]), 2));

	    double DR_medbjet_tightMu1 = TMath::Sqrt(TMath::Power((Jet_Eta[i] - MuonTight_Eta[0]), 2) + TMath::Power((Jet_Phi[i] - MuonTight_Phi[0]), 2));
            double DR_medbjet_tightMu2 = TMath::Sqrt(TMath::Power((Jet_Eta[i] - MuonTight_Eta[1]), 2) + TMath::Power((Jet_Phi[i] - MuonTight_Phi[1]), 2));
	    double DR_medbjet_tightMu3 = TMath::Sqrt(TMath::Power((Jet_Eta[i] - MuonTight_Eta[2]), 2) + TMath::Power((Jet_Phi[i] - MuonTight_Phi[2]), 2));
            double DR_medbjet_tightMu4 = TMath::Sqrt(TMath::Power((Jet_Eta[i] - MuonTight_Eta[3]), 2) + TMath::Power((Jet_Phi[i] - MuonTight_Phi[3]), 2));

	    drMedBjet4Mu.DR_MedBJet_LooseMu1.push_back(DR_medbjet_looseMu1);
            drMedBjet4Mu.DR_MedBJet_LooseMu2.push_back(DR_medbjet_looseMu2);
            drMedBjet4Mu.DR_MedBJet_LooseMu3.push_back(DR_medbjet_looseMu3);
            drMedBjet4Mu.DR_MedBJet_LooseMu4.push_back(DR_medbjet_looseMu4);
            drMedBjet4Mu.DR_MedBJet_TightMu1.push_back(DR_medbjet_tightMu1);
            drMedBjet4Mu.DR_MedBJet_TightMu2.push_back(DR_medbjet_tightMu2);
            drMedBjet4Mu.DR_MedBJet_TightMu3.push_back(DR_medbjet_tightMu3);
            drMedBjet4Mu.DR_MedBJet_TightMu4.push_back(DR_medbjet_tightMu4);
	    

	 }


	 // tight bjets
	 Bool_t isTight = ( jet_bTag & (1 << 2) );

	 if ( isTight == 1 ){

            jets.tight_bjet_pt.push_back(j_pt);
            jets.tight_bjet_eta.push_back(j_eta);
            jets.tight_bjet_phi.push_back(j_phi);
            nb_tight_bjets++;

	    double DR_tightbjet_looseMu1 = TMath::Sqrt(TMath::Power((Jet_Eta[i] - MuonLoose_Eta[0]), 2) + TMath::Power((Jet_Phi[i] - MuonLoose_Phi[0]), 2));
            double DR_tightbjet_looseMu2 = TMath::Sqrt(TMath::Power((Jet_Eta[i] - MuonLoose_Eta[1]), 2) + TMath::Power((Jet_Phi[i] - MuonLoose_Phi[1]), 2));
	    double DR_tightbjet_looseMu3 = TMath::Sqrt(TMath::Power((Jet_Eta[i] - MuonLoose_Eta[2]), 2) + TMath::Power((Jet_Phi[i] - MuonLoose_Phi[2]), 2));
            double DR_tightbjet_looseMu4 = TMath::Sqrt(TMath::Power((Jet_Eta[i] - MuonLoose_Eta[3]), 2) + TMath::Power((Jet_Phi[i] - MuonLoose_Phi[3]), 2));

	    double DR_tightbjet_tightMu1 = TMath::Sqrt(TMath::Power((Jet_Eta[i] - MuonTight_Eta[0]), 2) + TMath::Power((Jet_Phi[i] - MuonTight_Phi[0]), 2));
            double DR_tightbjet_tightMu2 = TMath::Sqrt(TMath::Power((Jet_Eta[i] - MuonTight_Eta[1]), 2) + TMath::Power((Jet_Phi[i] - MuonTight_Phi[1]), 2));
	    double DR_tightbjet_tightMu3 = TMath::Sqrt(TMath::Power((Jet_Eta[i] - MuonTight_Eta[2]), 2) + TMath::Power((Jet_Phi[i] - MuonTight_Phi[2]), 2));
            double DR_tightbjet_tightMu4 = TMath::Sqrt(TMath::Power((Jet_Eta[i] - MuonTight_Eta[3]), 2) + TMath::Power((Jet_Phi[i] - MuonTight_Phi[3]), 2));

	    drTightBjet4Mu.DR_TightBJet_LooseMu1.push_back(DR_tightbjet_looseMu1);
            drTightBjet4Mu.DR_TightBJet_LooseMu2.push_back(DR_tightbjet_looseMu2);
            drTightBjet4Mu.DR_TightBJet_LooseMu3.push_back(DR_tightbjet_looseMu3);
            drTightBjet4Mu.DR_TightBJet_LooseMu4.push_back(DR_tightbjet_looseMu4);
            drTightBjet4Mu.DR_TightBJet_TightMu1.push_back(DR_tightbjet_tightMu1);
            drTightBjet4Mu.DR_TightBJet_TightMu2.push_back(DR_tightbjet_tightMu2);
            drTightBjet4Mu.DR_TightBJet_TightMu3.push_back(DR_tightbjet_tightMu3);
            drTightBjet4Mu.DR_TightBJet_TightMu4.push_back(DR_tightbjet_tightMu4);
	    

	 }
	 
		   
      } // end loop over jets

      
       // Number of Jets / event
       Int_t nb_jets = Jet_size;

       jets.jet_size.push_back(nb_jets);
       jets.loose_bjet_size.push_back(nb_loose_bjets);
       jets.med_bjet_size.push_back(nb_med_bjets);
       jets.tight_bjet_size.push_back(nb_tight_bjets);


       
      //---------------Loop over MET----------------
       
      for (Int_t i = 0; i < MissingET_size; i++){
		  
	Float_t missingET = MissingET_MET[i];
	Float_t MET_eta   = MissingET_Eta[i];
        Float_t MET_phi   = MissingET_Phi[i];
	        
	met.met_MET.push_back(missingET);
	met.met_eta.push_back(MET_eta);
	met.met_phi.push_back(MET_phi);
	       
       }
       
       */ 
          
       //=====================================================================//
       //                        Start 4Muons , ZZ Selections                 //
       //=====================================================================//
	      
	    
       
       //cout << "Original Number of Muons: " << MuonTight_size << endl;
	    
       for ( Int_t i = 0; i < MuonTight_size; i++ ){
	 
	  Float_t muon_pt = MuonTight_PT[i];
	  Float_t muon_eta = MuonTight_Eta[i];
	  Float_t abseta = abs(muon_eta);
	  Float_t muon_phi = MuonTight_Phi[i];
			
 	  tight_mu.tmu_pt.push_back(muon_pt);
	  tight_mu.tmu_eta.push_back(muon_eta);
	  tight_mu.tmu_phi.push_back(muon_phi);
	  
		     
	  // 1st Selection: on pT, eta of reconstructed Muons (Object_Selection) 
			
	  if ( muon_pt > 5. && abseta < 2.7 ){

	   			
	    //NEvents[0]++;
					
	       int mu_idx = i;
	       v_muon_idx.push_back(mu_idx);
	       
           
	   
          } // end if muonL_pt > 5 & abs eta < 2.7
       } // end loop overall loose muons		       

       if ( MuonTight_PT[0] > 5. && MuonTight_PT[1] > 5.
	    && MuonTight_PT[2] > 5. && MuonTight_PT[3] > 5.){

	 NEvents[0]++;

         if ( MuonTight_Eta[0] < 2.7 && MuonTight_Eta[1] < 2.7
	    && MuonTight_Eta[2] < 2.7 && MuonTight_Eta[3] < 2.7 )NEvents[1]++;
       }
      
       // 2nd Selection: on number of Muons per event

      // if ( v_muon_idx.size() > 3 ){ // having at least 4 Muons per event
       if ( v_muon_idx.size() == 4 ) { // having exactly 4 Muons per eevent
			
         //cout << "muon indx vector size :- " << v_muon_idx.size() << endl;

	 NEvents[2]++;
			  

	 // TLorentzVector declarations 
         TLorentzVector mu1, mu2, mu3, mu4, Z1, Z2, b1, b2, h1, h2, H;
      

	 // Set Pt, eta, phi mass for 4 muons TLV

	 float muon_mass = 0.105658375;  // in GeV
         //double m_mu1, m_mu2, m_mu3, m_mu4;
         Float_t pt_mu1, pt_mu2, pt_mu3, pt_mu4;
         Float_t eta_mu1, eta_mu2, eta_mu3, eta_mu4;
         Float_t phi_mu1, phi_mu2, phi_mu3, phi_mu4;
	 Float_t m12, m34, m13, m24, m14, m23, mComb_1, mComb_2, mComb_3;

	 
	 // indicies of 4Muons
         int mu1_idx = v_muon_idx[0];
         int mu2_idx = v_muon_idx[1];
         int mu3_idx = v_muon_idx[2];
         int mu4_idx = v_muon_idx[3];


         // impact parameter in xy plane (distance from interaction point "Primary Vertex"
	 Float_t dxy_mu1 = abs( MuonTight_D0[mu1_idx] );
	 Float_t dxy_mu2 = abs( MuonTight_D0[mu2_idx] );
	 Float_t dxy_mu3 = abs( MuonTight_D0[mu3_idx] );
	 Float_t dxy_mu4 = abs( MuonTight_D0[mu4_idx] );

	 
	 // impact parameter component in z-direction along beam line
	 Float_t dz_mu1 = abs( MuonTight_DZ[mu1_idx] );
	 Float_t dz_mu2 = abs( MuonTight_DZ[mu2_idx] );
	 Float_t dz_mu3 = abs( MuonTight_DZ[mu3_idx] );
	 Float_t dz_mu4 = abs( MuonTight_DZ[mu4_idx] );

	 cout << "muon1_id | muon2_id | muon3_id | muon4_id : " << mu1_idx << " | " << mu2_idx << " | " << mu3_idx << " | " << mu4_idx << endl;
	
      
         // In Delphes tree Muon_PT are sorted from highest to least one
         mu1.SetPtEtaPhiM(MuonTight_PT[mu1_idx], MuonTight_Eta[mu1_idx], MuonTight_Phi[mu1_idx], muon_mass); 
         mu2.SetPtEtaPhiM(MuonTight_PT[mu2_idx], MuonTight_Eta[mu2_idx], MuonTight_Phi[mu2_idx], muon_mass);
         mu3.SetPtEtaPhiM(MuonTight_PT[mu3_idx], MuonTight_Eta[mu3_idx], MuonTight_Phi[mu3_idx], muon_mass);
         mu4.SetPtEtaPhiM(MuonTight_PT[mu4_idx], MuonTight_Eta[mu4_idx], MuonTight_Phi[mu4_idx], muon_mass);  
      
         // Leading Muon pT > 35 GeV (Muon with highest pT)
         pt_mu1 = mu1.Pt();
         eta_mu1 = mu1.Eta();
         phi_mu1 = mu1.Phi();
                
         // subleading Muon pT > 20 GeV (Muon with second-highest pT)    
         pt_mu2 = mu2.Pt();
         eta_mu2 = mu2.Eta();
         phi_mu2 = mu2.Phi();
                
	 // next subleading Muon 
	 pt_mu3 = mu3.Pt();
         eta_mu3 = mu3.Eta();
         phi_mu3 = mu3.Phi();
                
	 // next to next subleading Muon
	 pt_mu4 = mu4.Pt();
         eta_mu4 = mu4.Eta();
         phi_mu4 = mu4.Phi();



	 //cout << "mu1_pt | mu2_pt | mu3_pt | mu4_pt : " << pt_mu1 << " | " << pt_mu2 << " | " << pt_mu3 << " | " << pt_mu4  << endl;
	 //cout << "mu1_eta | mu2_eta | mu3_eta | mu4_eta : " << eta_mu1 << " | " << eta_mu2 << " | " << eta_mu3 << " | " << eta_mu4  << endl;
	 
	 
              
         // Delta R between 4 Muons
         Float_t DR_Tight_mu12, DR_Tight_mu34, DR_Tight_mu13, DR_Tight_mu24, DR_Tight_mu14, DR_Tight_mu23; 
      
         DR_Tight_mu12 = TMath::Sqrt(TMath::Power((MuonTight_Eta[mu1_idx] - MuonTight_Eta[mu2_idx]), 2) + TMath::Power((MuonTight_Phi[mu1_idx] - MuonTight_Phi[mu2_idx]), 2));
         DR_Tight_mu34 = TMath::Sqrt(TMath::Power((MuonTight_Eta[mu3_idx] - MuonTight_Eta[mu4_idx]), 2) + TMath::Power((MuonTight_Phi[mu3_idx] - MuonTight_Phi[mu4_idx]), 2));
         DR_Tight_mu13 = TMath::Sqrt(TMath::Power((MuonTight_Eta[mu1_idx] - MuonTight_Eta[mu3_idx]), 2) + TMath::Power((MuonTight_Phi[mu1_idx] - MuonTight_Phi[mu3_idx]), 2));
         DR_Tight_mu24 = TMath::Sqrt(TMath::Power((MuonTight_Eta[mu2_idx] - MuonTight_Eta[mu4_idx]), 2) + TMath::Power((MuonTight_Phi[mu2_idx] - MuonTight_Phi[mu4_idx]), 2));
         DR_Tight_mu14 = TMath::Sqrt(TMath::Power((MuonTight_Eta[mu1_idx] - MuonTight_Eta[mu4_idx]), 2) + TMath::Power((MuonTight_Phi[mu1_idx] - MuonTight_Phi[mu4_idx]), 2));
         DR_Tight_mu23 = TMath::Sqrt(TMath::Power((MuonTight_Eta[mu2_idx] - MuonTight_Eta[mu3_idx]), 2) + TMath::Power((MuonTight_Phi[mu2_idx] - MuonTight_Phi[mu3_idx]), 2));

         
         
	 // invariant mass for each muon pair
	 //	 double m12 = (mu1 + mu2).M();
	 // double m34 = (mu3 + mu4).M();
	 /* double m13 = (mu1 + mu3).M();
	 double m24 = (mu2 + mu4).M();
         double m14 = (mu1 + mu4).M();
	 double m23 = (mu2 + mu3).M(); */

	 
         
         // Nominal Z mass 
         float mZ = 91.1876; // in GeV

	 Float_t  mZ12, mZ34, mZ13, mZ24, mZ14, mZ23, mZ1, mZ2, ptZ1, ptZ2, etaZ1, etaZ2, phiZ1, phiZ2;
      
         
         // 3rd Selection: on 4Muons Total charge

	 // Selection on IP D0,DZ
	 if ( (dxy_mu1 < 0.12) && (dxy_mu2 < 0.12)
	       && (dxy_mu3 < 0.12) && (dxy_mu4 < 0.12) ){

            if ( (dz_mu1 < 0.47) && (dz_mu2 < 0.47)
		 && (dz_mu3 < 0.47) && (dz_mu4 < 0.47) ){

	        NEvents[3]++;

	
	       // Selection on Total Charge of 4 Muons should be = 0 
               if ( MuonTight_Charge[mu1_idx] + MuonTight_Charge[mu2_idx] + MuonTight_Charge[mu3_idx] + MuonTight_Charge[mu4_idx] == 0){ //Sure that muon pairs are of opposite signs 
       
                  NEvents[4]++;
		 
                 // Selection on Leading & subleading mu pT
		 if ( ( pt_mu1 > 20. ) && ( pt_mu2 > 10. ) ) {

		    NEvents[5]++;

                    cout << " is Case 1 "     << isCase1<< endl;
	            cout << " is Case 2 "     << isCase2<< endl;

	            if ( ( MuonTight_Charge[mu1_idx] + MuonTight_Charge[mu2_idx] == 0) &&    // mu 12&34 opp
				   
		        ( MuonTight_Charge[mu3_idx] + MuonTight_Charge[mu4_idx] == 0) ) {

	                NEvents[6]++;
	                isCase1 = true;


	            }

	            if ( ( MuonTight_Charge[mu1_idx] + MuonTight_Charge[mu2_idx] != 0) &&    // mu 12&34 same
				   
		       ( MuonTight_Charge[mu3_idx] + MuonTight_Charge[mu4_idx] != 0) ) {

	               NEvents[7]++;
	               isCase2 = true;


	            }


	    
	      //===========================================================================//
              //                                CASE1 Fullfilled                           //
	      //===========================================================================//
	      
	     cout << " is Category 1 " << catgr1 << endl;
             cout << " is Category 2 " << catgr2 << endl;
	     cout << " is Category 3 " << catgr3 << endl;
	     
	     if ( isCase1 ){

                //------------SAME CHARGE 13, 24------------
		
		// Working for combinations 1234, 1423

		if ( ( MuonTight_Charge[mu1_idx] + MuonTight_Charge[mu3_idx] != 0) &&  
				   
		   ( MuonTight_Charge[mu2_idx] + MuonTight_Charge[mu4_idx] != 0) ){  // 13, 24 are of same charge

		  //NEvents[4]++;


		  // Selection on IP
		  
		  // if ( (dxy_mu1 < 0.12) && (dxy_mu2 < 0.12)
		  // && (dxy_mu3 < 0.12) && (dxy_mu4 < 0.12) ){

		  //  if ( (dz_mu1 < 0.47) && (dz_mu2 < 0.47)
		  //  && (dz_mu3 < 0.47) && (dz_mu4 < 0.47) ){

		       
		  // NEvents[6]++;
			  
		  			  
                          //6th Selection: on Leading (highest pT) and Subleading (second highest pT) Muons
	       
		  // if ( ( pt_mu1 > 20. ) && ( pt_mu2 > 10. ) ){
			    //if ( ( pt_mu1 > 35. ) && ( pt_mu2 > 20. ) ){ 

		  //   NEvents[7]++;

			      // Define masses for muons

			      m12 = (mu1 + mu2).M();
			      m34 = (mu3 + mu4).M();
                              m14 = (mu1 + mu4).M();
                              m23 = (mu2 + mu3).M();

			      if ( m12 > m34 ) mComb_1 = m12;
			      if ( m12 < m34 ) mComb_1 = m34;

			      if ( m14 > m23 ) mComb_2 = m14;
			      if ( m14 < m23 ) mComb_2 = m23;

			      if ( mComb_1 > mComb_2 ) {

				cout << "< Fullfill Category 1 >  " << "mass_Comb_1  " << mComb_1 << "  mass_Comb_2  " << mComb_2 << endl; 
				catgr1 = true;
				cout << " is Category 1 " << catgr1 << endl;

			      }
			      
			      if ( mComb_1 < mComb_2 ) {

                                cout << "< Fullfill Category 2 >  " << "mass_Comb_1  " << mComb_1 << "  mass_Comb_2  " << mComb_2 << endl; 
				catgr2 = true;
				cout << "is Category 2 " << catgr2 << endl;

			      }

			      // Work on 1st Category:- 1234
			      
			      if ( catgr1 ){

				NEvents[8]++;
				cout << " << Working with Category1, muons 1234 >>" << endl;

				if ( mComb_1 == m12 ){

				  cout << " mass_Comb_1 = mass mu12 = " << mComb_1 << endl;

				  if ( m12 > 3. && m34 > 3.) {    // selection on mass of each mu pair

				    if ( DR_Tight_mu12 > 0.07 && DR_Tight_mu34 > 0.07 ) {     // selection on DR

				      // NEvents[7]++;

				      // start filling leptons vectors 1234

				      // Lepton1
				      lep1.push_back(pt_mu1);
				      lep1.push_back(eta_mu1);
				      lep1.push_back(phi_mu1);
				      lep1.push_back( mu1.M() );

				      // Lepton2
				      lep2.push_back(pt_mu2);
				      lep2.push_back(eta_mu2);
				      lep2.push_back(phi_mu2);
				      lep2.push_back( mu2.M() );

				      // Lepton3
				      lep3.push_back(pt_mu3);
				      lep3.push_back(eta_mu3);
				      lep3.push_back(phi_mu3);
				      lep3.push_back( mu3.M() );

				      // Lepton4
				      lep4.push_back(pt_mu4);
				      lep4.push_back(eta_mu4);
				      lep4.push_back(phi_mu4);
				      lep4.push_back( mu4.M() );

				      // Filling 4 Muons of signal Branches

				      // Muon1
				      muon1.pT.push_back(lep1[0]);
				      muon1.eta.push_back(lep1[1]);
				      muon1.phi.push_back(lep1[2]);
				      muon1.mass.push_back(lep1[3]);

				      // Muon2
				      muon2.pT.push_back(lep2[0]);
				      muon2.eta.push_back(lep2[1]);
				      muon2.phi.push_back(lep2[2]);
				      muon2.mass.push_back(lep2[3]);

				      // Muon3
				      muon3.pT.push_back(lep3[0]);
				      muon3.eta.push_back(lep3[1]);
				      muon3.phi.push_back(lep3[2]);
				      muon3.mass.push_back(lep3[3]);

				      // Muon4
				      muon4.pT.push_back(lep4[0]);
				      muon4.eta.push_back(lep4[1]);
				      muon4.phi.push_back(lep4[2]);
				      muon4.mass.push_back(lep4[3]);

				      
				      // Z mass distribution for each muon pair
				      mZ12 = (mu1 + mu2).M();
				      mZ34 = (mu3 + mu4).M();

				      mz.mZ_12.push_back(mZ12);
				      mz.mZ_34.push_back(mZ34);
				      

				      // On-Shell Z Signal
                                      Z1 = mu1 + mu2;
				      mZ1   = Z1.M();   // invariant mass of muons 1&2
		                      ptZ1  = Z1.Pt();
                                      etaZ1 = Z1.Eta();
		                      phiZ1 = Z1.Phi();


                                      // Off-Shell Z
				      Z2 = mu3 + mu4;
                                      mZ2   = Z2.M();   // invariant mass of muons 1&2
		                      ptZ2  = Z2.Pt();
                                      etaZ2 = Z2.Eta();
		                      phiZ2 = Z2.Phi();

				      // Filling Z1 Z2 MASS distributions for Signal
				      za_sig.Za_mass.push_back(mZ1);
				      zb_sig.Zb_mass.push_back(mZ2);

				    } // end DR

				  } // end mass of each mu pair > 3 GeV

				} // end if mComb_1 is m12

				else{  // 3412

				   cout << " mass_Comb_1 = mass mu34 = " << mComb_1 << endl;

				  if ( m12 > 3. && m34 > 3.) {    // selection on mass of each mu pair

				    if ( DR_Tight_mu12 > 0.07 && DR_Tight_mu34 > 0.07 ) {     // selection on DR

				      //NEvents[8]++;

				       // Lepton1
				       lep1.push_back(pt_mu3);
				       lep1.push_back(eta_mu3);
				       lep1.push_back(phi_mu3);
				       lep1.push_back( mu3.M() );

				       // Lepton2
				       lep2.push_back(pt_mu4);
				       lep2.push_back(eta_mu4);
				       lep2.push_back(phi_mu4);
				       lep2.push_back( mu4.M() );

				       // Lepton3
				       lep3.push_back(pt_mu1);
				       lep3.push_back(eta_mu1);
				       lep3.push_back(phi_mu1);
				       lep3.push_back( mu1.M() );

				       // Lepton4
				       lep4.push_back(pt_mu2);
				       lep4.push_back(eta_mu2);
				       lep4.push_back(phi_mu2);
				       lep4.push_back( mu2.M() );


				       // Filling 4 Muons of signal Branches
				       
				       // Muon1
				       muon1.pT.push_back(lep1[0]);
				       muon1.eta.push_back(lep1[1]);
				       muon1.phi.push_back(lep1[2]);
				       muon1.mass.push_back(lep1[3]);

				       // Muon2
				       muon2.pT.push_back(lep2[0]);
				       muon2.eta.push_back(lep2[1]);
				       muon2.phi.push_back(lep2[2]);
				       muon2.mass.push_back(lep2[3]);

				       // Muon3
				       muon3.pT.push_back(lep3[0]);
				       muon3.eta.push_back(lep3[1]);
				       muon3.phi.push_back(lep3[2]);
				       muon3.mass.push_back(lep3[3]);

				       // Muon4
				       muon4.pT.push_back(lep4[0]);
				       muon4.eta.push_back(lep4[1]);
				       muon4.phi.push_back(lep4[2]);
				       muon4.mass.push_back(lep4[3]);


				       // Z mass distribution for each muon pair
				       mZ12 = (mu1 + mu2).M();
				       mZ34 = (mu3 + mu4).M();

				       mz.mZ_12.push_back(mZ12);
				       mz.mZ_34.push_back(mZ34);
				      

				       // On-Shell Z
                                       Z1 = mu3 + mu4;
				       mZ1   = Z1.M();   // invariant mass of muons 1&2
		                       ptZ1  = Z1.Pt();
                                       etaZ1 = Z1.Eta();
		                       phiZ1 = Z1.Phi();


                                       // Off-Shell Z
				       Z2 = mu1 + mu2;
                                       mZ2   = Z2.M();   // invariant mass of muons 1&2
		                       ptZ2  = Z2.Pt();
                                       etaZ2 = Z2.Eta();
		                       phiZ2 = Z2.Phi();
				       

                                       // Filling Z1 Z2 MASS distributions for Signal
				       za_sig.Za_mass.push_back(mZ1);
				       zb_sig.Zb_mass.push_back(mZ2);
				       

				    } // end DR
				  } // end mass of each mu pair > 3GeV
				} // end mComb_1 is m34
				
			      } // end 1st category


			      if ( catgr2 ){  // 1423

				NEvents[9]++;

				cout << " << Working with Category2, muons 1423 >>" << endl;

				if ( mComb_2 == m14 ){    // 1423

                                  cout << " mass_Comb_2 = mass mu14 = " << mComb_2 << endl;
				  
				  if ( m14 > 3. && m23 > 3.) {    // selection on mass of each mu pair

				    if ( DR_Tight_mu14 > 0.07 && DR_Tight_mu23 > 0.07 ) {     // selection on DR

				      //NEvents[9]++;

				      // start filling leptons vectors 1423

				      // Lepton1
				      lep1.push_back(pt_mu1);
				      lep1.push_back(eta_mu1);
				      lep1.push_back(phi_mu1);
				      lep1.push_back( mu1.M() );

				      // Lepton2
				      lep2.push_back(pt_mu4);
				      lep2.push_back(eta_mu4);
				      lep2.push_back(phi_mu4);
				      lep2.push_back( mu4.M() );

				      // Lepton3
				      lep3.push_back(pt_mu2);
				      lep3.push_back(eta_mu2);
				      lep3.push_back(phi_mu2);
				      lep3.push_back( mu2.M() );

				      // Lepton4
				      lep4.push_back(pt_mu3);
				      lep4.push_back(eta_mu3);
				      lep4.push_back(phi_mu3);
				      lep4.push_back( mu3.M() );

				      // Filling 4 Muons of signal Branches

				      // Muon1
				      muon1.pT.push_back(lep1[0]);
				      muon1.eta.push_back(lep1[1]);
				      muon1.phi.push_back(lep1[2]);
				      muon1.mass.push_back(lep1[3]);

				      // Muon2
				      muon2.pT.push_back(lep2[0]);
				      muon2.eta.push_back(lep2[1]);
				      muon2.phi.push_back(lep2[2]);
				      muon2.mass.push_back(lep2[3]);

				      // Muon3
				      muon3.pT.push_back(lep3[0]);
				      muon3.eta.push_back(lep3[1]);
				      muon3.phi.push_back(lep3[2]);
				      muon3.mass.push_back(lep3[3]);

				      // Muon4
				      muon4.pT.push_back(lep4[0]);
				      muon4.eta.push_back(lep4[1]);
				      muon4.phi.push_back(lep4[2]);
				      muon4.mass.push_back(lep4[3]);

				      
				      // Z mass distribution for each muon pair
				      mZ14 = (mu1 + mu4).M();
				      mZ23 = (mu2 + mu3).M();

				      mz.mZ_14.push_back(mZ14);
				      mz.mZ_23.push_back(mZ23);
				      

				      // On-Shell Z Signal
                                      Z1 = mu1 + mu4;
				      mZ1   = Z1.M();   // invariant mass of muons 1&2
		                      ptZ1  = Z1.Pt();
                                      etaZ1 = Z1.Eta();
		                      phiZ1 = Z1.Phi();


                                      // Off-Shell Z
				      Z2 = mu2 + mu3;
                                      mZ2   = Z2.M();   // invariant mass of muons 1&2
		                      ptZ2  = Z2.Pt();
                                      etaZ2 = Z2.Eta();
		                      phiZ2 = Z2.Phi();

				      // Filling Z1 Z2 MASS distributions for Signal
				      za_sig.Za_mass.push_back(mZ1);
				      zb_sig.Zb_mass.push_back(mZ2);

				    } // end DR

				  } // end mass of each mu pair > 3 GeV

				} // end if mComb_2 is m14

				else{  // 2314

				   cout << " mass_Comb_2 = mass mu23 = " << mComb_2 << endl;

				  if ( m14 > 3. && m23 > 3.) {    // selection on mass of each mu pair

				    if ( DR_Tight_mu14 > 0.07 && DR_Tight_mu23 > 0.07 ) {     // selection on DR

				      //   NEvents[10]++;

				       // Lepton1
				       lep1.push_back(pt_mu2);
				       lep1.push_back(eta_mu2);
				       lep1.push_back(phi_mu2);
				       lep1.push_back( mu2.M() );

				       // Lepton2
				       lep2.push_back(pt_mu3);
				       lep2.push_back(eta_mu3);
				       lep2.push_back(phi_mu3);
				       lep2.push_back( mu3.M() );

				       // Lepton3
				       lep3.push_back(pt_mu1);
				       lep3.push_back(eta_mu1);
				       lep3.push_back(phi_mu1);
				       lep3.push_back( mu1.M() );

				       // Lepton4
				       lep4.push_back(pt_mu4);
				       lep4.push_back(eta_mu4);
				       lep4.push_back(phi_mu4);
				       lep4.push_back( mu4.M() );


				       // Filling 4 Muons of signal Branches
				       
				       // Muon1
				       muon1.pT.push_back(lep1[0]);
				       muon1.eta.push_back(lep1[1]);
				       muon1.phi.push_back(lep1[2]);
				       muon1.mass.push_back(lep1[3]);

				       // Muon2
				       muon2.pT.push_back(lep2[0]);
				       muon2.eta.push_back(lep2[1]);
				       muon2.phi.push_back(lep2[2]);
				       muon2.mass.push_back(lep2[3]);

				       // Muon3
				       muon3.pT.push_back(lep3[0]);
				       muon3.eta.push_back(lep3[1]);
				       muon3.phi.push_back(lep3[2]);
				       muon3.mass.push_back(lep3[3]);

				       // Muon4
				       muon4.pT.push_back(lep4[0]);
				       muon4.eta.push_back(lep4[1]);
				       muon4.phi.push_back(lep4[2]);
				       muon4.mass.push_back(lep4[3]);


				       // Z mass distribution for each muon pair
				       mZ14 = (mu1 + mu4).M();
				       mZ23 = (mu2 + mu3).M();

				       mz.mZ_14.push_back(mZ14);
				       mz.mZ_23.push_back(mZ23);
				      

				       // On-Shell Z
                                       Z1 = mu2 + mu3;
				       mZ1   = Z1.M();   // invariant mass of muons 1&2
		                       ptZ1  = Z1.Pt();
                                       etaZ1 = Z1.Eta();
		                       phiZ1 = Z1.Phi();


                                       // Off-Shell Z
				       Z2 = mu1 + mu4;
                                       mZ2   = Z2.M();   // invariant mass of muons 1&2
		                       ptZ2  = Z2.Pt();
                                       etaZ2 = Z2.Eta();
		                       phiZ2 = Z2.Phi();
				       

                                       // Filling Z1 Z2 MASS distributions for Signal
				       za_sig.Za_mass.push_back(mZ1);
				       zb_sig.Zb_mass.push_back(mZ2);
				       

				    } // end DR
				  } // end mass of each mu pair > 3GeV
				} // end mComb_2 is m23
				
			      } // end 2nd category

			      //  } // end if on muons pT 
					   //  } // end if dz                
			      //  } // end if dxy
	         } // end mu 13,24 same charge


                 //------------SAME CHARGE 14, 23------------
		
		// Working for combinations 1234, 1324

		if ( ( MuonTight_Charge[mu1_idx] + MuonTight_Charge[mu4_idx] != 0) &&  
				   
		   ( MuonTight_Charge[mu2_idx] + MuonTight_Charge[mu3_idx] != 0) ){  // 14, 23 are of same charge

		  //  NEvents[11]++;


		  // Selection on IP
		  
		  //if ( (dxy_mu1 < 0.12) && (dxy_mu2 < 0.12)
		  // && (dxy_mu3 < 0.12) && (dxy_mu4 < 0.12) ){

		  //  if ( (dz_mu1 < 0.47) && (dz_mu2 < 0.47)
		  //   && (dz_mu3 < 0.47) && (dz_mu4 < 0.47) ){

		       
			//	  NEvents[12]++;
			  
		  			  
                          //6th Selection: on Leading (highest pT) and Subleading (second highest pT) Muons
	       
		  //   if ( ( pt_mu1 > 20. ) && ( pt_mu2 > 10. ) ){
			    //if ( ( pt_mu1 > 35. ) && ( pt_mu2 > 20. ) ){

			    //       NEvents[13]++;

			      // Define masses for muons

			      m12 = (mu1 + mu2).M();
			      m34 = (mu3 + mu4).M();
                              m13 = (mu1 + mu3).M();
                              m24 = (mu2 + mu4).M();

			      if ( m12 > m34 ) mComb_1 = m12;
			      if ( m12 < m34 ) mComb_1 = m34;

			      if ( m13 > m24 ) mComb_3 = m13;
			      if ( m13 < m24 ) mComb_3 = m24;

			      if ( mComb_1 > mComb_3 ) {

                                cout << "< Fullfill Category 1 >  " << "mass_Comb_1  " << mComb_1 << "  mass_Comb_3  " << mComb_3 << endl; 
				catgr1 = true;
				cout << " is Category 1 " << catgr1 << endl;

			      }
			      
			      if ( mComb_1 < mComb_3 ) {

                                cout << "< Fullfill Category 3 >  " << "mass_Comb_1  " << mComb_1 << "  mass_Comb_3  " << mComb_3 << endl; 
				catgr3 = true;
				cout << " is Category 3 " << catgr3 << endl;

			      }

			      // Work on 1st Category:- 1234
			      
			      if ( catgr1 ){

				NEvents[8]++;

				cout << " << Working with Category1, muons 1234 >>" << endl;

				if ( mComb_1 == m12 ){

				  cout << " mass_Comb_1 = mass mu12 = " << mComb_1 << endl;

				  if ( m12 > 3. && m34 > 3.) {    // selection on mass of each mu pair

				    if ( DR_Tight_mu12 > 0.07 && DR_Tight_mu34 > 0.07 ) {     // selection on DR

				      //   NEvents[14]++;

				      // start filling leptons vectors 1234

				      // Lepton1
				      lep1.push_back(pt_mu1);
				      lep1.push_back(eta_mu1);
				      lep1.push_back(phi_mu1);
				      lep1.push_back( mu1.M() );

				      // Lepton2
				      lep2.push_back(pt_mu2);
				      lep2.push_back(eta_mu2);
				      lep2.push_back(phi_mu2);
				      lep2.push_back( mu2.M() );

				      // Lepton3
				      lep3.push_back(pt_mu3);
				      lep3.push_back(eta_mu3);
				      lep3.push_back(phi_mu3);
				      lep3.push_back( mu3.M() );

				      // Lepton4
				      lep4.push_back(pt_mu4);
				      lep4.push_back(eta_mu4);
				      lep4.push_back(phi_mu4);
				      lep4.push_back( mu4.M() );

				      // Filling 4 Muons of signal Branches

				      // Muon1
				      muon1.pT.push_back(lep1[0]);
				      muon1.eta.push_back(lep1[1]);
				      muon1.phi.push_back(lep1[2]);
				      muon1.mass.push_back(lep1[3]);

				      // Muon2
				      muon2.pT.push_back(lep2[0]);
				      muon2.eta.push_back(lep2[1]);
				      muon2.phi.push_back(lep2[2]);
				      muon2.mass.push_back(lep2[3]);

				      // Muon3
				      muon3.pT.push_back(lep3[0]);
				      muon3.eta.push_back(lep3[1]);
				      muon3.phi.push_back(lep3[2]);
				      muon3.mass.push_back(lep3[3]);

				      // Muon4
				      muon4.pT.push_back(lep4[0]);
				      muon4.eta.push_back(lep4[1]);
				      muon4.phi.push_back(lep4[2]);
				      muon4.mass.push_back(lep4[3]);

				      
				      // Z mass distribution for each muon pair
				      mZ12 = (mu1 + mu2).M();
				      mZ34 = (mu3 + mu4).M();

				      mz.mZ_12.push_back(mZ12);
				      mz.mZ_34.push_back(mZ34);
				      

				      // On-Shell Z Signal
                                      Z1 = mu1 + mu2;
				      mZ1   = Z1.M();   // invariant mass of muons 1&2
		                      ptZ1  = Z1.Pt();
                                      etaZ1 = Z1.Eta();
		                      phiZ1 = Z1.Phi();


                                      // Off-Shell Z
				      Z2 = mu3 + mu4;
                                      mZ2   = Z2.M();   // invariant mass of muons 1&2
		                      ptZ2  = Z2.Pt();
                                      etaZ2 = Z2.Eta();
		                      phiZ2 = Z2.Phi();

				      // Filling Z1 Z2 MASS distributions for Signal
				      za_sig.Za_mass.push_back(mZ1);
				      zb_sig.Zb_mass.push_back(mZ2);

				    } // end DR

				  } // end mass of each mu pair > 3 GeV

				} // end if mComb_1 is m12

				else{  // 3412

				  cout << " mass_Comb_1 = mass mu34 = " << mComb_1 << endl;

				  if ( m12 > 3. && m34 > 3.) {    // selection on mass of each mu pair

				    if ( DR_Tight_mu12 > 0.07 && DR_Tight_mu34 > 0.07 ) {     // selection on DR

				      // NEvents[15]++;

				       // Lepton1
				       lep1.push_back(pt_mu3);
				       lep1.push_back(eta_mu3);
				       lep1.push_back(phi_mu3);
				       lep1.push_back( mu3.M() );

				       // Lepton2
				       lep2.push_back(pt_mu4);
				       lep2.push_back(eta_mu4);
				       lep2.push_back(phi_mu4);
				       lep2.push_back( mu4.M() );

				       // Lepton3
				       lep3.push_back(pt_mu1);
				       lep3.push_back(eta_mu1);
				       lep3.push_back(phi_mu1);
				       lep3.push_back( mu1.M() );

				       // Lepton4
				       lep4.push_back(pt_mu2);
				       lep4.push_back(eta_mu2);
				       lep4.push_back(phi_mu2);
				       lep4.push_back( mu2.M() );


				       // Filling 4 Muons of signal Branches
				       
				       // Muon1
				       muon1.pT.push_back(lep1[0]);
				       muon1.eta.push_back(lep1[1]);
				       muon1.phi.push_back(lep1[2]);
				       muon1.mass.push_back(lep1[3]);

				       // Muon2
				       muon2.pT.push_back(lep2[0]);
				       muon2.eta.push_back(lep2[1]);
				       muon2.phi.push_back(lep2[2]);
				       muon2.mass.push_back(lep2[3]);

				       // Muon3
				       muon3.pT.push_back(lep3[0]);
				       muon3.eta.push_back(lep3[1]);
				       muon3.phi.push_back(lep3[2]);
				       muon3.mass.push_back(lep3[3]);

				       // Muon4
				       muon4.pT.push_back(lep4[0]);
				       muon4.eta.push_back(lep4[1]);
				       muon4.phi.push_back(lep4[2]);
				       muon4.mass.push_back(lep4[3]);


				       // Z mass distribution for each muon pair
				       mZ12 = (mu1 + mu2).M();
				       mZ34 = (mu3 + mu4).M();

				       mz.mZ_12.push_back(mZ12);
				       mz.mZ_34.push_back(mZ34);
				      

				       // On-Shell Z
                                       Z1 = mu3 + mu4;
				       mZ1   = Z1.M();   // invariant mass of muons 1&2
		                       ptZ1  = Z1.Pt();
                                       etaZ1 = Z1.Eta();
		                       phiZ1 = Z1.Phi();


                                       // Off-Shell Z
				       Z2 = mu1 + mu2;
                                       mZ2   = Z2.M();   // invariant mass of muons 1&2
		                       ptZ2  = Z2.Pt();
                                       etaZ2 = Z2.Eta();
		                       phiZ2 = Z2.Phi();
				       

                                       // Filling Z1 Z2 MASS distributions for Signal
				       za_sig.Za_mass.push_back(mZ1);
				       zb_sig.Zb_mass.push_back(mZ2);
				       

				    } // end DR
				  } // end mass of each mu pair > 3GeV
				} // end mComb_1 is m34
				
			      } // end 1st category


			      if ( catgr3 ){  // 1324

				NEvents[10]++;

				cout << " << Working with Category3, muons 1324 >>" << endl;

				if ( mComb_3 == m13 ){    // 1324

				  cout << " mass_Comb_3 = mass mu13 = " << mComb_3 << endl;

				  if ( m13 > 3. && m24 > 3.) {    // selection on mass of each mu pair

				    if ( DR_Tight_mu13 > 0.07 && DR_Tight_mu24 > 0.07 ) {     // selection on DR

				      // NEvents[16]++;

				      // start filling leptons vectors 1423

				      // Lepton1
				      lep1.push_back(pt_mu1);
				      lep1.push_back(eta_mu1);
				      lep1.push_back(phi_mu1);
				      lep1.push_back( mu1.M() );

				      // Lepton2
				      lep2.push_back(pt_mu3);
				      lep2.push_back(eta_mu3);
				      lep2.push_back(phi_mu3);
				      lep2.push_back( mu3.M() );

				      // Lepton3
				      lep3.push_back(pt_mu2);
				      lep3.push_back(eta_mu2);
				      lep3.push_back(phi_mu2);
				      lep3.push_back( mu2.M() );

				      // Lepton4
				      lep4.push_back(pt_mu4);
				      lep4.push_back(eta_mu4);
				      lep4.push_back(phi_mu4);
				      lep4.push_back( mu4.M() );

				      // Filling 4 Muons of signal Branches

				      // Muon1
				      muon1.pT.push_back(lep1[0]);
				      muon1.eta.push_back(lep1[1]);
				      muon1.phi.push_back(lep1[2]);
				      muon1.mass.push_back(lep1[3]);

				      // Muon2
				      muon2.pT.push_back(lep2[0]);
				      muon2.eta.push_back(lep2[1]);
				      muon2.phi.push_back(lep2[2]);
				      muon2.mass.push_back(lep2[3]);

				      // Muon3
				      muon3.pT.push_back(lep3[0]);
				      muon3.eta.push_back(lep3[1]);
				      muon3.phi.push_back(lep3[2]);
				      muon3.mass.push_back(lep3[3]);

				      // Muon4
				      muon4.pT.push_back(lep4[0]);
				      muon4.eta.push_back(lep4[1]);
				      muon4.phi.push_back(lep4[2]);
				      muon4.mass.push_back(lep4[3]);

				      
				      // Z mass distribution for each muon pair
				      mZ13 = (mu1 + mu3).M();
				      mZ24 = (mu2 + mu4).M();

				      mz.mZ_13.push_back(mZ13);
				      mz.mZ_24.push_back(mZ24);
				      

				      // On-Shell Z Signal
                                      Z1 = mu1 + mu3;
				      mZ1   = Z1.M();   // invariant mass of muons 1&2
		                      ptZ1  = Z1.Pt();
                                      etaZ1 = Z1.Eta();
		                      phiZ1 = Z1.Phi();


                                      // Off-Shell Z
				      Z2 = mu2 + mu4;
                                      mZ2   = Z2.M();   // invariant mass of muons 1&2
		                      ptZ2  = Z2.Pt();
                                      etaZ2 = Z2.Eta();
		                      phiZ2 = Z2.Phi();

				      // Filling Z1 Z2 MASS distributions for Signal
				      za_sig.Za_mass.push_back(mZ1);
				      zb_sig.Zb_mass.push_back(mZ2);

				    } // end DR

				  } // end mass of each mu pair > 3 GeV

				} // end if mComb_3 is m13

				else{  // 2413

				  cout << " mass_Comb_3 = mass mu24 = " << mComb_3 << endl;

				  if ( m13 > 3. && m24 > 3.) {    // selection on mass of each mu pair

				    if ( DR_Tight_mu13 > 0.07 && DR_Tight_mu24 > 0.07 ) {     // selection on DR

				      // NEvents[17]++;

				       // Lepton1
				       lep1.push_back(pt_mu2);
				       lep1.push_back(eta_mu2);
				       lep1.push_back(phi_mu2);
				       lep1.push_back( mu2.M() );

				       // Lepton2
				       lep2.push_back(pt_mu4);
				       lep2.push_back(eta_mu4);
				       lep2.push_back(phi_mu4);
				       lep2.push_back( mu4.M() );

				       // Lepton3
				       lep3.push_back(pt_mu1);
				       lep3.push_back(eta_mu1);
				       lep3.push_back(phi_mu1);
				       lep3.push_back( mu1.M() );

				       // Lepton4
				       lep4.push_back(pt_mu3);
				       lep4.push_back(eta_mu3);
				       lep4.push_back(phi_mu3);
				       lep4.push_back( mu3.M() );


				       // Filling 4 Muons of signal Branches
				       
				       // Muon1
				       muon1.pT.push_back(lep1[0]);
				       muon1.eta.push_back(lep1[1]);
				       muon1.phi.push_back(lep1[2]);
				       muon1.mass.push_back(lep1[3]);

				       // Muon2
				       muon2.pT.push_back(lep2[0]);
				       muon2.eta.push_back(lep2[1]);
				       muon2.phi.push_back(lep2[2]);
				       muon2.mass.push_back(lep2[3]);

				       // Muon3
				       muon3.pT.push_back(lep3[0]);
				       muon3.eta.push_back(lep3[1]);
				       muon3.phi.push_back(lep3[2]);
				       muon3.mass.push_back(lep3[3]);

				       // Muon4
				       muon4.pT.push_back(lep4[0]);
				       muon4.eta.push_back(lep4[1]);
				       muon4.phi.push_back(lep4[2]);
				       muon4.mass.push_back(lep4[3]);


				       // Z mass distribution for each muon pair
				       mZ13 = (mu1 + mu3).M();
				       mZ24 = (mu2 + mu4).M();

				       mz.mZ_13.push_back(mZ13);
				       mz.mZ_24.push_back(mZ24);
				      

				       // On-Shell Z
                                       Z1 = mu2 + mu4;
				       mZ1   = Z1.M();   // invariant mass of muons 1&2
		                       ptZ1  = Z1.Pt();
                                       etaZ1 = Z1.Eta();
		                       phiZ1 = Z1.Phi();


                                       // Off-Shell Z
				       Z2 = mu1 + mu3;
                                       mZ2   = Z2.M();   // invariant mass of muons 1&2
		                       ptZ2  = Z2.Pt();
                                       etaZ2 = Z2.Eta();
		                       phiZ2 = Z2.Phi();
				       

                                       // Filling Z1 Z2 MASS distributions for Signal
				       za_sig.Za_mass.push_back(mZ1);
				       zb_sig.Zb_mass.push_back(mZ2);
				       

				    } // end DR
				  } // end mass of each mu pair > 3GeV
				} // end mComb_3 is m24
				
			      } // end 3rd category

			      // } // end if on muons pT 
					   // } // end if dz                
			      //  } // end if dxy
	         } // end mu 14,23 same charge

	      } // end isCase1

	     
	     //===========================================================================//
             //                                CASE2 Fullfilled                           //
	     //===========================================================================//

	     if ( isCase2 ){    // we have only here 1324, 1423

	       // if ( (dxy_mu1 < 0.12) && (dxy_mu2 < 0.12)
	       // && (dxy_mu3 < 0.12) && (dxy_mu4 < 0.12) ){

	       //  if ( (dz_mu1 < 0.47) && (dz_mu2 < 0.47)
	       //	&& (dz_mu3 < 0.47) && (dz_mu4 < 0.47) ){

	       //  if ( ( pt_mu1 > 20. ) && ( pt_mu2 > 10. ) ){
		     // if ( ( pt_mu1 > 35. ) && ( pt_mu2 > 20. ) ){
			
                          m13 = (mu1 + mu3).M();
                          m24 = (mu2 + mu4).M();
			  m14 = (mu1 + mu4).M();
                          m23 = (mu2 + mu3).M();

			  if ( m13 > m24 ) mComb_3 = m13;
		          if ( m13 < m24 ) mComb_3 = m24;

	         	  if ( m14 > m23 ) mComb_2 = m14;
			  if ( m14 < m23 ) mComb_2 = m23;

			  if ( mComb_3 > mComb_2 ) {

                             cout << "< Fullfill Category 3 >  " << "mass_Comb_3  " << mComb_3 << "  mass_Comb_2  " << mComb_2 << endl; 
	          	     catgr3 = true;
		             cout << " is Category 3 " << catgr3 << endl;

			  }

                          if ( mComb_3 < mComb_2 ) {

                             cout << "< Fullfill Category 2 >  " << "mass_Comb_3  " << mComb_3 << "  mass_Comb_2  " << mComb_2 << endl; 
	          	     catgr2 = true;
		             cout << " is Category 2 " << catgr2 << endl;

			  }

                          if ( catgr3 ){  // 1324

				NEvents[10]++;

				cout << " << Working with Category3, muons 1324 >>" << endl;

				if ( mComb_3 == m13 ){    // 1324

				  cout << " mass_Comb_3 = mass mu13 = " << mComb_3 << endl;

				  if ( m13 > 3. && m24 > 3.) {    // selection on mass of each mu pair

				    if ( DR_Tight_mu13 > 0.07 && DR_Tight_mu24 > 0.07 ) {     // selection on DR

				      // NEvents[16]++;

				      // start filling leptons vectors 1423

				      // Lepton1
				      lep1.push_back(pt_mu1);
				      lep1.push_back(eta_mu1);
				      lep1.push_back(phi_mu1);
				      lep1.push_back( mu1.M() );

				      // Lepton2
				      lep2.push_back(pt_mu3);
				      lep2.push_back(eta_mu3);
				      lep2.push_back(phi_mu3);
				      lep2.push_back( mu3.M() );

				      // Lepton3
				      lep3.push_back(pt_mu2);
				      lep3.push_back(eta_mu2);
				      lep3.push_back(phi_mu2);
				      lep3.push_back( mu2.M() );

				      // Lepton4
				      lep4.push_back(pt_mu4);
				      lep4.push_back(eta_mu4);
				      lep4.push_back(phi_mu4);
				      lep4.push_back( mu4.M() );

				      // Filling 4 Muons of signal Branches

				      // Muon1
				      muon1.pT.push_back(lep1[0]);
				      muon1.eta.push_back(lep1[1]);
				      muon1.phi.push_back(lep1[2]);
				      muon1.mass.push_back(lep1[3]);

				      // Muon2
				      muon2.pT.push_back(lep2[0]);
				      muon2.eta.push_back(lep2[1]);
				      muon2.phi.push_back(lep2[2]);
				      muon2.mass.push_back(lep2[3]);

				      // Muon3
				      muon3.pT.push_back(lep3[0]);
				      muon3.eta.push_back(lep3[1]);
				      muon3.phi.push_back(lep3[2]);
				      muon3.mass.push_back(lep3[3]);

				      // Muon4
				      muon4.pT.push_back(lep4[0]);
				      muon4.eta.push_back(lep4[1]);
				      muon4.phi.push_back(lep4[2]);
				      muon4.mass.push_back(lep4[3]);

				      
				      // Z mass distribution for each muon pair
				      mZ13 = (mu1 + mu3).M();
				      mZ24 = (mu2 + mu4).M();

				      mz.mZ_13.push_back(mZ13);
				      mz.mZ_24.push_back(mZ24);
				      

				      // On-Shell Z Signal
                                      Z1 = mu1 + mu3;
				      mZ1   = Z1.M();   // invariant mass of muons 1&2
		                      ptZ1  = Z1.Pt();
                                      etaZ1 = Z1.Eta();
		                      phiZ1 = Z1.Phi();


                                      // Off-Shell Z
				      Z2 = mu2 + mu4;
                                      mZ2   = Z2.M();   // invariant mass of muons 1&2
		                      ptZ2  = Z2.Pt();
                                      etaZ2 = Z2.Eta();
		                      phiZ2 = Z2.Phi();

				      // Filling Z1 Z2 MASS distributions for Signal
				      za_sig.Za_mass.push_back(mZ1);
				      zb_sig.Zb_mass.push_back(mZ2);

				    } // end DR

				  } // end mass of each mu pair > 3 GeV

				} // end if mComb_3 is m13

				else{  // 2413

				  cout << " mass_Comb_3 = mass mu24 = " << mComb_3 << endl;

				  if ( m13 > 3. && m24 > 3.) {    // selection on mass of each mu pair

				    if ( DR_Tight_mu13 > 0.07 && DR_Tight_mu24 > 0.07 ) {     // selection on DR

				      // NEvents[17]++;

				       // Lepton1
				       lep1.push_back(pt_mu2);
				       lep1.push_back(eta_mu2);
				       lep1.push_back(phi_mu2);
				       lep1.push_back( mu2.M() );

				       // Lepton2
				       lep2.push_back(pt_mu4);
				       lep2.push_back(eta_mu4);
				       lep2.push_back(phi_mu4);
				       lep2.push_back( mu4.M() );

				       // Lepton3
				       lep3.push_back(pt_mu1);
				       lep3.push_back(eta_mu1);
				       lep3.push_back(phi_mu1);
				       lep3.push_back( mu1.M() );

				       // Lepton4
				       lep4.push_back(pt_mu3);
				       lep4.push_back(eta_mu3);
				       lep4.push_back(phi_mu3);
				       lep4.push_back( mu3.M() );


				       // Filling 4 Muons of signal Branches
				       
				       // Muon1
				       muon1.pT.push_back(lep1[0]);
				       muon1.eta.push_back(lep1[1]);
				       muon1.phi.push_back(lep1[2]);
				       muon1.mass.push_back(lep1[3]);

				       // Muon2
				       muon2.pT.push_back(lep2[0]);
				       muon2.eta.push_back(lep2[1]);
				       muon2.phi.push_back(lep2[2]);
				       muon2.mass.push_back(lep2[3]);

				       // Muon3
				       muon3.pT.push_back(lep3[0]);
				       muon3.eta.push_back(lep3[1]);
				       muon3.phi.push_back(lep3[2]);
				       muon3.mass.push_back(lep3[3]);

				       // Muon4
				       muon4.pT.push_back(lep4[0]);
				       muon4.eta.push_back(lep4[1]);
				       muon4.phi.push_back(lep4[2]);
				       muon4.mass.push_back(lep4[3]);


				       // Z mass distribution for each muon pair
				       mZ13 = (mu1 + mu3).M();
				       mZ24 = (mu2 + mu4).M();

				       mz.mZ_13.push_back(mZ13);
				       mz.mZ_24.push_back(mZ24);
				      

				       // On-Shell Z
                                       Z1 = mu2 + mu4;
				       mZ1   = Z1.M();   // invariant mass of muons 1&2
		                       ptZ1  = Z1.Pt();
                                       etaZ1 = Z1.Eta();
		                       phiZ1 = Z1.Phi();


                                       // Off-Shell Z
				       Z2 = mu1 + mu3;
                                       mZ2   = Z2.M();   // invariant mass of muons 1&2
		                       ptZ2  = Z2.Pt();
                                       etaZ2 = Z2.Eta();
		                       phiZ2 = Z2.Phi();
				       

                                       // Filling Z1 Z2 MASS distributions for Signal
				       za_sig.Za_mass.push_back(mZ1);
				       zb_sig.Zb_mass.push_back(mZ2);
				       

				    } // end DR
				  } // end mass of each mu pair > 3GeV
				} // end mComb_3 is m24
				
			      } // end 3rd category


			      if ( catgr2 ){  // 1423

				NEvents[9]++;

				cout << " << Working with Category2, muons 1423 >>" << endl;

				if ( mComb_2 == m14 ){    // 1423

                                  cout << " mass_Comb_2 = mass mu14 = " << mComb_2 << endl;
				  
				  if ( m14 > 3. && m23 > 3.) {    // selection on mass of each mu pair

				    if ( DR_Tight_mu14 > 0.07 && DR_Tight_mu23 > 0.07 ) {     // selection on DR

				      //NEvents[9]++;

				      // start filling leptons vectors 1423

				      // Lepton1
				      lep1.push_back(pt_mu1);
				      lep1.push_back(eta_mu1);
				      lep1.push_back(phi_mu1);
				      lep1.push_back( mu1.M() );

				      // Lepton2
				      lep2.push_back(pt_mu4);
				      lep2.push_back(eta_mu4);
				      lep2.push_back(phi_mu4);
				      lep2.push_back( mu4.M() );

				      // Lepton3
				      lep3.push_back(pt_mu2);
				      lep3.push_back(eta_mu2);
				      lep3.push_back(phi_mu2);
				      lep3.push_back( mu2.M() );

				      // Lepton4
				      lep4.push_back(pt_mu3);
				      lep4.push_back(eta_mu3);
				      lep4.push_back(phi_mu3);
				      lep4.push_back( mu3.M() );

				      // Filling 4 Muons of signal Branches

				      // Muon1
				      muon1.pT.push_back(lep1[0]);
				      muon1.eta.push_back(lep1[1]);
				      muon1.phi.push_back(lep1[2]);
				      muon1.mass.push_back(lep1[3]);

				      // Muon2
				      muon2.pT.push_back(lep2[0]);
				      muon2.eta.push_back(lep2[1]);
				      muon2.phi.push_back(lep2[2]);
				      muon2.mass.push_back(lep2[3]);

				      // Muon3
				      muon3.pT.push_back(lep3[0]);
				      muon3.eta.push_back(lep3[1]);
				      muon3.phi.push_back(lep3[2]);
				      muon3.mass.push_back(lep3[3]);

				      // Muon4
				      muon4.pT.push_back(lep4[0]);
				      muon4.eta.push_back(lep4[1]);
				      muon4.phi.push_back(lep4[2]);
				      muon4.mass.push_back(lep4[3]);

				      
				      // Z mass distribution for each muon pair
				      mZ14 = (mu1 + mu4).M();
				      mZ23 = (mu2 + mu3).M();

				      mz.mZ_14.push_back(mZ14);
				      mz.mZ_23.push_back(mZ23);
				      

				      // On-Shell Z Signal
                                      Z1 = mu1 + mu4;
				      mZ1   = Z1.M();   // invariant mass of muons 1&2
		                      ptZ1  = Z1.Pt();
                                      etaZ1 = Z1.Eta();
		                      phiZ1 = Z1.Phi();


                                      // Off-Shell Z
				      Z2 = mu2 + mu3;
                                      mZ2   = Z2.M();   // invariant mass of muons 1&2
		                      ptZ2  = Z2.Pt();
                                      etaZ2 = Z2.Eta();
		                      phiZ2 = Z2.Phi();

				      // Filling Z1 Z2 MASS distributions for Signal
				      za_sig.Za_mass.push_back(mZ1);
				      zb_sig.Zb_mass.push_back(mZ2);

				    } // end DR

				  } // end mass of each mu pair > 3 GeV

				} // end if mComb_2 is m14

				else{  // 2314

				   cout << " mass_Comb_2 = mass mu23 = " << mComb_2 << endl;

				  if ( m14 > 3. && m23 > 3.) {    // selection on mass of each mu pair

				    if ( DR_Tight_mu14 > 0.07 && DR_Tight_mu23 > 0.07 ) {     // selection on DR

				      //   NEvents[10]++;

				       // Lepton1
				       lep1.push_back(pt_mu2);
				       lep1.push_back(eta_mu2);
				       lep1.push_back(phi_mu2);
				       lep1.push_back( mu2.M() );

				       // Lepton2
				       lep2.push_back(pt_mu3);
				       lep2.push_back(eta_mu3);
				       lep2.push_back(phi_mu3);
				       lep2.push_back( mu3.M() );

				       // Lepton3
				       lep3.push_back(pt_mu1);
				       lep3.push_back(eta_mu1);
				       lep3.push_back(phi_mu1);
				       lep3.push_back( mu1.M() );

				       // Lepton4
				       lep4.push_back(pt_mu4);
				       lep4.push_back(eta_mu4);
				       lep4.push_back(phi_mu4);
				       lep4.push_back( mu4.M() );


				       // Filling 4 Muons of signal Branches
				       
				       // Muon1
				       muon1.pT.push_back(lep1[0]);
				       muon1.eta.push_back(lep1[1]);
				       muon1.phi.push_back(lep1[2]);
				       muon1.mass.push_back(lep1[3]);

				       // Muon2
				       muon2.pT.push_back(lep2[0]);
				       muon2.eta.push_back(lep2[1]);
				       muon2.phi.push_back(lep2[2]);
				       muon2.mass.push_back(lep2[3]);

				       // Muon3
				       muon3.pT.push_back(lep3[0]);
				       muon3.eta.push_back(lep3[1]);
				       muon3.phi.push_back(lep3[2]);
				       muon3.mass.push_back(lep3[3]);

				       // Muon4
				       muon4.pT.push_back(lep4[0]);
				       muon4.eta.push_back(lep4[1]);
				       muon4.phi.push_back(lep4[2]);
				       muon4.mass.push_back(lep4[3]);


				       // Z mass distribution for each muon pair
				       mZ14 = (mu1 + mu4).M();
				       mZ23 = (mu2 + mu3).M();

				       mz.mZ_14.push_back(mZ14);
				       mz.mZ_23.push_back(mZ23);
				      

				       // On-Shell Z
                                       Z1 = mu2 + mu3;
				       mZ1   = Z1.M();   // invariant mass of muons 1&2
		                       ptZ1  = Z1.Pt();
                                       etaZ1 = Z1.Eta();
		                       phiZ1 = Z1.Phi();


                                       // Off-Shell Z
				       Z2 = mu1 + mu4;
                                       mZ2   = Z2.M();   // invariant mass of muons 1&2
		                       ptZ2  = Z2.Pt();
                                       etaZ2 = Z2.Eta();
		                       phiZ2 = Z2.Phi();
				       

                                       // Filling Z1 Z2 MASS distributions for Signal
				       za_sig.Za_mass.push_back(mZ1);
				       zb_sig.Zb_mass.push_back(mZ2);
				       

				    } // end DR
				  } // end mass of each mu pair > 3GeV
				} // end mComb_2 is m23
				
			      } // end 2nd category
		      //}
	    	//  }
	      // }   
	     } // end case2 
	       

        
	     // Specifiy Z1,Z2 mass windows

	     if ( (mZ1 > 80.) && (mZ1 < 100.) ){

		if ( (mZ2 > 20.) && (mZ2 < 40.) ){

		    // if ( (mZ1 > 40.) && (mZ1 < 120.) ){

		    // if ( (mZ2 > 12.) && (mZ2 < 80.) ){

		    NEvents[11]++;

				    
				   
		   //====================================================//
                   //              Start b1, b2 jets Selections          //
                   //====================================================//                                        

				    
		   int nbjets = 0;            // total nb of b-jets found per event
                                       
                   // Loop overall Jets
                            
		   for (Int_t i = 0; i < Jet_size; i++){
	
	               int bj_indx = i;
	                            
	    	       // Jet_BTag stores the bit value for each jet

                       //------------------------------------
	 	       //   WP       BitNumber    Bit Values
                       //------------------------------------
		       // Loose          0          1,3,5,7
		       // Medium         1          2,3,6,7
		       // Tight          2          4,5,6,7
		       //------------------------------------

		       Float_t j_pt = Jet_PT[i];
	               Float_t j_eta = Jet_Eta[i];
	               Float_t j_phi = Jet_Phi[i];
		       UInt_t jet_bTag = Jet_BTag[i];

		       jets.jet_pt.push_back(j_pt);
                       jets.jet_eta.push_back(j_eta);
                       jets.jet_phi.push_back(j_phi);
                                        
	    	       // Select loose b-jets 
                      // Bool_t isLoose = ( jet_bTag & (1 << 0) );
                       Bool_t isMed = ( jet_bTag & (1 << 1) );

		       //if ( isLoose == 1){
			   if ( isMed == 1){	   

			  jets.loose_bjet_pt.push_back(j_pt);
                          jets.loose_bjet_eta.push_back(j_eta);
                          jets.loose_bjet_phi.push_back(j_phi);
					  
		          bjet_indx.push_back(bj_indx);
			  nbjets++;
									
		       }  // end isLoose 
	           }  // end loop overall jets

		   jets.jet_size.push_back(Jet_size);
		   jets.loose_bjet_size.push_back(nbjets);
      
					 
		   // 10th Selection: on having at least 2 b-jets/event 
 
		   if ( nbjets > 1 ){ 
								 
		     NEvents[12]++;
								 
	    	     cout << "========================================" << endl;
		     cout << " number of b-jets/event =  " << nbjets    << endl;
		     cout << "========================================" << endl;
		     cout << " b-jets vector size is  ;  " << bjet_indx.size() << endl;
		     cout << "========================================" << endl;
							 		 
		     // Loop overall loose b-jets/event 

		     for ( Int_t i = 0; i < bjet_indx.size(); i++ ) { 
									 
			 int bj_indx_AfterSelec = bjet_indx[i];

			 float bjet_pt = Jet_PT[bj_indx_AfterSelec];
			 float bjet_eta = Jet_Eta[bj_indx_AfterSelec];
			 float abs_bjet_eta = abs(bjet_eta);
			 float bjet_phi = Jet_Phi[bj_indx_AfterSelec];
				    			    
			 // get DR between loose b-jet & 4TightMuons

			 double DEta_b_mu1_sqr = TMath::Power(( bjet_eta - lep1[1] ), 2);
			 double DEta_b_mu2_sqr = TMath::Power(( bjet_eta - lep2[1] ), 2);
			 double DEta_b_mu3_sqr = TMath::Power(( bjet_eta - lep3[1] ), 2);
			 double DEta_b_mu4_sqr = TMath::Power(( bjet_eta - lep4[1] ), 2);
			 double DPhi_b_mu1_sqr = TMath::Power(( bjet_phi - lep1[2]), 2);
			 double DPhi_b_mu2_sqr = TMath::Power(( bjet_phi - lep2[2]), 2);
			 double DPhi_b_mu3_sqr = TMath::Power(( bjet_phi - lep3[2]), 2);
			 double DPhi_b_mu4_sqr = TMath::Power(( bjet_phi - lep4[2]), 2);
					        			 
			 double DR_b_mu1 = TMath::Sqrt( DEta_b_mu1_sqr + DPhi_b_mu1_sqr );
			 double DR_b_mu2 = TMath::Sqrt( DEta_b_mu2_sqr + DPhi_b_mu2_sqr );
			 double DR_b_mu3 = TMath::Sqrt( DEta_b_mu3_sqr + DPhi_b_mu3_sqr );
			 double DR_b_mu4 = TMath::Sqrt( DEta_b_mu4_sqr + DPhi_b_mu4_sqr );


			 // 11th Selection on pT & |eta| of loose b-jet

			 if ( (bjet_pt > 20.) && (abs_bjet_eta < 2.7) ){

			    //NEvents[13]++;

 			    // 12th Selection: on DR(bjet_loose, mu_tight)
 
			   if ( ( DR_b_mu1 > 0.4 ) && ( DR_b_mu2 > 0.4 )
				&& ( DR_b_mu3 > 0.4 ) && ( DR_b_mu4 > 0.4 ) ){

			      // NEvents[11]++;

			      // save b-jet index pass all above selections 

			      bjet_indx_AfterSel.push_back(bj_indx_AfterSelec);

		              /*drLooseBjet4Mu.DR_LooseBJet_TightMu1.clear();
                              drLooseBjet4Mu.DR_LooseBJet_TightMu2.clear();
                              drLooseBjet4Mu.DR_LooseBJet_TightMu3.clear();
                              drLooseBjet4Mu.DR_LooseBJet_TightMu4.clear();
                              */
						   
		           } // end if DR

			  } // end if pT & abs eta
				    
			} // end loop over bjet_index.size()

		   } // end nbjets > 1
                                   
                                
		   // 13th Selection: check if  we still have at least 2 bjets after above selections

		   cout << "b-jets nb after all selections :- " << bjet_indx_AfterSel.size() << endl;
		   
	           if ( bjet_indx_AfterSel.size() > 1 ) {
					  
		      NEvents[18]++;  // nb of events passed selections till bjet pt, eta, DR selections

		      // cout << "bjet_indx_AfterSel size is : " <<  bjet_indx_AfterSel.size() << endl; 
                                      
											
	   	      // Select the signal 2 b-jets as the two with highest pT, no scores in Delphes

		      int b1_idx = bjet_indx_AfterSel[0];
		      int b2_idx = bjet_indx_AfterSel[1];

					  
		      // b1 jet of signal kinematics
		      Float_t b1_pt   = Jet_PT[b1_idx];
		      Float_t b1_eta  = Jet_Eta[b1_idx];
                      Float_t b1_phi  = Jet_Phi[b1_idx];
		      Int_t b1_charge = Jet_Charge[b1_idx];


		      // b2 jet of signal kinematics
		      Float_t b2_pt   = Jet_PT[b2_idx];
		      Float_t b2_eta  = Jet_Eta[b2_idx];
                      Float_t b2_phi  = Jet_Phi[b2_idx];
		      Int_t b2_charge = Jet_Charge[b2_idx];

					  
		      v_b1_jet_Charge.push_back(b1_charge);
		      v_b2_jet_Charge.push_back(b2_charge);
				    						    
	  	      cout << "2 B-jets of signal are ( " << b1_idx << ", " << b2_idx 
			   << " ) with Highest PT values " << b1_pt << ", " << b2_pt
			   << " GeV respectively!" << endl;
							          
                      // 14th Selection: on Total charge of 2 bjets 
					  
		      if ( b1_charge + b2_charge == 0){

			 NEvents[14]++;
		   	      
	                // Set TLorentzVectors for selected 2 b-jets of signal

				             b1.SetPtEtaPhiM( Jet_PT[b1_idx], Jet_Eta[b1_idx], Jet_Phi[b1_idx], Jet_Mass[b1_idx] );		                     
				             b2.SetPtEtaPhiM( Jet_PT[b2_idx], Jet_Eta[b2_idx], Jet_Phi[b2_idx], Jet_Mass[b2_idx] );
							     
							     
				             Float_t bj1_mass =  b1.M();
			                     Float_t bj1_pt   =  b1.Pt();
				             Float_t bj1_eta  =  b1.Eta();
				             Float_t bj1_phi  =  b1.Phi();
							     
				             Float_t bj2_mass =  b2.M();
				             Float_t bj2_pt   =  b2.Pt();
				             Float_t bj2_eta  =  b2.Eta();
				             Float_t bj2_phi  =  b2.Phi();

					     
                                             // DR between 2 bjets of signal b1,b2
					     Float_t DeltaEta_b1b2_sqr = TMath::Power( ( bj1_eta - bj2_eta ), 2); 
				             Float_t DeltaPhi_b1b2_sqr = TMath::Power( ( bj1_phi - bj2_phi ), 2); 
				             Float_t DR_b1b2 = TMath::Sqrt( DeltaEta_b1b2_sqr + DeltaPhi_b1b2_sqr );

				             dr_b1b2.push_back(DR_b1b2);

					     
					     // DR ( b1 of signal, 4TightMuons of signal )
					     double DEta_b1_mu1_sqr = TMath::Power((bj1_eta - eta_mu1), 2);
				             double DEta_b1_mu2_sqr = TMath::Power((bj1_eta - eta_mu2), 2);
				             double DEta_b1_mu3_sqr = TMath::Power((bj1_eta - eta_mu3), 2);
			                     double DEta_b1_mu4_sqr = TMath::Power((bj1_eta - eta_mu4), 2);
			                     double DPhi_b1_mu1_sqr = TMath::Power((bj1_phi - phi_mu1), 2);
			                     double DPhi_b1_mu2_sqr = TMath::Power((bj1_phi - phi_mu2), 2);
				             double DPhi_b1_mu3_sqr = TMath::Power((bj1_phi - phi_mu3), 2);
				             double DPhi_b1_mu4_sqr = TMath::Power((bj1_phi - phi_mu4), 2);
					        			 
				             double dR_b1_mu1 = TMath::Sqrt( DEta_b1_mu1_sqr + DPhi_b1_mu1_sqr );
				             double dR_b1_mu2 = TMath::Sqrt( DEta_b1_mu2_sqr + DPhi_b1_mu2_sqr );
				             double dR_b1_mu3 = TMath::Sqrt( DEta_b1_mu3_sqr + DPhi_b1_mu3_sqr );
				             double dR_b1_mu4 = TMath::Sqrt( DEta_b1_mu4_sqr + DPhi_b1_mu4_sqr );

					     b1_sig.DR_b1_mu1.push_back(dR_b1_mu1);
                                             b1_sig.DR_b1_mu2.push_back(dR_b1_mu2);
                                             b1_sig.DR_b1_mu3.push_back(dR_b1_mu3);
                                             b1_sig.DR_b1_mu4.push_back(dR_b1_mu4);


					     // DR ( b2 of signal, 4TightMuons of signal )
					     double DEta_b2_mu1_sqr = TMath::Power((bj2_eta - eta_mu1), 2);
				             double DEta_b2_mu2_sqr = TMath::Power((bj2_eta - eta_mu2), 2);
				             double DEta_b2_mu3_sqr = TMath::Power((bj2_eta - eta_mu3), 2);
			                     double DEta_b2_mu4_sqr = TMath::Power((bj2_eta - eta_mu4), 2);
			                     double DPhi_b2_mu1_sqr = TMath::Power((bj2_phi - phi_mu1), 2);
			                     double DPhi_b2_mu2_sqr = TMath::Power((bj2_phi - phi_mu2), 2);
				             double DPhi_b2_mu3_sqr = TMath::Power((bj2_phi - phi_mu3), 2);
				             double DPhi_b2_mu4_sqr = TMath::Power((bj2_phi - phi_mu4), 2);
					        			 
				             double dR_b2_mu1 = TMath::Sqrt( DEta_b2_mu1_sqr + DPhi_b2_mu1_sqr );
				             double dR_b2_mu2 = TMath::Sqrt( DEta_b2_mu2_sqr + DPhi_b2_mu2_sqr );
				             double dR_b2_mu3 = TMath::Sqrt( DEta_b2_mu3_sqr + DPhi_b2_mu3_sqr );
				             double dR_b2_mu4 = TMath::Sqrt( DEta_b2_mu4_sqr + DPhi_b2_mu4_sqr );

					     b2_sig.DR_b2_mu1.push_back(dR_b2_mu1);
                                             b2_sig.DR_b2_mu2.push_back(dR_b2_mu2);
                                             b2_sig.DR_b2_mu3.push_back(dR_b2_mu3);
                                             b2_sig.DR_b2_mu4.push_back(dR_b2_mu4);

					     
                                             // filling 4Muons kinematics after full selection
                                             tight_mu.muon1_pt.push_back(  MuonTight_PT[mu1_idx]  );
                                             tight_mu.muon2_pt.push_back(  MuonTight_PT[mu2_idx]  );
                                             tight_mu.muon3_pt.push_back(  MuonTight_PT[mu3_idx]  );
                                             tight_mu.muon4_pt.push_back(  MuonTight_PT[mu4_idx]  );
                                             tight_mu.muon1_eta.push_back( MuonTight_Eta[mu1_idx] );
                                             tight_mu.muon2_eta.push_back( MuonTight_Eta[mu2_idx] );
                                             tight_mu.muon3_eta.push_back( MuonTight_Eta[mu3_idx] );
                                             tight_mu.muon4_eta.push_back( MuonTight_Eta[mu4_idx] );
                                             tight_mu.muon1_phi.push_back( MuonTight_Phi[mu1_idx] );
                                             tight_mu.muon2_phi.push_back( MuonTight_Phi[mu2_idx] );
                                             tight_mu.muon3_phi.push_back( MuonTight_Phi[mu3_idx] );
                                             tight_mu.muon4_phi.push_back( MuonTight_Phi[mu4_idx] );
				             drMuonsTight.v_DR_mu1mu2.push_back(DR_Tight_mu12);
                                             drMuonsTight.v_DR_mu3mu4.push_back(DR_Tight_mu34);
                                             tight_mu.tight_mu12_mass.push_back(m12);
                                             tight_mu.tight_mu34_mass.push_back(m34);


					     // filling  Za, Zb of signal Vectors after full selection
				             za_sig.Za_mass.push_back(mZ1);
			                     za_sig.Za_pt.push_back(ptZ1);
			                     za_sig.Za_eta.push_back(etaZ1);
			                     zb_sig.Zb_mass.push_back(mZ2);
			                     zb_sig.Zb_pt.push_back(ptZ2);
			                     zb_sig.Zb_eta.push_back(etaZ2);

				             ZZ->Fill(mZ2,mZ1);


                                             // filling b1, b2 of signal Vectors after full selection
				             b1_sig.b1_mass.push_back(bj1_mass);
				             b1_sig.b1_pt.push_back(bj1_pt);
				             b1_sig.b1_eta.push_back(bj1_eta);	            
				             b2_sig.b2_mass.push_back(bj2_mass);
				             b2_sig.b2_pt.push_back(bj2_pt);
				             b2_sig.b2_eta.push_back(bj2_eta);
							            

                                    
				             //========================================//
		                             //       Reconstruct h1 from Z1, Z2       //
		                             //========================================//

				             Float_t mh1_ZaZb, pt_h1_ZaZb, eta_h1_ZaZb, phi_h1_ZaZb;
				    
                                             h1 = Z1 + Z2;
                                             mh1_ZaZb = h1.M(); // get invariant mass of 4Muons
			                     pt_h1_ZaZb = h1.Pt();
                                             eta_h1_ZaZb = h1.Eta();
                                             phi_h1_ZaZb = h1.Phi();

				             smhiggs1.v_h1_mass.push_back(mh1_ZaZb);
				             smhiggs1.v_h1_pt.push_back(pt_h1_ZaZb);
				             smhiggs1.v_h1_eta.push_back(eta_h1_ZaZb);
                                 

				             //========================================//
		                             //       Reconstruct h2 from b1, b2       //
		                             //========================================//
							     
				             // TLorentzVector for 2nd SM higgs
				             h2 = b1 + b2;
							     
				             Float_t h2_b1b2_mass =  h2.M();
				             Float_t h2_b1b2_pt   =  h2.Pt();
				             Float_t h2_b1b2_eta  =  h2.Eta();
				             Float_t h2_b1b2_phi  =  h2.Phi();
				    
				             smhiggs2.v_h2_mass.push_back(h2_b1b2_mass);
				             smhiggs2.v_h2_pt.push_back(h2_b1b2_pt);
				             smhiggs2.v_h2_eta.push_back(h2_b1b2_eta);


					     // DR between 2 SM Higgs of signal
				             Float_t DeltaEta_h1h2_sqr = TMath::Power( ( eta_h1_ZaZb - h2_b1b2_eta ), 2); 
				             Float_t DeltaPhi_h1h2_sqr = TMath::Power( ( phi_h1_ZaZb - h2_b1b2_phi ), 2);
				             Float_t DR_h1h2 = TMath::Sqrt( DeltaEta_h1h2_sqr + DeltaPhi_h1h2_sqr );

                                             dr_h1h2.push_back(DR_h1h2);

				             hh->Fill(mh1_ZaZb, h2_b1b2_mass);
				             
					     
                                             //===========================================//
		                             //      Reconstruct BSM H from SM h1, h2     //
		                             //===========================================//
                
                                             // TLorentzVector for BSM Heavy Higgs 
				             H = h1 + h2;
                                 
                                             Float_t H_bb4Mu_mass =  H.M();
                                             Float_t H_bb4Mu_pt   =  H.Pt();
                                             Float_t H_bb4Mu_eta  =  H.Eta();
                                             Float_t H_bb4Mu_phi  =  H.Phi(); 
                                 
                                             heavyHiggs.v_H_mass.push_back(H_bb4Mu_mass); 
				             heavyHiggs.v_H_pt.push_back(H_bb4Mu_pt); 
				             heavyHiggs.v_H_eta.push_back(H_bb4Mu_eta);

				             /* H_mass->Fill( H_bb4Mu_mass, wt );
				             H_pt->Fill( H_bb4Mu_pt, wt );
				             H_eta->Fill( H_bb4Mu_eta, wt );
				             H_phi->Fill( H_bb4Mu_phi, wt );
                                             */
				    
                                             /*  H_mass->Fill(H_bb4Mu_mass);
				             H_pt->Fill(H_bb4Mu_pt);
				             H_eta->Fill(H_bb4Mu_eta);
				             H_phi->Fill(H_bb4Mu_phi); */
		     
				     				       
				            } // end if on 2b-jets charge 				       
				         } // end if bjet_indx_AfterSel.size() > 1 
					   } //mZ2
				    } // mZ1
			     } // pT lead & sublead
	          } // total charge of 4 mu = 0 
		   } // IP dz
	    } // IP d0     
      }  // // number of tight muons = 4 
           
   
       // new_tree->SetWeight(wt);
    
      cout << "Filling tree" << endl;
      new_tree->Fill();


   } // end loop overall events
   
   
   // Efficiency = nSelectedEvents/nentries;
   
   cout << "***Analysis Loop Ends!***" << endl; 
   
   cout << " Total Number of Events is : " << nentries << " Events" << endl;
   cout << "     " << endl;  
   
   cout << "================================================================================"    << endl;
   cout << "                         Cut                             |  NEvents PASS Cut    "    << endl;
   cout << "================================================================================"    << endl;
   cout << " Total Number of Generated Events                        |   " <<  nentries          << endl;
   cout << "================================================================================"    << endl;
   cout << " mu pT > 5GeV                                            |   " <<  NEvents[0]        << endl;
   cout << "================================================================================"    << endl;
   cout << " mu eta < 2.7                                            |   " <<  NEvents[1]        << endl;
   cout << "================================================================================"    << endl;
   cout << " number of muons/event = 4                               |   " <<  NEvents[2]        << endl;
   cout << "================================================================================"    << endl;
   cout << " IP dxy & dz for 4 muons                                 |   " <<  NEvents[3]        << endl; 
   cout << "================================================================================"    << endl;
   cout << " Total charge of 4 muons = 0                             |   " <<  NEvents[4]        << endl; 
   cout << "================================================================================"    << endl;
   cout << " pT lead & sublead > 20, 10 GeV                          |   " <<  NEvents[5]        << endl; 
   cout << "================================================================================"    << endl;
   cout << " isCase1 True                                            |   " <<  NEvents[6]        << endl; 
   cout << "================================================================================"    << endl;
   cout << " isCase2 True                                            |   " <<  NEvents[7]        << endl; 
   cout << "================================================================================"    << endl;
   cout << " Category 1; 1234                                        |   " <<  NEvents[8]        << endl; 
   cout << "================================================================================"    << endl;
   cout << " Category 2; 1423                                        |   " <<  NEvents[9]        << endl; 
   cout << "================================================================================"    << endl;
   cout << " Category 3; 1324                                        |   " <<  NEvents[10]       << endl; 
   cout << "================================================================================"    << endl;
   cout << " Z1,Z2 mass windows                                      |   " <<  NEvents[11]        << endl; 
   cout << "================================================================================"    << endl;
   cout << " nbjets > 1                                              |   " <<  NEvents[12]        << endl; 
   cout << "================================================================================"    << endl;
   cout << " bjet pT, eta, nbjets_Sel > 1                            |   " <<  NEvents[18]        << endl; 
   cout << "================================================================================"    << endl;
   // cout << " bjet pT, eta, nbjets_Sel > 1     new                   |   " <<  NEvents_13        << endl; 
   //cout << "================================================================================"    << endl;
   cout << " Charge of 2 b-jets of signal                            |   " <<  NEvents[14]        << endl; 
   cout << "================================================================================"    << endl;
   /* cout << " mass of mu pairs 12 & 34 > 3                            |   " <<  NEvents[7]        << endl; 
   cout << "================================================================================"    << endl;
   cout << "  Z1,Z2 mass windows                                     |   " <<  NEvents[8]        << endl; 
   cout << "================================================================================"    << endl;
   cout << " At least 2 loose b-jets/event  (mu selec + loose selec) |   " <<  NEvents[9]        << endl; 
   cout << "================================================================================"    << endl;
   cout << " bjet_pt > 20 && abs_bjet _eta < 2.7                     |   " <<  NEvents[10]       << endl; 
   cout << "================================================================================"    << endl;
   cout << " DR loose bjets & 4Muons of signal                       |   " <<  NEvents[11]       << endl; 
   cout << "================================================================================"    << endl;
   cout << " At least 2 b-jets/event  after full bjet selec          |   " <<  NEvents[12]       << endl; 
   cout << "================================================================================"    << endl;
   cout << " On Total charge of 2 bjets                              |   " <<  NEvents[13]       << endl; 
   cout << "================================================================================"    << endl;
   cout << " at least 2 b-jets/event  -  no cuts on b-jets           |   " <<  NEvents[14]       << endl; 
   cout << "================================================================================"    << endl;
   cout << " DeltaR(b-jet,lepton) of ZZ candidates > 0.3             |   " <<  NEvents[15]       << endl; 
   cout << "================================================================================"    << endl;
   cout << " b-jet pT > 20 GeV                                       |   " <<  NEvents[16]       << endl; 
   cout << "================================================================================"    << endl;
   cout << " absolute b-jet Eta < 2.4                                |   " <<  NEvents[17]       << endl; 
   cout << "================================================================================"    << endl;
   cout << " at least 2 b-jets/event after applying above selections |   " <<  NEvents[18]       << endl; 
   cout << "================================================================================"    << endl;   
   */
   
  /* cout << " Total Number of Events = " << nentries 
        << ", Number of Selected Events = " << nSelectedEvents
        << ", Efficiency = Nsel/Ngen = " << Efficiency << endl;
  */    
   
   cout << "Writing tree!" << endl;  
   new_tree->Write();

   /*  H_mass->Write();
   H_pt->Write();
   H_eta->Write();
   H_phi->Write();*/
     
   cout << "Writing output file! " << endl;
   out_sig->Write();
 
    //out_TT_CMS_Snowmass->Write();
   //  out_ttbar->Write();
   // out_ZZ4Mu->Write();
   //  out_DY->Write(); 
   // out_ZZbb2Mu->Write();  
   // out_ZZZbb4Mu->Write();
   // out_ZWpbbMuNL->Write(); 
    /* ZWp2MuMuNL->Write(); 
    ZWmbbMuNL->Write();  
    ZWm2MuMuNL->Write(); 
    WW2Mu2NL->Write(); 
    smHZZ4Mu->Write();  */
   //out_sm_HTobb->Write();
   // out_sm_H4Mu->Write();
   
   cout << "saving..." << endl;    
   out_sig->Close();   

   //out_TT_CMS_Snowmass->Close();
   // out_ttbar->Close();
   //  out_ZZ4Mu->Close();
   //  out_DY->Close(); 
   // out_ZZbb2Mu->Close();  
   // out_ZZZbb4Mu->Close();
   // out_ZWpbbMuNL->Close(); 
    /* ZWp2MuMuNL->Close(); 
    ZWmbbMuNL->Close();  
    ZWm2MuMuNL->Close(); 
    WW2Mu2NL->Close(); 
    smHZZ4Mu->Close();  */
   // out_sm_HTobb->Close();
   // out_sm_H4Mu->Close();
   
   cout << "---DONE---" << endl;
   // cout << "ROOT file: " << out_sig << " has been created sucessfully!" << endl;
   
   
}
