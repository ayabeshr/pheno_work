#define HHBBWWjj2l_cxx
#include "HHBBWWjj2l.h"
#include <TROOT.h>
#include "TH1.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include "THStack.h"
#include "TFile.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>

void HHBBWWjj2l::Loop()
{

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   
   TFile *inputfile = new TFile("/home/aya/Delphes/Delphes-3.4.2/HHBBWW_qq~l-vl~_2HDM_h2-400GeV_13TeV_10000events_SIM.root", "READ");
   //TFile *outputfile = new TFile("/home/aya/ROOTCERN/buildroot/macros/Output.root", "RECREATE" );
  // TFile *outputfile = new TFile("/home/aya/ROOTCERN/buildroot/macros/Histograms after cut.root", "RECREATE" );
   TFile *outputfile = new TFile("/home/aya/ROOTCERN/buildroot/macros/Reconstructed Histograms.root", "RECREATE");
   
   float pts[10]={0., 25., 50. , 75. ,100. ,125. ,150. ,175. ,200. ,225.};
   float Ms[9]={0., 100., 200., 300., 400., 500., 600., 700., 800.};
   int nbins=15 ,npt=9 ,nMs=8;
   
   float Ldata= 3e+5;
   float Lmc=(10000/519.7);
   float wt=Ldata/Lmc ;
   
   TCanvas *c = new TCanvas("c","c",600,400);
   gStyle->SetOptStat("emr");                                      
                                        //Histograms for MET before applying cuts//
                                            
                                            
   TH1D *MET = new TH1D("MET", "MET Distribution", 20, 0, 200);
   TH1D *MET_Eta = new TH1D("MET_Eta", "Eta Distribution for MET", 20, -10, 10);
   TH1D *MET_Phi = new TH1D("MET_Phi", "Phi Distribution for MET", 20, -10, 10);
   
   MET->GetXaxis()->SetTitle("MET (GeV)");
   MET->GetYaxis()->SetTitle("Events");
   MET->SetMarkerStyle(20);
   MET_Eta->GetXaxis()->SetTitle("Eta");
   MET_Eta->GetYaxis()->SetTitle("Events");
   MET_Eta->SetMarkerStyle(20);
   MET_Phi->GetXaxis()->SetTitle("Phi");
   MET_Phi->GetYaxis()->SetTitle("Events");
   MET_Phi->SetMarkerStyle(20);
   
   
                                      //Histograms for MET after applying cuts//
                                   
   TH1D *MET_afterCut = new TH1D("MET_afterCut", "MET Distribution after applying cut", 20, 0, 200);
   TH1D *MET_Eta_afterCut = new TH1D("MET_Eta_afterCut", "Eta Distribution for MET after applying cut", 20, -10, 10);
   TH1D *MET_Phi_afterCut = new TH1D("MET_Phi_afterCut", "Phi Distribution for MET after applying cut", 20, -10, 10);
   
   MET_afterCut->GetXaxis()->SetTitle("MET (GeV)");
   MET_afterCut->GetYaxis()->SetTitle("Events");
   MET_afterCut->SetMarkerStyle(20);
   MET_Eta_afterCut->GetXaxis()->SetTitle("MET {Eta}");
   MET_Eta_afterCut->GetYaxis()->SetTitle("Events");
   MET_Eta_afterCut->SetMarkerStyle(20);
   MET_Phi_afterCut->GetXaxis()->SetTitle("MET {Phi}");
   MET_Phi_afterCut->GetYaxis()->SetTitle("Events");
   MET_Phi_afterCut->SetMarkerStyle(20);
                               
                          
  
                                   //Histograms for Electrons before applying cuts//
   
   
    TH1D *Electron_Pt = new TH1D("Electron_Pt", "Pt Distribution for all Electrons", 20, 0, 100);
    TH1D *Electron_ETA = new TH1D("Electron_ETA", "Eta Distribution for all Electrons", 20, -7, 7);
    TH1D *Electron_pt = new TH1D("Electron_pt", "Pt Distribution for Leading Electron", 20, 0, 100);
    TH1D *Electron_eta = new TH1D("Electron_eta", "Eta Distribution for Leading Electron", 20, -7, 7);
    
    Electron_Pt->GetXaxis()->SetTitle("pT (GeV/c)");
    Electron_Pt->GetYaxis()->SetTitle("Events");
    Electron_Pt->SetMarkerStyle(20);
    Electron_ETA->GetXaxis()->SetTitle("Eta");
    Electron_ETA->GetYaxis()->SetTitle("Events");
    Electron_ETA->SetMarkerStyle(20);
    Electron_pt->GetXaxis()->SetTitle("pT (GeV/c)");
    Electron_pt->GetYaxis()->SetTitle("Events");
    Electron_pt->SetMarkerStyle(20);
    Electron_eta->GetXaxis()->SetTitle("Eta");
    Electron_eta->GetYaxis()->SetTitle("Events");
    Electron_eta->SetMarkerStyle(20);
    
    
                                    //Histograms for Electrons after applying cuts//
   
   
    TH1D *Electron_Pt_afterCut = new TH1D("Electron_Pt_afterCut", "Pt Distribution for all Electrons after applying cut", 20, 0, 500);
    TH1D *Electron_ETA_afterCut = new TH1D("Electron_ETA_afterCut", "Eta Distribution for all Electrons after applying cut", 20, -7, 7);
    TH1D *Electron_pt_afterCut = new TH1D("Electron_pt_afterCut", "Pt Distribution for Leading Electron after applying cut", 20, 0, 500);
    TH1D *Electron_eta_afterCut = new TH1D("Electron_eta_afterCut", "Eta Distribution for Leading Electron after applying cut", 20, -7, 7);
    
    Electron_Pt_afterCut->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    Electron_Pt_afterCut->GetYaxis()->SetTitle("Events");
    Electron_Pt_afterCut->SetMarkerStyle(20);
    Electron_ETA_afterCut->GetXaxis()->SetTitle("Electron {Eta}");
    Electron_ETA_afterCut->GetYaxis()->SetTitle("Events");
    Electron_ETA_afterCut->SetMarkerStyle(20);
    Electron_pt_afterCut->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    Electron_pt_afterCut->GetYaxis()->SetTitle("Events");
    Electron_pt_afterCut->SetMarkerStyle(20);
    Electron_eta_afterCut->GetXaxis()->SetTitle("Electron {Eta}");
    Electron_eta_afterCut->GetYaxis()->SetTitle("Events");
    Electron_eta_afterCut->SetMarkerStyle(20); 
    
   
                                     //Histograms for Muons before applying cuts//
   
   
    TH1D *Muon_Pt = new TH1D("Muon_Pt", "Pt Distribution for all Muons", 20, 0, 200);
    TH1D *Muon_ETA = new TH1D("Muon_ETA", "Eta Distribution for all Muons", 20, -7, 7);
    TH1D *Muon_pt = new TH1D("Muon_pt", "Pt Distribution for Leading Muon", 20, 0, 200);
    TH1D *Muon_eta = new TH1D("Muon_eta", "Eta Distribution for Leading Muon", 20, -7, 7);
    
    Muon_Pt->GetXaxis()->SetTitle("pT (GeV/c)");
    Muon_Pt->GetYaxis()->SetTitle("Events");
    Muon_Pt->SetMarkerStyle(20);
    Muon_ETA->GetXaxis()->SetTitle("Eta");
    Muon_ETA->GetYaxis()->SetTitle("Events");
    Muon_ETA->SetMarkerStyle(20);
    Muon_pt->GetXaxis()->SetTitle("pT (GeV/c)");
    Muon_pt->GetYaxis()->SetTitle("Events");
    Muon_pt->SetMarkerStyle(20);
    Muon_eta->GetXaxis()->SetTitle("Eta");
    Muon_eta->GetYaxis()->SetTitle("Events");
    Muon_eta->SetMarkerStyle(20);
    
    
                                      //Histograms for Muons after applying cuts//
   
   
    TH1D *Muon_Pt_afterCut = new TH1D("Muon_Pt_afterCut", "Pt Distribution for all Muons after applying cut", 20, 0, 200);
    TH1D *Muon_ETA_afterCut = new TH1D("Muon_ETA_afterCut", "Eta Distribution for all Muons after applying cut", 20, -7, 7);
    TH1D *Muon_pt_afterCut = new TH1D("Muon_pt_afterCut", "Pt Distribution for Leading Muon after applying cut", 20, 0, 200);
    TH1D *Muon_eta_afterCut = new TH1D("Muon_eta_afterCut", "Eta Distribution for Leading Muon after applying cut", 20, -7, 7);
    
    Muon_Pt_afterCut->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    Muon_Pt_afterCut->GetYaxis()->SetTitle("Events");
    Muon_Pt_afterCut->SetMarkerStyle(20);
    Muon_ETA_afterCut->GetXaxis()->SetTitle("Muon {Eta}");
    Muon_ETA_afterCut->GetYaxis()->SetTitle("Events");
    Muon_ETA_afterCut->SetMarkerStyle(20);
    Muon_pt_afterCut->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    Muon_pt_afterCut->GetYaxis()->SetTitle("Events");
    Muon_pt_afterCut->SetMarkerStyle(20);
    Muon_eta_afterCut->GetXaxis()->SetTitle("Muon {Eta}");
    Muon_eta_afterCut->GetYaxis()->SetTitle("Events");
    Muon_eta_afterCut->SetMarkerStyle(20);
   
   
                                        //Histograms for Jets before applying cuts//
                                              
     
    TH1D *Jet_Pt = new TH1D("Jet_Pt", "Pt Distribution for all Jets", 20, 0, 300);  
    TH1D *Jet_ETA = new TH1D("Jet_ETA", "Eta Distribution for all Jets", 20, -10, 10);
    TH1D *Jet0_Pt = new TH1D("Jet0_Pt", "Pt Distribution for 1st Jet", 20, 0, 300);   
    TH1D *Jet0_Eta = new TH1D("Jet0_Eta", "Eta Distribution for 1st Jet", 20, -7, 7);                                     
    TH1D *Jet1_Pt = new TH1D("Jet1_Pt", "Pt Distribution for 2nd Jet", 20, 0, 300);   
    TH1D *Jet1_Eta = new TH1D("Jet1_Eta", "Eta Distribution for 2nd Jet", 20, -7, 7); 
    TH1D *Jet2_Pt = new TH1D("Jet2_Pt", "Pt Distribution for 3rd Jet", 20, 0, 300);   
    TH1D *Jet2_Eta = new TH1D("Jet2_Eta", "Eta Distribution for 3rd Jet", 20, -7, 7);  
    TH1D *Jet3_Pt = new TH1D("Jet3_Pt", "Pt Distribution for 4th Jet", npt, pts);   
    TH1D *Jet3_Eta = new TH1D("Jet3_Eta", "Eta Distribution for 4th Jet", 20, -7, 7);
    TH1D *BTag0_Pt = new TH1D("BTag0_Pt", "Pt Distribution for 1st assumed BTag", 20, 0, 300);
    TH1D *BTag0_Eta = new TH1D("BTag0_Eta", "Eta Distribution for 1st assumed BTag", 20, -7, 7);    
    TH1D *BTag1_Pt = new TH1D("BTag1_Pt", "Pt Distribution for 2nd assumed BTag", 20, 0, 300);
    TH1D *BTag1_Eta = new TH1D("BTag1_Eta", "Eta Distribution for 2nd assumed BTag", 20, -7, 7);                                                                                                       
    TH1D *BTag0Gen_Pt = new TH1D("BTag0Gen_Pt", "Pt Distribution for 1st assumed generated BTag", npt, pts);
    TH1D *BTag0Gen_Eta= new TH1D("BTag0Gen_Eta", "Eta Distribution for 1st assumed generated BTag", 20, -7, 7);
    TH1D *BTag1Gen_Pt = new TH1D("BTag1Gen_Pt", "Pt Distribution for 2nd assumed generated BTag", npt, pts);
    TH1D *BTag1Gen_Eta= new TH1D("BTag1Gen_Eta", "Eta Distribution for 2nd assumed generated BTag", 20, -7, 7);
   
    Jet_Pt->GetXaxis()->SetTitle("pT (GeV/c)");
    Jet_Pt->GetYaxis()->SetTitle("Events");
    Jet_Pt->SetMarkerStyle(20);
    
    Jet_ETA->GetXaxis()->SetTitle("Eta");
    Jet_ETA->GetYaxis()->SetTitle("Events");
    Jet_ETA->SetMarkerStyle(20);
    
    Jet0_Pt->GetXaxis()->SetTitle("pT (GeV/c)");
    Jet0_Pt->GetYaxis()->SetTitle("Events");
    Jet0_Pt->SetMarkerStyle(20);
    
    Jet0_Eta->GetXaxis()->SetTitle("Eta");
    Jet0_Eta->GetYaxis()->SetTitle("Events");
    Jet0_Eta->SetMarkerStyle(20);
    
    Jet1_Pt->GetXaxis()->SetTitle("pT (GeV/c)");
    Jet1_Pt->GetYaxis()->SetTitle("Events");
    Jet1_Pt->SetMarkerStyle(20);
    //Jet1_Pt->SetLineColor(kRed);
    
    Jet1_Eta->GetXaxis()->SetTitle("Eta");
    Jet1_Eta->GetYaxis()->SetTitle("Events");
    Jet1_Eta->SetMarkerStyle(20);
    
    Jet2_Pt->GetXaxis()->SetTitle("pT (GeV/c)");
    Jet2_Pt->GetYaxis()->SetTitle("Events");
    Jet2_Pt->SetMarkerStyle(20);
    
    Jet2_Eta->GetXaxis()->SetTitle("Eta");
    Jet2_Eta->GetYaxis()->SetTitle("Events");
    Jet2_Eta->SetMarkerStyle(20);
    
    Jet3_Pt->GetXaxis()->SetTitle("pT (GeV/c)");
    Jet3_Pt->GetYaxis()->SetTitle("Events");
    Jet3_Pt->SetMarkerStyle(20);
    
    Jet3_Eta->GetXaxis()->SetTitle("Eta");
    Jet3_Eta->GetYaxis()->SetTitle("Events");
    Jet3_Eta->SetMarkerStyle(20);
    
    BTag0_Pt->GetXaxis()->SetTitle("pT (GeV/c)");
    BTag0_Pt->GetYaxis()->SetTitle("Events");
    BTag0_Pt->SetMarkerStyle(20);
    
    BTag0_Eta->GetXaxis()->SetTitle("Eta");
    BTag0_Eta->GetYaxis()->SetTitle("Events"); 
    BTag0_Eta->SetMarkerStyle(20);
    
    BTag1_Pt->GetXaxis()->SetTitle("pT (GeV/c)");
    BTag1_Pt->GetYaxis()->SetTitle("Events");
    BTag1_Pt->SetMarkerStyle(20);
    
    BTag1_Eta->GetXaxis()->SetTitle("Eta");
    BTag1_Eta->GetYaxis()->SetTitle("Events"); 
    BTag1_Eta->SetMarkerStyle(20);
    
    BTag0Gen_Pt->GetXaxis()->SetTitle("pT (GeV/c)");
    BTag0Gen_Pt->GetYaxis()->SetTitle("Events");
    BTag0Gen_Pt->SetMarkerStyle(20);
    
    BTag0Gen_Eta->GetXaxis()->SetTitle("Eta");
    BTag0Gen_Eta->GetYaxis()->SetTitle("Events");
    BTag0Gen_Eta->SetMarkerStyle(20);
     
    BTag1Gen_Pt->GetXaxis()->SetTitle("pT (GeV/c)");
    BTag1Gen_Pt->GetYaxis()->SetTitle("Events");
    BTag1Gen_Pt->SetMarkerStyle(20);
    
    BTag1Gen_Eta->GetXaxis()->SetTitle("Eta");
    BTag1Gen_Eta->GetYaxis()->SetTitle("Events"); 
    BTag1Gen_Eta->SetMarkerStyle(20);
   
   
                                      //Histograms for Jets after applying cuts//
                                              
     
    TH1D *Jet_Pt_afterCut = new TH1D("Jet_Pt_afterCut", "Pt Distribution for all Jets after applying cut", 20, 0, 300);  
    TH1D *Jet_ETA_afterCut = new TH1D("Jet_ETA_afterCut", "Eta Distribution for all Jets after applying cut", 20, -10, 10);
    TH1D *Jet0_Pt_afterCut = new TH1D("Jet0_Pt_afterCut", "Pt Distribution for 1st Jet after applying cut", 20, 0, 300);   
    TH1D *Jet0_Eta_afterCut = new TH1D("Jet0_Eta_afterCut", "Eta Distribution for 1st Jet after applying cut", 20, -7, 7);                                     
    TH1D *Jet1_Pt_afterCut = new TH1D("Jet1_Pt_afterCut", "Pt Distribution for 2nd Jet after applying cut", 20, 0, 300);   
    TH1D *Jet1_Eta_afterCut = new TH1D("Jet1_Eta_afterCut", "Eta Distribution for 2nd Jet after applying cut", 20, -7, 7); 
    TH1D *Jet2_Pt_afterCut = new TH1D("Jet2_Pt_afterCut", "Pt Distribution for 3rd Jet after applying cut", 20, 0, 300);   
    TH1D *Jet2_Eta_afterCut = new TH1D("Jet2_Eta_afterCut", "Eta Distribution for 3rd Jet after applying cut", 20, -7, 7);  
    TH1D *Jet3_Pt_afterCut = new TH1D("Jet3_Pt_afterCut", "Pt Distribution for 4th Jet after applying cut", npt, pts);   
    TH1D *Jet3_Eta_afterCut = new TH1D("Jet3_Eta_afterCut", "Eta Distribution for 4th Jet after applying cut", 20, -7, 7);
    TH1D *BTag0_Pt_afterCut = new TH1D("BTag0_Pt_afterCut", "Pt Distribution for 1st assumed BTag after applying cut", 20, 0, 300);
    TH1D *BTag0_Eta_afterCut = new TH1D("BTag0_Eta_afterCut", "Eta Distribution for 1st assumed BTag after applying cut", 20, -7, 7);    
    TH1D *BTag1_Pt_afterCut = new TH1D("BTag1_Pt_afterCut", "Pt Distribution for 2nd assumed BTag after applying cut", 20, 0, 300);
    TH1D *BTag1_Eta_afterCut = new TH1D("BTag1_Eta_afterCut", "Eta Distribution for 2nd assumed BTag after applying cut", 20, -7, 7); 
    
    
    Jet_Pt_afterCut->GetXaxis()->SetTitle("Jet p_{T} (GeV/c)");
    Jet_Pt_afterCut->GetYaxis()->SetTitle("Events");
    Jet_Pt_afterCut->SetMarkerStyle(20);
    
    Jet_ETA_afterCut->GetXaxis()->SetTitle("Jet {Eta}");
    Jet_ETA_afterCut->GetYaxis()->SetTitle("Events");
    Jet_ETA_afterCut->SetMarkerStyle(20);
    
    Jet0_Pt_afterCut->GetXaxis()->SetTitle("Jet0 p_{T} (GeV/c)");
    Jet0_Pt_afterCut->GetYaxis()->SetTitle("Events");
    Jet0_Pt_afterCut->SetMarkerStyle(20);
    
    Jet0_Eta_afterCut->GetXaxis()->SetTitle("Jet0 {Eta}");
    Jet0_Eta_afterCut->GetYaxis()->SetTitle("Events");
    Jet0_Eta_afterCut->SetMarkerStyle(20);
    
    Jet1_Pt_afterCut->GetXaxis()->SetTitle("Jet1 p_{T} (GeV/c)");
    Jet1_Pt_afterCut->GetYaxis()->SetTitle("Events");
    Jet1_Pt_afterCut->SetMarkerStyle(20);
    //Jet1_Pt->SetLineColor(kRed);
    
    Jet1_Eta_afterCut->GetXaxis()->SetTitle("Jet1 {Eta}");
    Jet1_Eta_afterCut->GetYaxis()->SetTitle("Events");
    Jet1_Eta_afterCut->SetMarkerStyle(20);
    
    Jet2_Pt_afterCut->GetXaxis()->SetTitle("Jet2 p_{T} (GeV/c)");
    Jet2_Pt_afterCut->GetYaxis()->SetTitle("Events");
    Jet2_Pt_afterCut->SetMarkerStyle(20);
    
    Jet2_Eta_afterCut->GetXaxis()->SetTitle("Jet2 {Eta}");
    Jet2_Eta_afterCut->GetYaxis()->SetTitle("Events");
    Jet2_Eta_afterCut->SetMarkerStyle(20);
    
    Jet3_Pt_afterCut->GetXaxis()->SetTitle("Jet3 p_{T} (GeV/c)");
    Jet3_Pt_afterCut->GetYaxis()->SetTitle("Events");
    Jet3_Pt_afterCut->SetMarkerStyle(20);
    
    Jet3_Eta_afterCut->GetXaxis()->SetTitle("Jet3 {Eta}");
    Jet3_Eta_afterCut->GetYaxis()->SetTitle("Events");
    Jet3_Eta_afterCut->SetMarkerStyle(20);
    
    BTag0_Pt_afterCut->GetXaxis()->SetTitle("BTag0 p_{T} (GeV/c)");
    BTag0_Pt_afterCut->GetYaxis()->SetTitle("Events");
    BTag0_Pt_afterCut->SetMarkerStyle(20);
    
    BTag0_Eta_afterCut->GetXaxis()->SetTitle("BTag0 {Eta}");
    BTag0_Eta_afterCut->GetYaxis()->SetTitle("Events"); 
    BTag0_Eta_afterCut->SetMarkerStyle(20);
    
    BTag1_Pt_afterCut->GetXaxis()->SetTitle("BTag1 p_{T} (GeV/c)");
    BTag1_Pt_afterCut->GetYaxis()->SetTitle("Events");
    BTag1_Pt_afterCut->SetMarkerStyle(20);
    
    BTag1_Eta_afterCut->GetXaxis()->SetTitle("BTag1 {Eta}");
    BTag1_Eta_afterCut->GetYaxis()->SetTitle("Events"); 
    BTag1_Eta_afterCut->SetMarkerStyle(20);
   
                             
                         //Histograms for reconstructing 1st Higgs from b-bbarjets with cuts applied//   
                          
   
   TH1D *Higgs1_mass = new TH1D("1st_Higgs_mass", "Invariant mass distribution for 1st Higgs", 20, 50, 170);
   TH1D *Higgs1_pt = new TH1D("1st_Higgs_pt", "Pt distribution for 1st Higgs", npt, pts);
   TH1D *Higgs1_eta = new TH1D("1st_Higgs_eta", "Eta distribution for 1st Higgs", 20, -7, 7);
   
   Higgs1_mass->GetXaxis()->SetTitle("Higgs Invariant Mass (GeV)");
   Higgs1_mass->GetXaxis()->CenterTitle();
   Higgs1_mass->GetYaxis()->SetTitle("Events"); 
   Higgs1_mass->GetYaxis()->CenterTitle();
   Higgs1_mass->SetMarkerStyle(20);
   Higgs1_pt->GetXaxis()->SetTitle("Higgs p_{T} (GeV/c)");
   Higgs1_pt->GetXaxis()->CenterTitle();
   Higgs1_pt->GetYaxis()->SetTitle("Events");
   Higgs1_pt->GetYaxis()->CenterTitle();
   Higgs1_pt->SetMarkerStyle(20);
   Higgs1_eta->GetXaxis()->SetTitle("Higgs Eta");
   Higgs1_eta->GetXaxis()->CenterTitle();
   Higgs1_eta->GetYaxis()->SetTitle("Events");
   Higgs1_eta->GetYaxis()->CenterTitle();
   Higgs1_eta->SetMarkerStyle(20);
 
 
                                 //Histograms for reconstructing W+ boson from dijets//
                                 
   TH1D *Wplus_pt = new TH1D("Wplus_pt", "Pt distribution for W+", npt, pts);
   TH1D *Wplus_eta = new TH1D("Wplus_eta", "Eta distribution for W+", 20, -7, 7);
   TH1D *Wplus_mass = new TH1D("Wplus_mass", "Invariant Mass distribution for W+", 20, 30, 150);
   
   Wplus_pt->GetXaxis()->SetTitle("W+ p_{T} (GeV/c)");
   Wplus_pt->GetXaxis()->CenterTitle();
   Wplus_pt->GetYaxis()->SetTitle("Events");
   Wplus_pt->GetYaxis()->CenterTitle();
   Wplus_pt->SetMarkerStyle(20);
   
   Wplus_eta->GetXaxis()->SetTitle("W+ Eta");
   Wplus_eta->GetXaxis()->CenterTitle();
   Wplus_eta->GetYaxis()->SetTitle("Events");
   Wplus_eta->GetYaxis()->CenterTitle();
   Wplus_eta->SetMarkerStyle(20);
                                       
   Wplus_mass->GetXaxis()->SetTitle("W+ Invariant Mass (GeV)");
   Wplus_mass->GetXaxis()->CenterTitle();
   Wplus_mass->GetYaxis()->SetTitle("Events");
   Wplus_mass->GetYaxis()->CenterTitle();
   Wplus_mass->SetMarkerStyle(20);
   
                                       //***********************************************//
                                       //Histograms for reconstructing W- from dileptons//
                                       //***********************************************//
                                       
   TH1F *Wm_mass = new TH1F("Wm_mass", "Transverse mass distribution for dileptons", 50, 0, 200);                                   
   Wm_mass->GetXaxis()->SetTitle("W- {mass} - (GeV/c^2)");
   Wm_mass->GetYaxis()->SetTitle("Events");
   
   TH1F *Wm_eta = new TH1F("Wm_eta", "Pseudorapidity distribution for dileptons", 20, -7., 7.);
   Wm_eta->GetXaxis()->SetTitle("W- {eta}");
   Wm_eta->GetYaxis()->SetTitle("Events");
   
   TH1F *Wm_pT = new TH1F("Wm_pT", "Transverse mom. distribution for dileptons", npt, pts);
   Wm_pT->GetXaxis()->SetTitle("W- {pT} - (GeV/c)");
   Wm_pT->GetYaxis()->SetTitle("Events");
                                     
                                       
                                    //**********************************************//
                                   //Histograms for reconstructing 2nd Higgs from Ws//   
                                   //************************************************//
   
   TH1F *Higgs2_mass = new TH1F("Higgs2_mass", "Invariant mass distribution for W bosons", nMs, Ms);
   Higgs2_mass->GetXaxis()->SetTitle("b bbar invariant mass - (GeV)");
   Higgs2_mass->GetYaxis()->SetTitle("Events");
   
   TH1F *Higgs2_eta = new TH1F("Higgs2_eta", "Pseudorapidity distribution for W bosons", 20, -7., 7.);
   Higgs2_eta->GetXaxis()->SetTitle("eta");
   Higgs2_eta->GetYaxis()->SetTitle("Events");
   
   TH1F *Higgs2_pT = new TH1F("Higgs2_pT", "Transverse mom. disrtibution for W bosons", npt, pts);
   Higgs2_pT->GetXaxis()->SetTitle("pT (GeV/c)");
   Higgs2_pT->GetYaxis()->SetTitle("Events");
   
   
   Long64_t nbytes = 0, nb = 0;
   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      
    /*  for (Int_t i = 0; i < MissingET_size; i++)
   
   {
	   MET->Fill(MissingET_MET[i], wt);
	   
	   MET_Eta->Fill(MissingET_Eta[i], wt);
	   
	   MET_Phi->Fill(MissingET_Phi[i], wt); 
	    
    }
    
    
   
   for (Int_t j = 0; j < Electron_size; j++)
    {
		Electron_Pt->Fill(Electron_PT[j], wt);
		
		Electron_ETA->Fill(Electron_Eta[j], wt);
		
	}
   
    Electron_pt->Fill(Electron_PT[0], wt);
    Electron_eta->Fill(Electron_Eta[0], wt); 
   
   
   
   for (Int_t k = 0; k < Muon_size; k++)
    {
		Muon_Pt->Fill(Muon_PT[k], wt);
		
		Muon_ETA->Fill(Muon_Eta[k], wt);
		
	}
   
    Muon_pt->Fill(Muon_PT[0], wt);
    Muon_eta->Fill(Muon_Eta[0], wt); 
    
   
   
   for (Int_t l = 0; l < Jet_size; l++)
   {
	   Jet_Pt->Fill(Jet_PT[l], wt);
	   
	   Jet_ETA->Fill(Jet_Eta[l], wt);
	      
   }
   
   
   Jet0_Pt->Fill(Jet_PT[0], wt);
   Jet0_Eta->Fill(Jet_Eta[0], wt);
   Jet1_Pt->Fill(Jet_PT[1], wt);
   Jet1_Eta->Fill(Jet_Eta[1], wt);
   Jet2_Pt->Fill(Jet_PT[2], wt);
   Jet2_Eta->Fill(Jet_Eta[2], wt);
   Jet3_Pt->Fill(Jet_PT[3], wt);
   Jet3_Eta->Fill(Jet_Eta[3], wt);
   
   
   if (Jet_BTag[0] == 1) 
   {
	  BTag0_Pt->Fill(Jet_PT[0], wt);
	 
	  BTag0_Eta->Fill(Jet_Eta[0], wt);
	  
   }
   
   else 
   {
	  BTag0_Pt->Fill(GenJet_PT[0], wt);
	 
      BTag0_Eta->Fill(GenJet_Eta[0], wt);
	    
   }
   
   
   if (Jet_BTag[1] == 1) 
   {
	  BTag1_Pt->Fill(Jet_PT[1], wt);
	 
	  BTag1_Eta->Fill(Jet_Eta[1], wt);
	  
   }
   
   else 
   {
	  BTag1_Pt->Fill(GenJet_PT[1], wt);
	 
      BTag1_Eta->Fill(GenJet_Eta[1], wt);
	    
   }
   
   
   if (GenJet_BTag[0] == 1)
   {
	   BTag0Gen_Pt->Fill(GenJet_PT[0], wt);
	   
	   BTag0Gen_Eta->Fill(GenJet_Eta[0], wt);
	   
   }
   
   
   if (GenJet_BTag[1] == 1)
   {
	   BTag1Gen_Pt->Fill(GenJet_PT[1], wt);
	   
	   BTag1Gen_Eta->Fill(GenJet_Eta[1], wt);
	   
   }*/
      
   
                                                  //Applying Cuts//
                                                 
  ///////////////////////////MET CUTS////////////////////////////////////
   
  /* for (Int_t n = 0; n < MissingET_size; n++)
   {
       if (MissingET_MET[n] > 20)
       {
         MET_afterCut->Fill(MissingET_MET[n], wt);
	   
	     MET_Eta_afterCut->Fill(MissingET_Eta[n], wt);
	   
	     MET_Phi_afterCut->Fill(MissingET_Phi[n], wt); 
	       
       }
   }
   
  //////////////////////////ELECTRON CUTS////////////////////////////////
   
   for (Int_t p = 0; p < Electron_size; p++)
   {
	   if ((Electron_PT[p] > 27) && (fabs(Electron_Eta[p]) < 2.4))   
	   {
	      Electron_Pt_afterCut->Fill(Electron_PT[p], wt);
		
		  Electron_ETA_afterCut->Fill(Electron_Eta[p], wt);
		 
       }
   }
   
   
   if ((Electron_PT[0] > 27) && (fabs(Electron_Eta[0]) < 2.4))
   {
	   Electron_pt_afterCut->Fill(Electron_PT[0], wt);
    
       Electron_eta_afterCut->Fill(Electron_Eta[0], wt); 
	 
   }   
   
  //////////////////////////////MUON CUTS////////////////////////////////
   
   for (Int_t q = 0; q < Muon_size; q++)
   { 
		if ((Muon_PT[q] > 24) && (fabs(Muon_Eta[q]) < 2.4))
		{
			Muon_Pt_afterCut->Fill(Muon_PT[q], wt);
		
		    Muon_ETA_afterCut->Fill(Muon_Eta[q], wt);
		}
	}
   
   if ((Muon_PT[0] > 24) && (fabs(Muon_Eta[0]) < 2.4))
   { 
       Muon_pt_afterCut->Fill(Muon_PT[0], wt);
      
       Muon_eta_afterCut->Fill(Muon_Eta[0], wt); 
    
   }
   
  ///////////////////////////////JET CUTS///////////////////////////////
   
   for (Int_t r = 0; r < Jet_size; r++)
   {
	   if ((Jet_PT[r] > 30) && (fabs(Jet_Eta[r]) < 2.5))
	   {
		   Jet_Pt_afterCut->Fill(Jet_PT[r], wt);
	   
	       Jet_ETA_afterCut->Fill(Jet_Eta[r], wt);   
	   }
	   
   }
   
   
    if ((Jet_PT[0] > 30) && (fabs(Jet_Eta[0]) < 2.5))
    {
		Jet0_Pt_afterCut->Fill(Jet_PT[0], wt);
   
        Jet0_Eta_afterCut->Fill(Jet_Eta[0], wt);
		
	}
    
    
    if ((Jet_PT[1] > 30) && (fabs(Jet_Eta[1]) < 2.5))
    {
		Jet1_Pt_afterCut->Fill(Jet_PT[1], wt);
   
        Jet1_Eta_afterCut->Fill(Jet_Eta[1], wt);
		
	}
    
    
    if ((Jet_PT[2] > 30) && (fabs(Jet_Eta[2]) < 2.5))
    {
		Jet2_Pt_afterCut->Fill(Jet_PT[2], wt);
   
        Jet2_Eta_afterCut->Fill(Jet_Eta[2], wt);
		
	}
    
    
    if ((Jet_PT[3] > 30) && (fabs(Jet_Eta[3]) < 2.5))
    {
		Jet3_Pt_afterCut->Fill(Jet_PT[3], wt);
   
        Jet3_Eta_afterCut->Fill(Jet_Eta[3], wt);
		
	}
    
    
    if (Jet_BTag[0] == 1) 
   {
	  if ((Jet_PT[0] > 30) && (fabs(Jet_Eta[0]) < 2.5))
	  {
	     BTag0_Pt_afterCut->Fill(Jet_PT[0], wt);
	 
	     BTag0_Eta_afterCut->Fill(Jet_Eta[0], wt);
	  }
	  
   }
    
    
    if (Jet_BTag[1] == 1) 
   {
	  if ((Jet_PT[1] > 30) && (fabs(Jet_Eta[1]) < 2.5))
	  {
	     BTag1_Pt_afterCut->Fill(Jet_PT[1], wt);
	 
	     BTag1_Eta_afterCut->Fill(Jet_Eta[1], wt);
	  }
	  
   }*/
    
    
    if ((Jet_PT[0] > 30) && (fabs(Jet_Eta[0]) < 2.5)) continue;
    if ((Jet_PT[1] > 30) && (fabs(Jet_Eta[1]) < 2.5)) continue;
    if ((Jet_PT[2] > 30) && (fabs(Jet_Eta[2]) < 2.5)) continue;
    if ((Jet_PT[3] > 30) && (fabs(Jet_Eta[3]) < 2.5)) continue;   
      
     
     //___________________________//
    //TLorentz Vector decleration//
     //___________________________//
     
     TLorentzVector B1, B2, J1, J2, diBjet, diJet;
    
     float mu_mass, b_mass;
    
     
     mu_mass = 0.1057;
     b_mass = 4.18;
  
     
     B1.SetPtEtaPhiM(Jet_PT[0], Jet_Eta[0], Jet_Phi[0], b_mass);
     B2.SetPtEtaPhiM(Jet_PT[1], Jet_Eta[1], Jet_Phi[1], b_mass);
     
     diBjet = B1+B2;
     
     J1.SetPtEtaPhiM(Jet_PT[2], Jet_Eta[2], Jet_Phi[2], Jet_Mass[2]);
     J2.SetPtEtaPhiM(Jet_PT[3], Jet_Eta[3], Jet_Phi[3], Jet_Mass[3]);
     
     diJet = J1 + J2;
     
     //mu.SetPtEtaPhiE(Muon_PT[0], Muon_Eta[0], Muon_Phi[0], mu_mass);
     
    
     
     
     
    // double dR1, dR2;
     
     
    // dR1 = TMath::Sqrt(TMath::Power(( bjet_eta[0] -  bjet_eta[1]),2) + TMath::Power(( bjet_phi[0] -  bjet_phi[1]),2));
    // dR2 = TMath::Sqrt(TMath::Power((jet_eta[0] - jeta_eta[1]),2) + TMath::Power((jet_phi[0] - jeta_phi[1]),2));
     
    // float Higgs1_Pt , Higgs1_Eta, Higgs1_Mass, Wplus_Pt, Wplus_Eta, Wplus_Mass;
     
    float Higgs1_Pt = diBjet.Pt();            
    float Higgs1_Eta = diBjet.Eta(); 
    float Higgs1_Mass = diBjet.M();             
                
    float Wplus_Pt =  diJet.Pt();          
    float Wplus_Eta = diJet.Eta();
    float Wplus_Mass = diJet.M();
                     
     /*Wmin_mass = (mu+ neutrino).Mt();
     Wmin_pt = (mu + neutrino).Pt();
     Wmin_eta = (mu + neutrino).Eta();
          
     higgs2_mass = (Wplus_mass + Wmin_mass).M();
     higgs2_pt = (Wplus_pt + Wmin_pt).Pt();
     higgs2_eta = (Wplus_eta + Wmin_eta).Eta();*/
   
   
                                                    //Filling Histograms//
   
   Higgs1_mass->Fill(Higgs1_Mass, wt);
   Higgs1_pt->Fill(Higgs1_Pt, wt);
   Higgs1_eta->Fill(Higgs1_Eta, wt);
   Wplus_pt->Fill(Wplus_Pt, wt);
   Wplus_eta->Fill(Wplus_Eta, wt);
   Wplus_mass->Fill(Wplus_Mass, wt);
   
   
   
   
   }
   
                                               

                                                 //Writing Histograms//
   /*
   // MET->Write();
	//MET_Eta->Write();
   // MET_Phi->Write();
    MET_afterCut->Write();
	MET_Eta_afterCut->Write();
    MET_Phi_afterCut->Write();  
   // Electron_Pt->Write();
	//Electron_ETA->Write();
    //Electron_pt->Write();
    //Electron_eta->Write(); 
    Electron_Pt_afterCut->Write();
	Electron_ETA_afterCut->Write();
    Electron_pt_afterCut->Write();
    Electron_eta_afterCut->Write(); 
   // Muon_Pt->Write();
	//Muon_ETA->Write();
	//Muon_pt->Write();
   // Muon_eta->Write(); 
    Muon_Pt_afterCut->Write();
	Muon_ETA_afterCut->Write();
	Muon_pt_afterCut->Write();
    Muon_eta_afterCut->Write(); 
   // Jet_Pt->Write();
	//Jet_ETA->Write();
	//Jet0_Pt->Write();
    //Jet0_Eta->Write();
    //Jet1_Pt->Write();
   // Jet1_Eta->Write();
   // Jet2_Pt->Write();
   // Jet2_Eta->Write();
   // Jet3_Pt->Write();
   // Jet3_Eta->Write();
   // BTag0_Pt->Write();
   // BTag0_Eta->Write();
	//BTag1_Pt->Write();
	//BTag1_Eta->Write();
	Jet_Pt_afterCut->Write();
	Jet_ETA_afterCut->Write();
	Jet0_Pt_afterCut->Write();
    Jet0_Eta_afterCut->Write();
    Jet1_Pt_afterCut->Write();
    Jet1_Eta_afterCut->Write();
    Jet2_Pt_afterCut->Write();
    Jet2_Eta_afterCut->Write();
    Jet3_Pt_afterCut->Write();
    Jet3_Eta_afterCut->Write();
    BTag0_Pt_afterCut->Write();
    BTag0_Eta_afterCut->Write();
	BTag1_Pt_afterCut->Write();
	BTag1_Eta_afterCut->Write();
    BTag0Gen_Pt->Write();
	BTag0Gen_Eta->Write();
    BTag1Gen_Pt->Write();
	BTag1Gen_Eta->Write();*/
	Higgs1_mass->Write();
    Higgs1_pt->Write();
    Higgs1_eta->Write();
    Wplus_pt->Write();
    Wplus_eta->Write();
    Wplus_mass->Write();
   
	
	
	                                          //Drawing Histograms//
	/*   
	//MET->Draw("COLZ");
	//MET_Eta->Draw("COLZ");
    //MET_Phi->Draw("COLZ"); 
    MET_afterCut->Draw("COLZ");
	MET_Eta_afterCut->Draw("COLZ");
    MET_Phi_afterCut->Draw("COLZ"); 
   // Electron_Pt->Draw("COLZ");
	//Electron_ETA->Draw("COLZ");
   // Electron_pt->Draw("COLZ");
    //Electron_eta->Draw("COLZ"); 
    Electron_Pt_afterCut->Draw("COLZ");
	Electron_ETA_afterCut->Draw("COLZ");
    Electron_pt_afterCut->Draw("COLZ");
    Electron_eta_afterCut->Draw("COLZ"); 
   // Muon_Pt->Draw("COLZ");
	//Muon_ETA->Draw("COLZ");
	//Muon_pt->Draw("COLZ");
   // Muon_eta->Draw("COLZ"); 
    Muon_Pt_afterCut->Draw("COLZ");
	Muon_ETA_afterCut->Draw("COLZ");
	Muon_pt_afterCut->Draw("COLZ");
    Muon_eta_afterCut->Draw("COLZ");
    //Jet_Pt->Draw("COLZ");
	//Jet_ETA->Draw("COLZ");
	//Jet0_Pt->Draw("COLZ");
   // Jet0_Eta->Draw("COLZ");
   // Jet1_Pt->Draw("COLZ");
   // Jet1_Eta->Draw("COLZ");
   // Jet2_Pt->Draw("COLZ");
   // Jet2_Eta->Draw("COLZ");
   // Jet3_Pt->Draw("COLZ");
   // Jet3_Eta->Draw("COLZ");
   // BTag0_Pt->Draw("COLZ");
   // BTag0_Eta->Draw("COLZ");
	//BTag1_Pt->Draw("COLZ");
	//BTag1_Eta->Draw("COLZ");
	Jet_Pt_afterCut->Draw("COLZ");
	Jet_ETA_afterCut->Draw("COLZ");
	Jet0_Pt_afterCut->Draw("COLZ");
    Jet0_Eta_afterCut->Draw("COLZ");
    Jet1_Pt_afterCut->Draw("COLZ");
    Jet1_Eta_afterCut->Draw("COLZ");
    Jet2_Pt_afterCut->Draw("COLZ");
    Jet2_Eta_afterCut->Draw("COLZ");
    Jet3_Pt_afterCut->Draw("COLZ");
    Jet3_Eta_afterCut->Draw("COLZ");
    BTag0_Pt_afterCut->Draw("COLZ");
    BTag0_Eta_afterCut->Draw("COLZ");
	BTag1_Pt_afterCut->Draw("COLZ");
	BTag1_Eta_afterCut->Draw("COLZ");
    BTag0Gen_Pt->Draw("COLZ");
	BTag0Gen_Eta->Draw("COLZ");
    BTag1Gen_Pt->Draw("COLZ");
	BTag1Gen_Eta->Draw("COLZ");*/
	Higgs1_mass->Draw("COLZ");
    Higgs1_pt->Draw("COLZ");
    Higgs1_eta->Draw("COLZ");
    Wplus_pt->Draw("COLZ");
    Wplus_eta->Draw("COLZ");
    Wplus_mass->Draw("COLZ");
	
	outputfile->Print();      
	outputfile->Close();   
	    

}
