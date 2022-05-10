//============================================
//  Header file defining Tree and branches
//           for Hhhbb4M
//============================================

#ifndef newtree_h
#define newtree_h

#include <TROOT.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>


struct GeneratedMuons{
	
	
   vector<Float_t> gen_muon_pt;
   vector<Float_t> gen_muon_eta;
   vector<Float_t> gen_muon_phi;
   

};

GeneratedMuons genMuons;

/////////////////////////////////////////////////////

struct GeneratedBHadrons{
	
   vector<Float_t> gen_bjet_pt;
   vector<Float_t> gen_bjet_eta;
   vector<Float_t> gen_bjet_phi;
	
	
};
GeneratedBHadrons gen_bjet;

//////////////////////////////////////////////////////

struct LooseMuons{
	
   vector<Float_t> lmu_pt;
   vector<Float_t> lmu_eta;
   vector<Float_t> lmu_phi;
   vector<Int_t> lmu_size;
   vector<Float_t> muon1_pt;
   vector<Float_t> muon2_pt;
   vector<Float_t> muon3_pt;
   vector<Float_t> muon4_pt;
   vector<Float_t> muon1_eta;	
   vector<Float_t> muon2_eta;
   vector<Float_t> muon3_eta;	
   vector<Float_t> muon4_eta;	
   vector<Float_t> muon1_phi;
   vector<Float_t> muon2_phi;
   vector<Float_t> muon3_phi;
   vector<Float_t> muon4_phi;
   vector<Float_t> dxy_mu1;
   vector<Float_t> dxy_mu2;
   vector<Float_t> dxy_mu3;
   vector<Float_t> dxy_mu4;
   vector<Double_t> loose_mu12_mass;
   vector<Double_t> loose_mu34_mass;
   vector<Double_t> loose_mu13_mass;
   vector<Double_t> loose_mu24_mass;	
   vector<Double_t> loose_mu14_mass;
   vector<Double_t> loose_mu23_mass;
  
};
LooseMuons loose_mu;

//////////////////////////////////////////////////////

struct DeltaR_MuonsLoose{
	
   vector<Float_t> v_DR_mu1mu2;
   vector<Float_t> v_DR_mu3mu4;
   vector<Float_t> v_DR_mu1mu3;
   vector<Float_t> v_DR_mu2mu4;
   vector<Float_t> v_DR_mu1mu4;	
   vector<Float_t> v_DR_mu2mu3;	
	
};
DeltaR_MuonsLoose drMuonsLoose;

//////////////////////////////////////////////////////


struct TightMuons{
	
   vector<Float_t> tmu_pt;
   vector<Float_t> tmu_eta;
   vector<Float_t> tmu_phi;
   vector<Int_t> tmu_size;
  /* vector<Float_t> muon1_pt;
   vector<Float_t> muon2_pt;
   vector<Float_t> muon3_pt;
   vector<Float_t> muon4_pt;
   vector<Float_t> muon1_eta;	
   vector<Float_t> muon2_eta;
   vector<Float_t> muon3_eta;	
   vector<Float_t> muon4_eta;	
   vector<Float_t> muon1_phi;
   vector<Float_t> muon2_phi;
   vector<Float_t> muon3_phi;
   vector<Float_t> muon4_phi;
   vector<Float_t> dxy_mu1;
   vector<Float_t> dxy_mu2;
   vector<Float_t> dxy_mu3;
   vector<Float_t> dxy_mu4;
   vector<Float_t> dz_mu1;
   vector<Float_t> dz_mu2;
   vector<Float_t> dz_mu3;
   vector<Float_t> dz_mu4;
   vector<Double_t> tight_mu12_mass;
   vector<Double_t> tight_mu34_mass;
   vector<Double_t> tight_mu13_mass;
   vector<Double_t> tight_mu24_mass;	
   vector<Double_t> tight_mu14_mass;
   vector<Double_t> tight_mu23_mass;*/
  	

};
TightMuons tight_mu;

//////////////////////////////////////////////////////


struct DeltaR_MuonsTight{
	
   vector<Float_t> DR_lep_12;
   vector<Float_t> DR_lep_34;
   vector<Float_t> DR_lep_13;
   vector<Float_t> DR_lep_24;
   vector<Float_t> DR_lep_14;	
   vector<Float_t> DR_lep_23;	
	
};
DeltaR_MuonsTight drMuonsTight;


//////////////////////////////////////////////////////

struct Muon1{

   vector<Float_t> pT;
   vector<Float_t> eta;
   vector<Float_t> phi;
   vector<Float_t> mass;

};
Muon1 muon1;

/////////////////////////////////////////////////
struct Muon2{
   
   vector<Float_t> pT;
   vector<Float_t> eta;
   vector<Float_t> phi;
   vector<Float_t> mass;

};
Muon2 muon2;

/////////////////////////////////////////////////
struct Muon3{
   
   vector<Float_t> pT;
   vector<Float_t> eta;
   vector<Float_t> phi;
   vector<Float_t> mass;

};
Muon3 muon3;

/////////////////////////////////////////////////
struct Muon4{
   
   vector<Float_t> pT;
   vector<Float_t> eta;
   vector<Float_t> phi;
   vector<Float_t> mass;

};
Muon4 muon4;

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

struct ZcombMass{

   vector<Float_t> mZ_12;
   vector<Float_t> mZ_34;
   vector<Float_t> mZ_13;
   vector<Float_t> mZ_24;
   vector<Float_t> mZ_14;
   vector<Float_t> mZ_23;

};

ZcombMass mz;
////////////////////////////////////////////////////////////


struct Jets{
	
   vector<Float_t> jet_pt;
   vector<Float_t> jet_eta;
   vector<Float_t> jet_phi;
   vector<Float_t> loose_bjet_pt;
   vector<Float_t> loose_bjet_eta;
   vector<Float_t> loose_bjet_phi; 
   vector<Float_t> med_bjet_pt;
   vector<Float_t> med_bjet_eta;
   vector<Float_t> med_bjet_phi;
   vector<Float_t> tight_bjet_pt;
   vector<Float_t> tight_bjet_eta;
   vector<Float_t> tight_bjet_phi;
   vector<Int_t> jet_size;
   vector<Int_t> loose_bjet_size; 
   vector<Int_t> med_bjet_size;
   vector<Int_t> tight_bjet_size;
  
   vector<Float_t> v_all_loose_bjets_pt;
   vector<Float_t> v_all_loose_bjets_eta; 
   
	

};
Jets jets;

///////////////////////////////////////////////////////


struct DRLooseBjets4Muons{

   vector<Double_t> DR_LooseBJet_LooseMu1;
   vector<Double_t> DR_LooseBJet_LooseMu2;
   vector<Double_t> DR_LooseBJet_LooseMu3;
   vector<Double_t> DR_LooseBJet_LooseMu4;
   vector<Double_t> DR_LooseBJet_TightMu1;
   vector<Double_t> DR_LooseBJet_TightMu2;
   vector<Double_t> DR_LooseBJet_TightMu3;
   vector<Double_t> DR_LooseBJet_TightMu4;

};

DRLooseBjets4Muons drLooseBjet4Mu;

//////////////////////////////////////////////////////

struct DRMedBjets4Muons{

   vector<Double_t> DR_MedBJet_LooseMu1;
   vector<Double_t> DR_MedBJet_LooseMu2;
   vector<Double_t> DR_MedBJet_LooseMu3;
   vector<Double_t> DR_MedBJet_LooseMu4;
   vector<Double_t> DR_MedBJet_TightMu1;
   vector<Double_t> DR_MedBJet_TightMu2;
   vector<Double_t> DR_MedBJet_TightMu3;
   vector<Double_t> DR_MedBJet_TightMu4;

};

DRMedBjets4Muons drMedBjet4Mu;

///////////////////////////////////////////////////////


struct DRTightBjets4Muons{

   vector<Double_t> DR_TightBJet_LooseMu1;
   vector<Double_t> DR_TightBJet_LooseMu2;
   vector<Double_t> DR_TightBJet_LooseMu3;
   vector<Double_t> DR_TightBJet_LooseMu4;
   vector<Double_t> DR_TightBJet_TightMu1;
   vector<Double_t> DR_TightBJet_TightMu2;
   vector<Double_t> DR_TightBJet_TightMu3;
   vector<Double_t> DR_TightBJet_TightMu4;

};

DRTightBjets4Muons drTightBjet4Mu;


///////////////////////////////////////////////////////

struct MET{
	
   vector<Float_t> met_MET;
   vector<Float_t> met_eta;
   vector<Float_t> met_phi;

		
};
MET met;

//////////////////////////////////////////////////////

/*struct DeltaR_Muons{
	
   vector<Float_t> v_DR_mu1mu2;
   vector<Float_t> v_DR_mu3mu4;
   vector<Float_t> v_DR_mu1mu3;
   vector<Float_t> v_DR_mu2mu4;
   vector<Float_t> v_DR_mu1mu4;	
   vector<Float_t> v_DR_mu2mu3;	
	
};
DeltaR_Muons drMuons; */

///////////////////////////////////////////////////////

struct FourMuons_BeforeCUT{
	
   vector<Float_t> muon1_pt;
   vector<Float_t> muon2_pt;
   vector<Float_t> muon3_pt;
   vector<Float_t> muon4_pt;
   vector<Float_t> muon1_eta;	
   vector<Float_t> muon2_eta;
   vector<Float_t> muon3_eta;	
   vector<Float_t> muon4_eta;	
   vector<Float_t> muon1_phi;
   vector<Float_t> muon2_phi;
   vector<Float_t> muon3_phi;
   vector<Float_t> muon4_phi;	
	
};
FourMuons_BeforeCUT four_muons_beforeCut;

///////////////////////////////////////////////////////

struct DeltaR_Muons_bjet{
	
   vector<Double_t> v_DR_bmu1;
   vector<Double_t> v_DR_bmu2;
   vector<Double_t> v_DR_bmu3;
   vector<Double_t> v_DR_bmu4;
 		
};
DeltaR_Muons_bjet dR_bjets_mu_NoCUT;


//////////////////////////////////////////////////////


struct FourMuons_after_FullSelection{
	
   vector<Float_t> muon1_pt;
   vector<Float_t> muon2_pt;
   vector<Float_t> muon3_pt;
   vector<Float_t> muon4_pt;
   vector<Float_t> muon1_eta;	
   vector<Float_t> muon2_eta;
   vector<Float_t> muon3_eta;	
   vector<Float_t> muon4_eta;	
   vector<Float_t> muon1_phi;
   vector<Float_t> muon2_phi;
   vector<Float_t> muon3_phi;
   vector<Float_t> muon4_phi;	
	
};
FourMuons_after_FullSelection four_muons_afterFullSel;


////////////////////////////////////////////////////////


struct DeltaR_4Muons_after_FullSelection{
	
   vector<Float_t> v_DR_mu1mu2;
   vector<Float_t> v_DR_mu3mu4;
   vector<Float_t> v_DR_mu1mu3;
   vector<Float_t> v_DR_mu2mu4;
   vector<Float_t> v_DR_mu1mu4;	
   vector<Float_t> v_DR_mu2mu3;	
	
};
DeltaR_4Muons_after_FullSelection drMuons_FullSel;


////////////////////////////////////////////////////////

struct Za_OfSignal{
	
   vector<Float_t> Za_mass;
   vector<Float_t> Za_pt;
   vector<Float_t> Za_eta;
   	
	
};
Za_OfSignal za_sig;

////////////////////////////////////////////////////////



struct Zb_OfSignal{
	
   vector<Float_t> Zb_mass;
   vector<Float_t> Zb_pt;
   vector<Float_t> Zb_eta;
   	
	
};
Zb_OfSignal zb_sig;

////////////////////////////////////////////////////////

struct b1_OfSignal{
	
   vector<Float_t> b1_mass;
   vector<Float_t> b1_pt;
   vector<Float_t> b1_eta;
   vector<Float_t> DR_b1_mu1;
   vector<Float_t> DR_b1_mu2;
   vector<Float_t> DR_b1_mu3;
   vector<Float_t> DR_b1_mu4;
   	
	
};
b1_OfSignal b1_sig;

////////////////////////////////////////////////////////

struct b2_OfSignal{
	
   vector<Float_t> b2_mass;
   vector<Float_t> b2_pt;
   vector<Float_t> b2_eta;
   vector<Float_t> DR_b2_mu1;
   vector<Float_t> DR_b2_mu2;
   vector<Float_t> DR_b2_mu3;
   vector<Float_t> DR_b2_mu4;
   	
	
};
b2_OfSignal b2_sig;

//////////////////////////////////////////////////////////

struct FirstHiggsOfSignal{
	
   vector<Float_t> v_h1_mass;
   vector<Float_t> v_h1_pt;
   vector<Float_t> v_h1_eta;
		
}; 
FirstHiggsOfSignal smhiggs1;

///////////////////////////////////////////////////////////

struct SecondHiggsOfSignal{
	
   vector<Float_t> v_h2_mass;
   vector<Float_t> v_h2_pt;
   vector<Float_t> v_h2_eta;
		
}; 
SecondHiggsOfSignal smhiggs2;

//////////////////////////////////////////////////////////

struct BSM_HiggsOfSignal{
	
   vector<Float_t> v_H_mass;
   vector<Float_t> v_H_pt;
   vector<Float_t> v_H_eta;
		
}; 
BSM_HiggsOfSignal heavyHiggs;

//////////////////////////////////////////////////////////


TTree* new_tree  = new TTree("output_demo", "output_demo");


TBranch*   mass_mu12pair_noCut;
TBranch*   mass_mu34pair_noCut;
TBranch*   mass_mu13pair_noCut;
TBranch*   mass_mu24pair_noCut;
TBranch*   mass_mu14pair_noCut;
TBranch*   mass_mu23pair_noCut;
TBranch*   all_bjet_loose_pT_noCut;
TBranch*   all_bjet_loose_eta_noCut;
TBranch*   b1_jet_loose_pT_afterCut;
TBranch*   b1_jet_loose_eta_afterCut;
TBranch*   b1_jet_loose_Charge_afterCut;
TBranch*   b2_jet_loose_pT_afterCut;
TBranch*   b2_jet_loose_eta_afterCut;
TBranch*   b2_jet_loose_Charge_afterCut;
TBranch*   all_bjet_tight_pT_noCut;
TBranch*   all_bjet_tight_eta_noCut;
TBranch*   b1_jet_tight_pT_afterCut;
TBranch*   b1_jet_tight_eta_afterCut;
TBranch*   b1_jet_tight_Charge_afterCut;
TBranch*   b2_jet_tight_pT_afterCut;
TBranch*   b2_jet_tight_eta_afterCut;
TBranch*   b2_jet_tight_Charge_afterCut;
TBranch*   b_Gen_Muon;
TBranch*   b_Gen_Bjets;
TBranch*   b_LooseMuons;
TBranch*   b_deltaR_muons_loose;
TBranch*   b_TightMuons;
TBranch*   b_deltaR_muons_tight;
TBranch*   b_muon1;
TBranch*   b_muon2;
TBranch*   b_muon3;
TBranch*   b_muon4;
TBranch*   b_ZcombMass;
TBranch*   b_Jets;
TBranch*   b_DR_LooseBJets_4Mu;
TBranch*   b_DR_MedBJets_4Mu;
TBranch*   b_DR_TightBJets_4Mu;
TBranch*   b_MET;
TBranch*   b_4Muons_beforeCut;
TBranch*   b_deltaR_muons;
TBranch*   b_dR_bjet_mu_NoCUT;
TBranch*   b_4Muons_FullSel;
TBranch*   b_dR_4Muons_FullSel;
TBranch*   b_Za_signal;
TBranch*   b_Zb_signal;
TBranch*   b_b1_signal;
TBranch*   b_b2_signal;
TBranch*   b_1st_Higgs_Ofsignal;
TBranch*   b_2nd_Higgs_Ofsignal;
TBranch*   b_BSM_Higgs_Ofsignal;
TBranch*   b_dr_b1b2;
TBranch*   b_dr_h1h2;
TBranch*   b_mass_4l_afterZ1Z2windowcut;
TBranch*   b_mass_4l_afterfullbjetcuts;
TBranch*   b_za_sig_Nozzwin;
TBranch*   b_zb_sig_Nozzwin;
TBranch*   b_za_sig_z1win_cut;
TBranch*   b_zb_sig_z1win_cut;
TBranch*   b_za_sig_zzwin_cut;
TBranch*   b_zb_sig_zzwin_cut;
TBranch*   b_za_sig_m4Mu_cut;
TBranch*   b_zb_sig_m4Mu_cut;
TBranch*   b_nLoose_bjets;
TBranch*   b_nMed_bjets;
TBranch*   b_m_4l_m4Mu_cut;
TBranch*   b_m_4l_tmp1_bjet_cut;
TBranch*   b_m_2bj_tmp1_bjet_cut;

#endif
