#include "TCanvas.h"
#include "TPad.h"
#include "TPaveLabel.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TFrame.h"

const Int_t NMAX = 20;
Int_t NLOOP;
Float_t SIGN[NMAX], INTG_LUM[NMAX];

void Lum_calc(Float_t, Float_t, Float_t, Float_t, Float_t, Float_t);

//__________________________________________________________________
void LumForDiscovery()
{

   Float_t lum_mc;
   Float_t nsg;
   Float_t nbg;
   Float_t sig1;
   Float_t sig2;
   Float_t delp;

   // Create a new canvas.
   TCanvas *c1 = new TCanvas("LumForDiscovery","Monte Carlo Study of Heavy Higgs discovery @ 14 TeV",10,40,800,600);
   c1->Range(0,0,25,18);
   c1->SetFillColor(0);
   
   TLatex *t = new TLatex();
   t->SetTextFont(32);
   t->SetTextColor(1);
   t->SetTextSize(0.03);
   t->SetTextAlign(12);


   TPad *pad1 = new TPad("pad1","This is pad1",0.0001,0.0001,0.99,0.99,0);
  // TPad *pad2 = new TPad("pad2","This is pad2",0.52,0.02,0.99,0.83,33);
   pad1->Draw();
   //pad2->Draw();

   
   //====================
   //       BMP1 - 2022/09/15
   //====================
   
   //lum_mc = 362.97/0.12099;         // weighted Ngen   w = lumi_data * xs/Ngen = 3000 * 0.12099/1.e+05 = 3.6297.e-03 
   lum_mc = 1610/0.5367;
   //nsg  = 65.1095;   // @ old code @ xs = 0.12099 fb @ LO @ wt = 0.0036
   nsg = 291.1842;      // @ old code @ xs = 0.5367 fb @ NNLO with k-factor corrections @ wt = 0.0161 
  //nsg = 261.8504;     // @ new code with modified slection of b-jets from Loose @ xs = 0.5367 fb @ NNLO with k-factor (= 4.427) corrections @ wt = 0.0161 
   nbg  = 18.72596;  //true value
   sig1 = 1.;
   sig2 = 7.;
   delp = 0.3;
   Lum_calc(lum_mc, nsg, nbg, sig1, sig2, delp);
   pad1->cd();
   pad1->Range(-0.255174,-19.25,2.29657,-6.75);
   pad1->SetGridx(1);
   pad1->SetGridy(1);
   pad1->SetTickx(1);
   pad1->SetTicky(1);
   pad1->SetLogx(0);
   pad1->SetLogy(1);

   // create a 2-d histogram to define the range
   pad1->DrawFrame(2,1e02,8,1e05);
   pad1->GetFrame()->SetFillColor(0);
   t = new TLatex();
   t->SetNDC();
   t->SetTextFont(72);
   t->SetTextColor(1);
   t->SetTextSize(0.04);
   t->SetTextAlign(13);
   //t->DrawLatex(0.2,0.8,"H discovery potential @ 14 TeV");
   t->DrawLatex(0.2,0.8,"HL-LHC Simulation(Delphes)");
   t->DrawLatex(0.2,0.7,"#sqrt{s} = 14 TeV");
   t->DrawLatex(0.2,0.6,"gg #rightarrow H #rightarrow hh #rightarrow b#bar{b} 4#mu");

   //t->SetTextSize(0.04);
   t->SetTextColor(1);
   t->SetTextAngle(90);
   t->DrawLatex(0.03,0.50,"Luminosity [fb^{-1}]");

   //t->SetTextSize(0.04);
   t->SetTextAngle(0);
   t->SetTextColor(kBlue);
   t->DrawLatex(0.7,0.4,"BP1");
   t->SetTextColor(kMagenta);
   t->DrawLatex(0.7,0.35,"BP2");
   t->SetTextColor(kGreen+2);
   t->DrawLatex(0.7,0.3,"BP3");
   //t->SetTextColor(kGreen+2);
   //t->DrawLatex(0.7,0.15,"BP4"); 
 

   //t->SetTextSize(0.04);
   t->SetTextColor(1);
   t->DrawLatex(0.68,0.05,"Significance");


   TGraph *gr1 = new TGraph(NLOOP,SIGN,INTG_LUM);

   gr1->SetMaximum(1.e+10);
   gr1->SetMinimum(1.e+01);
   
   gr1->SetLineColor(kBlue);
   gr1->SetLineWidth(2);
   gr1->SetMarkerColor(kBlue);
   gr1->SetMarkerStyle(21);
   gr1->SetMarkerSize(1.1);
   gr1->Draw("LP");
  
   

   //====================
   //       BMP2 - 2022/09/15
   //====================
   
   //lum_mc = 349.89/0.11663;
   lum_mc = 1548.963/0.5163;  // @ old code @ xs = 0.5163 fb @ NNLO with k-factor corrections @ wt = 0.015489 
   //nsg  = 59.232;   // old value must be modfied to new one (regenerate sample)
   nsg = 262.0784;     // @ old code @ xs =0.5163 fb @ NNLO with k-factor corrections @ wt = 0.015489 
   nbg  = 18.72596;
   
   Lum_calc(lum_mc,nsg, nbg, sig1, sig2, delp);


   TGraph *gr2 = new TGraph(NLOOP,SIGN,INTG_LUM);

   gr2->SetMaximum(1.e+09);
   gr2->SetMinimum(1.e+03);
   
   gr2->SetLineColor(kMagenta);
   gr2->SetLineWidth(2);
   gr2->SetMarkerColor(kMagenta);
   gr2->SetMarkerStyle(29);
   gr2->SetMarkerSize(1.5);
   gr2->Draw("LP");
   


   //====================
   //       BMP3 - 2022/09/15
   //====================
   
   //lum_mc = 351.9/0.1173;
   lum_mc = 1560/0.5193; // @ old code @ xs = 0.5193 fb @ NNLO with k-factor corrections @ wt = 0.0156 
  // nsg  = 60.709;
   nsg = 269.0878;     // @ old code @ xs = 0.5193 fb @ NNLO with k-factor corrections @ wt = 0.0156
   nbg  = 18.72596;
   
   Lum_calc(lum_mc, nsg, nbg, sig1, sig2, delp);

   TGraph *gr3 = new TGraph(NLOOP,SIGN,INTG_LUM);

   gr3->SetMaximum(1.e+09);
   gr3->SetMinimum(1.e+03);

   gr3->SetLineColor(kGreen+2);
   gr3->SetLineWidth(2);
   gr3->SetMarkerColor(kGreen+2);
   gr3->SetMarkerStyle(22);
   gr3->SetMarkerSize(1.1);
   gr3->Draw("LP");
   
   

   c1->Modified();
   c1->Update();
}


//Lum[j][k] = (LUM_MC*(pow(Sigma[j],2)))/(2*( ((sig_n[k]+BG_tot)*log(1+ sig_n[k]/BG_tot )) - sig_n[k]));
void Lum_calc(Float_t LUM_MC, Float_t NSG, Float_t NBG, Float_t Sig1, Float_t Sig2, Float_t DELP)
{
  Int_t I;

  Float_t Nbg, Nsg, Lum_mc;
  Float_t S;
  Float_t N;

  N = NSG + NBG;
  Nbg= NBG;
  Nsg= NSG;
  Lum_mc= LUM_MC;
  NLOOP = (Sig2-(Sig1 -1 ))/DELP;

  for (I=0; I<NLOOP;I++) {
     SIGN[I]=Sig1+I*DELP;

     //S    = sqrt( 2*( ((NSG+NBG)*log(1+ NSG/NBG )) - NSG));

      S  = sqrt( (2*N*log(1+ (NSG/NBG) )) - (2*NSG) );  // applied significance formula   for NBG < 100      @ 2022/9/15
     
     
     // S    = ( NSG) / (sqrt(NSG));
     //S    = ( NSG) / (sqrt(NBG));
     
     INTG_LUM[I]=(Lum_mc*(pow(SIGN[I],2)))/S;
     
     cout << "--> Integrated lum.  = " << INTG_LUM[I]<< " @ Significance = " << SIGN[I]<<" and Ns = "<< NSG << endl;
  }
}

