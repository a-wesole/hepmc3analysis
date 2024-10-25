#include <iostream>
#include "TRandom.h"
#include "TFile.h"
#include <iostream>
#include <random>
#include "TFile.h"
#include "TNtuple.h"
#include "TTree.h"
#include <TBranch.h>
#include <vector>
#include <TRandom3.h>
#include "TROOT.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TMath.h"

//parameters for D0 candidates, manually fitted suning temp.C output 
const float D0_signal_parameters[] = {3.0e+04, 1.86502e+00, 2.84804e-02, 1.19824e-02, 2.40622e-01}; //p0, p1, p2, p3, p4
const float D0_swap_paremeters[] = {6.33946e+03, 1.87550e+00, 1.16808e-01}; //const, mean, sigmacfor gaus
const float D0_bkg_parameters[] = {-365392, 1.54298e+06};  // p0 and p1 after fitting with pol1
//const float D0_bkg_parameters[] = {-365392, 2.1*365392};  // p0 and p1 after fitting with pol1
//const float D0_signal_parameters[] = {1.0, 1.86502e+00, 2.84804e-02, 1.19824e-02, 2.40622e-01}; //p0, p1, p2, p3, p4
//const float D0_swap_paremeters[] = {1.0, 1.87550e+00, 1.16808e-01}; //const, mean, sigmacfor gaus
//const float D0_bkg_parameters[] = { -365392, 0.0};  // p0 and p1 after fitting with pol1

//for Dbar candidates, manually fitted usuing temp.C
const float Dbar_signal_parameters[] = {3.1e+04, 1.86495e+00, 2.72819e-02, 1.16558e-02, 2.33265e-01}; // p0, p1, p2, p3, p4 fitted with double gaus as in data   
const float Dbar_swap_parameters[] = {6.25532e+03, 1.87519e+00, 1.19091e-01}; //const, mean and sigma for gaus 
const float Dbar_bkg_parameters[] = { -353807,  1.52099e+06};// p0 and p1 after fitting with pol1
//const float Dbar_bkg_parameters[] = { -353807,  2.1*353807};// p0 and p1 after fitting with pol1
//const float Dbar_signal_parameters[] = {1.0, 1.86495e+00, 2.72819e-02, 1.16558e-02, 2.33265e-01}; // p0, p1, p2, p3, p4 fitted with double gaus as in data   
//const float Dbar_swap_parameters[] = {1.0, 1.87519e+00, 1.19091e-01}; //const, mean and sigma for gaus 
//const float Dbar_bkg_parameters[] = {-353807, 0.0};// p0 and p1 after fitting with pol1

double histo_entries = 1.0e+05;
const float nbins = 16.0, xmin = 1.7, xmax=2.0;
const double avg_mass = 1.865, std_mass = 1.902e-02, fg_lower = avg_mass-3*std_mass, fg_upper=avg_mass+3*std_mass; // July 30 updates
const double fit_range_low = 1.6, fit_range_high = 2.1, D0_mass = 1.8648;


using namespace std;

void FitMassPlotsTF1(TF1 *fitFunction, TH1D *MassPlot) {
    for (int i = 0; i < 20; i++) {
        MassPlot->Fit(fitFunction, "Lq");
    }
    MassPlot->Fit(fitFunction, "q R");
}


void FitMassPlots(TF1 *f, TH1D *SignalMass, TH1D *SignalAndSwapMass, TH1D *MassPlot, TF1 *swap_fit, TH1D *swap_mass) {

  swap_fit->SetLineColor(kCyan);
  TCanvas *c1 = new TCanvas("c1", "c1", 1000, 1000);
  c1->Divide(2,2);
  c1->cd(1);
  f->SetLineColor(2);
  f->SetLineWidth(2);
  f->SetParNames("normalization", "mean mass", "sigma 1", "sigma 2", "gaus ratio", "signalfrac", "smearing", "sigma swap", "const", "x", "x*x");
  f->SetParName(11, "x*x*x");
  f->SetParName(12, "raio guas 2/3");
  f->SetParName(13, "swap gaus ratio");
  f->SetParName(14, "sigma 3");
  f->SetParName(15, "sigma2 swap");
  f->SetParameter(0, 100);    // signal plus swap
  f->SetParameter(1, D0_mass); // mean mass
  f->SetParameter(2, 0.01);    // sigma 1
  f->SetParameter(3, 0.005);   // sigma 2
  f->SetParameter(14, 0.03); //sigma 3 signal 
  f->SetParameter(4, 0.5);     // ratio between first and second guassian
  f->SetParameter(12, 0.5); //ratio gaus 2&3 

  f->FixParameter(5, 1.0); // fraction that is signal
  f->FixParameter(6, 0); // scaling factor----------always 0 in MC
   f->FixParameter(7,1.0); //sigma of swap-------------- does not really mater here as yield is fix to 0
  // f->SetParameter(8,D0_mass);//mean of swap
  f->FixParameter(8, 0);  // polynomial
  f->FixParameter(9, 0);  // polynomial
  f->FixParameter(10, 0); // polynomail
  f->FixParameter(11, 0); // poly
  f->FixParameter(13, 0); //ratio btw gaussian swap 
  f->FixParameter(15, 1); //sigma2 swap

  //f->SetParLimits(2, 0.01, 0.1);   // sigma 1
  f->SetParLimits(1, 1.85, 1.88);  // sigma 1
  //f->SetParLimits(3, 0.001, 0.05); // sigma 2
  f->SetParLimits(4, 0, 1);        // amplitude/ratio
  f->SetParLimits(12, 0, 1);        // amplitude/ratio

  f->SetParameter(1, 1.8648); // mean mass

  for (int i=0; i<29; i++){
  SignalMass->Fit(f, "q", "", fit_range_low, fit_range_high);
  }
  SignalMass->Fit(f, "M", "", fit_range_low, fit_range_high);
  SignalMass->Draw("ep");
  c1->cd(2);
  swap_mass->Draw();
  swap_fit->Draw("same");
  c1->cd(3);

  f->FixParameter(1, f->GetParameter(1)); // mean mass
  f->FixParameter(2, f->GetParameter(2)); // sigma 1
  f->FixParameter(3, f->GetParameter(3)); // sigma 2
  f->FixParameter(4, f->GetParameter(4)); // ratio/amplitude
  f->FixParameter(12, f->GetParameter(12)); // ratio guaes 2/3
  f->FixParameter(14, f->GetParameter(14)); // sigma 3

  /*
  f->ReleaseParameter(5); // fraction of signal
  f->ReleaseParameter(7); // sigma swap1
  f->ReleaseParameter(15); // sigma swap2
  f->ReleaseParameter(13); // ratio
  f->SetParLimits(5, 0, 1);        // fraction of signal
  f->SetParLimits(13, 0, 1);        // fraction of signal

  f->SetParameter(7, 0.1); // sigma swap
  f->SetParameter(15, 0.01); // sigma swap
  f->SetParameter(13, 0.7); //ratio
  f->SetParLimits(7, 0.0, 1.9); //sigma1
  f->SetParLimits(12, 0, 0.5); //sigma2
  //f->SetParameter(5, 0.99999); // sigma swap
  */
  /*
  f->SetParameter(7, swap_fit->GetParameter(2)); //sigma 1 
  f->SetParameter(15, swap_fit->GetParameter(3)); //sigma 2 
  f->SetParameter(13, swap_fit->GetParameter(4)); //raio 

  f->ReleaseParameter(5);
  f->ReleaseParameter(7);
  f->ReleaseParameter(15);
  f->ReleaseParameter(13);
  f->SetParLimits(5, 0, 1);
  f->SetParLimits(7, 0.05, 0.2); //sigma1
  //f->SetParLimits(13, swap_fit->GetParameter(4), 1);
  */
  for (int i = 0; i < 29; i++)
  {
    SignalAndSwapMass->Fit(f, "q", "", fit_range_low, fit_range_high);
  }
  SignalAndSwapMass->Fit(f, "M", "", fit_range_low, fit_range_high);
  
  
  TF1 *swap = new TF1("swap", "[0]*(1-[5])*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6])) + (1-[4])*(TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6]))))", fit_range_low, fit_range_high);
  swap->FixParameter(0, f->GetParameter(0)); //norm
  swap->FixParameter(1, f->GetParameter(1)); //mean mass
  swap->FixParameter(2, f->GetParameter(7));  //sigma1
  swap->FixParameter(3, f->GetParameter(15)); //sigma2
  swap->FixParameter(4, f->GetParameter(13)); //ratio
  swap->FixParameter(5, f->GetParameter(5)); //signal fraction
  swap->FixParameter(6, f->GetParameter(6)); //smearing
  swap->SetLineColor(8);

  SignalAndSwapMass->Draw("ep");
  swap_fit->Draw("same");
  swap->Draw("same");
  c1->Update();
  c1->cd(4);
  

  f->FixParameter(0, f->GetParameter(0)); // fraction of signal
  f->FixParameter(5, f->GetParameter(5)); // fraction of signal
  f->FixParameter(15, f->GetParameter(15)); // fraction of signal
  f->FixParameter(13, f->GetParameter(13)); // fraction of signal
  f->FixParameter(7, f->GetParameter(7)); // sigma swap

  f->ReleaseParameter(8);  // poly
  f->ReleaseParameter(9);  // poly
  f->ReleaseParameter(10); // poly
  //f->ReleaseParameter(11); // 12

  MassPlot->Fit(f, "q", "", fit_range_low, fit_range_high);
  MassPlot->Fit(f, "q", "", fit_range_low, fit_range_high);
  f->ReleaseParameter(1);    // mean mass ------allow data to have different mass peak mean than MC
  f->ReleaseParameter(6);    // ratio ------- allow data to have different peak width than MC
  f->SetParameter(6, 0);     // scaling factor
  f->SetParLimits(6, -1, 1); // scaling factor
  for (int i=0; i<10; i++){
  MassPlot->Fit(f, "L q", "", fit_range_low, fit_range_high);
  }
  MassPlot->Fit(f, "L m", "", fit_range_low, fit_range_high); ///

  TF1 *signal = new TF1("signal", "([6]*[7]*([0]*TMath::Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2]) + [5]*((1-[0]))*TMath::Gaus(x,[1],[3])/(sqrt(2*3.14159)*[3]) + (1-[5]-[0])*TMath::Gaus(x,[1],[4])/(sqrt(2*3.14159)*[4])))", 1.7, 2.1);
  signal->FixParameter(0, f->GetParameter(4)); //ratio btw 1 and 2
  signal->FixParameter(1, f->GetParameter(1)); //mean
  signal->FixParameter(2, f->GetParameter(2)); //sigma 1
  signal->FixParameter(3, f->GetParameter(3)); //sihgma 2
  signal->FixParameter(4, f->GetParameter(14)); //sigma 3
  signal->FixParameter(5, f->GetParameter(12));//ratio btw gaus 2 & 3
  signal->FixParameter(6, f->GetParameter(0)); //scaling
  signal->FixParameter(7, f->GetParameter(5)); //signal fraction

  

  TF1 *background = new TF1("background", "[8] + [9]*x + [10]*x*x + [11]*x*x*x", fit_range_low, fit_range_high);
  background->FixParameter(8, f->GetParameter(8));
  background->FixParameter(9, f->GetParameter(9));
  background->FixParameter(10, f->GetParameter(10));
  background->FixParameter(11, f->GetParameter(11));

  MassPlot->SetMinimum(0.5);
  MassPlot->Draw("EP");
  swap_fit->Draw("same");
  signal->SetLineColor(4);
  signal->Draw("same");
  swap->SetLineColor(8);
  swap->Draw("same");
  background->SetLineColor(1);
  background->Draw("same");
  f->SetLineColor(2);
  f->Draw("same");

  c1->Print("output.pdf"); // Print to the already opened PDF
  c1->Delete();
}


void updateCreateTemplates(){

  TString outfile = TString("TF1_outputs_skippingswap.root");
  //TString outfile = TString("garbage.root");
  TFile *results = new TFile(outfile, "recreate");

  // Create a temporary canvas to open and close the PDF
  TCanvas *c_temp = new TCanvas();

  // Open the PDF file using the temporary canvas
  c_temp->Print("output.pdf["); // This opens the PDF

  double F1F2Entries, F1SB2Entries, SB1F2Entries, SB1SB2Entries;
  cout << "working to open file..." << endl;
  //TFile *infile = TFile::Open("/home/awesole/pythia8/phi_corr_out_11July_newProductionFinal.root", "READ");
  //TString inputfile1 = "TH2F_output_1k_noswap.root";
  TString inputfile1 = "data_files/TH2F_skippingswap.root";
  TFile *inf1 = TFile::Open(inputfile1);

  for (int a=6; a<7; a++) {
  TH2D *M1M2Mass = (TH2D *)inf1->Get(Form("M1M2Mass_%d",a)); // Full Range of F1F2 Mass
  TH1D *M1M2hx = M1M2Mass->ProjectionX("M1M2hx");
  TH1D *M1M2hy = M1M2Mass->ProjectionY("M1M2hy");

  TH2D *S1S2Mass = (TH2D *)inf1->Get(Form("S1S2Mass_%d",a)); // Full Range of F1F2 Mass
  S1S2Mass->SetMinimum(0.5);
  TH1D *S1S2hx = S1S2Mass->ProjectionX("S1S2hx");
  TH1D *S1S2hy = S1S2Mass->ProjectionY("S1S2hy");

  TH2D *B1SW2Mass = (TH2D *)inf1->Get(Form("B1SW2Mass_%d",a));
  TH2D *B1S2Mass = (TH2D *)inf1->Get(Form("B1S2Mass_%d",a));
  TH2D *SW1S2Mass = (TH2D *)inf1->Get(Form("SW1S2Mass_%d",a)); // Full Range of F1F2 Mass
  SW1S2Mass->SetMinimum(0.5);
  SW1S2Mass->RebinX(1);
  SW1S2Mass->RebinY(1);
  TH2D *SW1B2Mass = (TH2D *)inf1->Get(Form("SW1B2Mass_%d",a));
  TH2D *S1B2Mass = (TH2D *)inf1->Get(Form("S1B2Mass_%d",a));
  TH2D *S1SW2Mass = (TH2D *)inf1->Get(Form("S1SW2Mass_%d",a)); // Full Range of F1F2 Mass
  S1SW2Mass->SetMinimum(0.5);
  S1SW2Mass->RebinX(1);
  S1SW2Mass->RebinY(1);

  S1S2hx->Add(S1B2Mass->ProjectionX());
  S1S2hx->Add(S1SW2Mass->ProjectionX());
  S1S2hy->Add(B1S2Mass->ProjectionY());
  S1S2hy->Add(SW1S2Mass->ProjectionY());

  TH2D *SignalSwap12Mass = (TH2D *)inf1->Get(Form("SignalSwapMass_%d",a)); // Full Range of F1F2 Mass
  SignalSwap12Mass->SetMinimum(0.5);
  TH1D *SSWhx = SignalSwap12Mass->ProjectionX("SSWhx");
  TH1D *SSWhy = SignalSwap12Mass->ProjectionY("SSWhy");

  SSWhx->Add(SW1B2Mass->ProjectionX());
  SSWhx->Add(S1B2Mass->ProjectionX());
  SSWhy->Add(B1SW2Mass->ProjectionY());
  SSWhy->Add(B1S2Mass->ProjectionY());

  TH2D *SwapOnlyMass = (TH2D *)inf1->Get(Form("SW1SW2Mass_%d",a)); // Full Range of F1F2 Mass
  SwapOnlyMass->SetMinimum(0.5);
  SwapOnlyMass->RebinX(1);
  SwapOnlyMass->RebinY(1);
  SwapOnlyMass->SetTitle("SW1SW2 Mass");
  TH1D *SWhx = SwapOnlyMass->ProjectionX("SWhx");
  TH1D *SWhy = SwapOnlyMass->ProjectionY("SWhy");

  SWhx->Add(SW1B2Mass->ProjectionX());
  SWhx->Add(SW1S2Mass->ProjectionX());
  SWhy->Add(B1SW2Mass->ProjectionY());
  SWhy->Add(S1SW2Mass->ProjectionY());

  TH2D *BkgOnlyMass = (TH2D *)inf1->Get(Form("B1B2Mass_%d",a)); // Full Range of F1F2 Mass
  BkgOnlyMass->SetMinimum(0.5);
  cout << "File Successfully Opened!" << endl;

  TString name1 = TString::Format("F1_%d", a);
  TString name2 = TString::Format("F2_%d", a);

  TF1 *SwapSW1 = new TF1("SwapSW1", "[0]*([4]*TMath::Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2]) + (1-[4])*TMath::Gaus(x,[1],[3])/(sqrt(2*3.14159)*[3]))", fit_range_low, fit_range_high);
  // TF1 *SwapSW1 = new TF1("SwapSW1", "[0]* TMath::Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])", 1.7, 2.1);
  SwapSW1->SetLineColor(2);
  // SwapSW1->SetLineWidth(1);
  SwapSW1->SetParameter(0, 100);     // normalization
  SwapSW1->SetParameter(1, D0_mass); // mean
  SwapSW1->SetParameter(2, 0.1);     // sigma1
  SwapSW1->SetParameter(3, 0.01);    //sigma 2
  SwapSW1->SetParameter(4, 0.7);     //ratio
  SwapSW1->SetParLimits(1, 1.8, 1.9); //mass
  SwapSW1->SetParLimits(2, 0.0, 1.9); //sigma1
  SwapSW1->SetParLimits(3, 0.0, 0.5); //sigma2
  SwapSW1->SetParLimits(4, 0.0, 1.0); //ratio
  FitMassPlotsTF1(SwapSW1, SWhx);

  TF1 *SwapSW2 = new TF1("SwapSW2", "[0]*([4]*TMath::Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2]) + (1-[4])*TMath::Gaus(x,[1],[3])/(sqrt(2*3.14159)*[3]))", fit_range_low, fit_range_high);
  // TF1 *SwapSW1 = new TF1("SwapSW1", "[0]* TMath::Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])", 1.7, 2.1);
  SwapSW2->SetLineColor(2);
  // SwapSW1->SetLineWidth(1);
  SwapSW2->SetParameter(0, 100);     // normalization
  SwapSW2->SetParameter(1, D0_mass); // mean
  SwapSW2->SetParameter(2, 0.1);     // sigma1
  SwapSW2->SetParameter(3, 0.01);
  SwapSW2->SetParameter(4, 0.7);
  SwapSW2->SetParLimits(1, 1.8, 1.9);
  SwapSW2->SetParLimits(2, 0.0, 1.9);
  SwapSW2->SetParLimits(3, 0.0, 0.5);
  SwapSW2->SetParLimits(4, 0.0, 1.0);
  FitMassPlotsTF1(SwapSW2, SWhy);

  //TF1 *F1 = new TF1(name1, "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 + [6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])))+(1-[5])*TMath::Gaus(x,[1],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6])))+ [8] + [9]*x + [10]*x*x + [11]*x*x*x   ", fit_range_low, fit_range_high);
  TF1 *F1 = new TF1(name1, "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 + [6]))+(1-[4])*[12]*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])) + (1-[4])*(1-[12])*TMath::Gaus(x,[1],[14]*(1.0 +[6]))/(sqrt(2*3.14159)*[14]*(1.0 +[6])))+(1-[5])*(([13]*TMath::Gaus(x,[1],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) + (1-[13])*(TMath::Gaus(x,[1],[15]*(1.0 +[6]))/(sqrt(2*3.14159)*[15]*(1.0 +[6])))))+ [8] + [9]*x + [10]*x*x + [11]*x*x*x", fit_range_low, fit_range_high);

  //  TF1 *F1 = new TF1(name1, "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 + [6]))+(1-[4])*[12]*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])) + (1-[4]-[12])*TMath::Gaus(x,[1],[14]*(1.0 +[6]))/(sqrt(2*3.14159)*[14]*(1.0 +[6])))+(1-[5])*([13]*TMath::Gaus(x,[1],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) + (1-[13])*(TMath::Gaus(x,[1],[15]*(1.0 +[6]))/(sqrt(2*3.14159)*[15]*(1.0 +[6]))))+ [8] + [9]*x + [10]*x*x + [11]*x*x*x", fit_range_low, fit_range_high);

  //TF1 *F1 = new TF1(name1, "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 + [6]))+(1-[4])*[12]*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])) + (1-[4]-[12])*TMath::Gaus(x,[1],[14]*(1.0 +[6]))/(sqrt(2*3.14159)*[14]*(1.0 +[6]))))+(1-[5])*([13]*TMath::Gaus(x,[1],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) + (1-[13])*(TMath::Gaus(x,[1],[15]*(1.0 +[6]))/(sqrt(2*3.14159)*[15]*(1.0 +[6])))+ [8] + [9]*x + [10]*x*x + [11]*x*x*x   ", fit_range_low, fit_range_high);
  //TF1 *F1 = new TF1(name1, "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 + [6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6]))) +(1-[5])*([13]*TMath::Gaus(x,[1],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) + (1-[13])*(TMath::Gaus(x,[1],[12]*(1.0 +[6]))/(sqrt(2*3.14159)*[12]*(1.0 +[6]))))+ [8] + [9]*x + [10]*x*x + [11]*x*x*x   ", fit_range_low, fit_range_high);
  FitMassPlots(F1, S1S2hx, SSWhx, M1M2hx, SwapSW1, SWhx);

  TF1 *F2 = new TF1(name2, "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 + [6]))+(1-[4])*[12]*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])) + (1-[4])*(1-[12])*TMath::Gaus(x,[1],[14]*(1.0 +[6]))/(sqrt(2*3.14159)*[14]*(1.0 +[6])))+(1-[5])*(([13]*TMath::Gaus(x,[1],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) + (1-[13])*(TMath::Gaus(x,[1],[15]*(1.0 +[6]))/(sqrt(2*3.14159)*[15]*(1.0 +[6])))))+ [8] + [9]*x + [10]*x*x + [11]*x*x*x", fit_range_low, fit_range_high);
  //TF1 *F2 = new TF1(name2, "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 + [6]))+(1-[4])*[12]*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])) + (1-[4]-[12])*TMath::Gaus(x,[1],[14]*(1.0 +[6]))/(sqrt(2*3.14159)*[14]*(1.0 +[6]))))+(1-[5])*([13]*TMath::Gaus(x,[1],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) + (1-[13])*(TMath::Gaus(x,[1],[15]*(1.0 +[6]))/(sqrt(2*3.14159)*[15]*(1.0 +[6])))+ [8] + [9]*x + [10]*x*x + [11]*x*x*x   ", fit_range_low, fit_range_high);
  //TF1 *F2 = new TF1(name2, "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 + [6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])))+(1-[5])*TMath::Gaus(x,[1],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6])))+ [8] + [9]*x + [10]*x*x + [11]*x*x*x   ", fit_range_low, fit_range_high);
  //TF1 *F2 = new TF1(name2, "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 + [6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6]))) +(1-[5])*([13]*TMath::Gaus(x,[1],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) + (1-[13])*(TMath::Gaus(x,[1],[12]*(1.0 +[6]))/(sqrt(2*3.14159)*[12]*(1.0 +[6]))))+ [8] + [9]*x + [10]*x*x + [11]*x*x*x   ", fit_range_low, fit_range_high);
  FitMassPlots(F2, S1S2hy, SSWhy, M1M2hy, SwapSW2, SWhy);

  results->cd();
  F1->Write();
  F2->Write();
  }

}



