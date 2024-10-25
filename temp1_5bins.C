#include <iostream>
#include <fstream>
#include <random>
#include <Pythia8/Pythia.h>
#include "TFile.h"
#include "TNtuple.h"
#include <TTree.h>
#include "TTree.h"
#include <TBranch.h>
#include <vector>
#include <TH1F.h>
#include <TH2F.h>
#include <TRandom3.h>
#include "TROOT.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TH1F.h"
#include <TF1.h>
#include <TArrayF.h>


const float upper_limit = 3.14159, lower_limit = -3.14159;                  // bounds to transition delta phi from (-2pi, 2pi) to (-0.5pi, 1.5*pi)
const float nbins_mass = 1000.0, xmin_mass = 1.6, xmax_mass = 2.1;          // for mass plots
const float nbins = 200.0, xmin = -2.0, xmax = 5.0, ymin = 0.0, ymax = 1.5; // for dphi plots
// const double avg_mass = 1.865, std_mass = 1.902e-02; // July 30 updates
// const double avg_mass = 1.86502, std_mass = 0.0147844;                       // July 30 updates
// const double avg_mass = 1.8648, std_mass = 1.40e-02; // July 30 updates
const float pT_minimum = 0.5; //   GeV/c
// the values below are corresponding to calculations on F1 histograms when finding bounds to include 99% of integral, this correspinds to 99.6% for F1 and 99.7% for F2
const float fbound_min = 1.7774, fbound_max = 1.95572, fbound_counts = 113883.0, finclusive_counts = 114315.0;
// the total number of counts we desire in sidebands -- will be updated after creating F1/F2Histograms to equal 99% of total entries
float SB_count = 0.99 * finclusive_counts;
double step_size = 1.0e-6;
bool fisB1, fisB2, fisS1, fisS2, fisF1, fisF2, fisSB1, fisSB2;
bool gisB1, gisB2, gisS1, gisS2, gisF1, gisF2, gisSB1, gisSB2;
const double fit_range_low = 1.6, fit_range_high = 2.1, D0_mass = 1.8648;


const float prob_cutoff = 0.1;

using namespace std;

void FitMassPlotsTF1(TF1 *fitFunction, TH1D *MassPlot)
{
    for (int i = 0; i < 20; i++)
    {
        MassPlot->Fit(fitFunction, "Lq0");
    }
    MassPlot->Fit(fitFunction, "M R0");
}

//int temp1_5bins()
int main()
{
    TFile *outfile = new TFile("siidebandslimits_99_skippingswap.root", "recreate");

    // TString inputfile1 = "/home/awesole/pythia8/TH2F_massPlots.root";
    cout << "working to open file..." << endl;
    // TString inputfile1 = "TH2F_massPlots_withBins.root"; // actual one we need
    //TString inputfile1 = "TH2F_output_withbins_30Sept.root"; // actual one we need
    TString inputfile1 = "data_files/TH2F_skippingswap.root"; // actual one we need
    TFile *inf1 = TFile::Open(inputfile1);
    cout << "File Successfully Opened!" << endl;
    // TFile *infile = TFile::Open("/home/awesole/pythia8/SingleParticleMass.root", "READ");i
    //TFile *inf3 = TFile::Open("TF1_outputs_1k.root");
    TFile *inf3 = TFile::Open("TF1_outputs_skippingswap.root");

    TCanvas *c1 = new TCanvas("c1", "c1", 1000, 500);
    c1->Divide(3, 2);
    c1->cd(1);
    TCanvas *c2 = new TCanvas("c2", "c2", 1000, 500);
    c2->Divide(3, 2);
    c2->cd(1);
    TCanvas *c3 = new TCanvas("c3", "c3", 1500, 1000);
    c3->Divide(3, 2);
    c3->cd(1);
    TCanvas *c4 = new TCanvas("c4", "c4", 1500, 1000);
    c4->Divide(3, 2);
    c4->cd(1);

    

    TNtuple *D_mass_values = new TNtuple("D_mass_values", "D_mass_values", "l1:l2:l3:l4");
    TNtuple *Dbar_mass_values = new TNtuple("Dbar_mass_values", "Dbar_mass_values", "l1:l2:l3:l4");

    for (int j = 6; j < 7; j++)
    {
        c1->cd(j);
        c2->cd(j);
        c3->cd(j);
        c4->cd(j);

        /*
        TH2F *S1S2Mass = nullptr; // Background of D0 candidates
        TH1F *S1Mass = nullptr;   // Background of D0 candidates
        */

        TH1D *SWB1Mass = new TH1D("SWB1Mass", "", 1000, 1.6, 2.1); // Background of D0 candidates
        TH1D *SWB2Mass = new TH1D("SWB2Mass", "", 1000, 1.6, 2.1); // Background of D0 candidates
        TH1D *M1Mass = nullptr; // Background of D0 candidates
        TH1D *M2Mass = nullptr; // Background of D0 candidates
        TH2F *M1M2Mass = dynamic_cast<TH2F *>(inf1->Get(Form("M1M2Mass_%d", j)));
        M1Mass = M1M2Mass->ProjectionX();
        M1Mass->SetMinimum(0.5);
        M2Mass = M1M2Mass->ProjectionY();
        M2Mass->SetMinimum(0.5);

        TH1D *B1Mass = nullptr;  // Background of D0 candidates
        TH1D *B2Mass = nullptr;  // Background of D0 candidates
        TH1D *SW1Mass = nullptr; // Background of D0 candidates
        TH1D *SW2Mass = nullptr; // Background of D0 candidates
        TH1D *S1Mass = nullptr;  // Background of D0 candidates
        TH1D *S2Mass = nullptr;  // Background of D0 candidates
        TH2F *S1S2Mass = dynamic_cast<TH2F *>(inf1->Get(Form("S1S2Mass_%d", j)));
        S1Mass = S1S2Mass->ProjectionX();
        S1Mass->SetMinimum(0.5);
        S2Mass = S1S2Mass->ProjectionY();
        S2Mass->SetMinimum(0.5);

        TF1 *F1 = (TF1 *)inf3->Get(Form("F1_%d", j));
        TF1 *F2 = (TF1 *)inf3->Get(Form("F2_%d", j));

        TH2D *SW1SW2Mass = (TH2D *)inf1->Get(Form("SW1SW2Mass_%d", j)); // Full Range of F1F2 Mass
        TH2D *SW1B2Mass = (TH2D *)inf1->Get(Form("SW1B2Mass_%d", j));   // Full Range of F1F2 Mass
        TH2D *SW1S2Mass = (TH2D *)inf1->Get(Form("SW1S2Mass_%d", j));   // Full Range of F1F2 Mass
        TH2D *S1SW2Mass = (TH2D *)inf1->Get(Form("S1SW2Mass_%d", j));   // Full Range of F1F2 Mass
        TH2D *S1B2Mass = (TH2D *)inf1->Get(Form("S1B2Mass_%d", j));     // Full Range of F1F2 Mass
        TH2D *B1SW2Mass = (TH2D *)inf1->Get(Form("B1SW2Mass_%d", j));   // Full Range of F1F2 Mass
        TH2D *B1B2Mass = (TH2D *)inf1->Get(Form("B1B2Mass_%d", j));     // Full Range of F1F2 Mass
        TH2D *B1S2Mass = (TH2D *)inf1->Get(Form("B1S2Mass_%d", j));     // Full Range of F1F2 Mass

        S1Mass->Add(S1B2Mass->ProjectionX());
        S1Mass->Add(S1SW2Mass->ProjectionX());
        S2Mass->Add(B1S2Mass->ProjectionY());
        S2Mass->Add(SW1S2Mass->ProjectionY());

        SW1Mass = SW1SW2Mass->ProjectionX();
        SW1Mass->Add(SW1B2Mass->ProjectionX());
        SW1Mass->Add(SW1S2Mass->ProjectionX());

        SW2Mass = SW1SW2Mass->ProjectionY();
        SW2Mass->Add(B1SW2Mass->ProjectionY());
        SW2Mass->Add(S1SW2Mass->ProjectionY());

        B1Mass = B1B2Mass->ProjectionX();
        B1Mass->Add(B1S2Mass->ProjectionX());
        B1Mass->Add(B1SW2Mass->ProjectionX());

        B2Mass = B1B2Mass->ProjectionY();
        B2Mass->Add(S1B2Mass->ProjectionY());
        B2Mass->Add(SW1B2Mass->ProjectionY());

        SWB1Mass->Add(SW1Mass);
        SWB1Mass->Add(B1Mass);
        SWB2Mass->Add(SW2Mass);
        SWB2Mass->Add(B2Mass);

        TNtuple *Mass_nt = dynamic_cast<TNtuple *>(inf1->Get(Form("MassTuple_%d", j)));

        TF1 *signal1 = new TF1("signal1", "([6]*[7]*([0]*TMath::Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2]) + [5]*((1-[0]))*TMath::Gaus(x,[1],[3])/(sqrt(2*3.14159)*[3]) + (1-[5]) * (1-[0])*TMath::Gaus(x,[1],[4])/(sqrt(2*3.14159)*[4])))", 1.7, 2.1);
        signal1->FixParameter(0, F1->GetParameter(4));  // ratio btw 1 and 2
        signal1->FixParameter(1, F1->GetParameter(1));  // mean
        signal1->FixParameter(2, F1->GetParameter(2));  // sigma 1
        signal1->FixParameter(3, F1->GetParameter(3));  // sihgma 2
        signal1->FixParameter(4, F1->GetParameter(14)); // sigma 3
        signal1->FixParameter(5, F1->GetParameter(12)); // ratio btw gaus 2 & 3
        signal1->FixParameter(6, F1->GetParameter(0));  // scaling
        signal1->FixParameter(7, F1->GetParameter(5));  // signal fraction

        TF1 *signal2 = new TF1("signal2", "([6]*[7]*([0]*TMath::Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2]) + [5]*((1-[0]))*TMath::Gaus(x,[1],[3])/(sqrt(2*3.14159)*[3]) + (1-[5]) * (1-[0])*TMath::Gaus(x,[1],[4])/(sqrt(2*3.14159)*[4])))", 1.7, 2.1);
        signal2->FixParameter(0, F2->GetParameter(4));  // ratio btw 1 and 2
        signal2->FixParameter(1, F2->GetParameter(1));  // mean
        signal2->FixParameter(2, F2->GetParameter(2));  // sigma 1
        signal2->FixParameter(3, F2->GetParameter(3));  // sihgma 2
        signal2->FixParameter(4, F2->GetParameter(14)); // sigma 3
        signal2->FixParameter(5, F2->GetParameter(12)); // ratio btw gaus 2 & 3
        signal2->FixParameter(6, F2->GetParameter(0));  // scaling
        signal2->FixParameter(7, F2->GetParameter(5));  // signal fraction


        TF1 *background1 = new TF1("background1", "[8] + [9]*x + [10]*x*x + [11]*x*x*x", fit_range_low, fit_range_high);
        background1->FixParameter(8, F1->GetParameter(8));
        background1->FixParameter(9, F1->GetParameter(9));
        background1->FixParameter(10, F1->GetParameter(10));
        background1->FixParameter(11, F1->GetParameter(11));
        // B1SW2Mass->ProjectionX()->Draw();
        // background->Draw("same");

        TF1 *swap1 = new TF1("swap1", "[0]*(1-[5])*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6])) + (1-[4])*(TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6]))))", fit_range_low, fit_range_high);
        swap1->SetLineColor(kRed);
        swap1->FixParameter(0, F1->GetParameter(0));  // norm
        swap1->FixParameter(1, F1->GetParameter(1));  // mean mass
        swap1->FixParameter(2, F1->GetParameter(7));  // sigma1
        swap1->FixParameter(3, F1->GetParameter(15)); // sigma2
        swap1->FixParameter(4, F1->GetParameter(13)); // ratio
        swap1->FixParameter(5, F1->GetParameter(5));  // signal fraction
        swap1->FixParameter(6, F1->GetParameter(6));  // smearing

        
        TF1 *swap2 = new TF1("swap2", "[0]*(1-[5])*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6])) + (1-[4])*(TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6]))))", fit_range_low, fit_range_high);
        swap2->SetLineColor(kRed);
        swap2->FixParameter(0, F2->GetParameter(0));  // norm
        swap2->FixParameter(1, F2->GetParameter(1));  // mean mass
        swap2->FixParameter(2, F2->GetParameter(7));  // sigma1
        swap2->FixParameter(3, F2->GetParameter(15)); // sigma2
        swap2->FixParameter(4, F2->GetParameter(13)); // ratio
        swap2->FixParameter(5, F2->GetParameter(5));  // signal fraction
        swap2->FixParameter(6, F2->GetParameter(6));  // smearing

        
        // swap->Draw("same");
        // c1->cd(3);

        TF1 *background2 = new TF1("background2", "[8] + [9]*x + [10]*x*x + [11]*x*x*x", fit_range_low, fit_range_high);
        background2->FixParameter(8, F2->GetParameter(8));
        background2->FixParameter(9, F2->GetParameter(9));
        background2->FixParameter(10, F2->GetParameter(10));
        background2->FixParameter(11, F2->GetParameter(11));
        // B1SW2Mass->ProjectionX()->Draw();
        // background->Draw("same");

        TF1 *SwapBkg1 = new TF1("SwapBkg1", " ([0]*(1-[5])*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6])) + (1-[4])*(TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6]))))) + ([7] + [8]*x + [9]*x*x + [10]*x*x*x) ", fit_range_low, fit_range_high);
        SwapBkg1->FixParameter(0, F1->GetParameter(0));
        SwapBkg1->FixParameter(1, F1->GetParameter(1));
        SwapBkg1->FixParameter(2, F1->GetParameter(7));
        SwapBkg1->FixParameter(3, F1->GetParameter(15));
        SwapBkg1->FixParameter(4, F1->GetParameter(13));
        SwapBkg1->FixParameter(5, F1->GetParameter(5));
        SwapBkg1->FixParameter(6, F1->GetParameter(6));
        SwapBkg1->FixParameter(7, F1->GetParameter(8));
        SwapBkg1->FixParameter(8, F1->GetParameter(9));
        SwapBkg1->FixParameter(9, F1->GetParameter(10));
        SwapBkg1->FixParameter(10, F1->GetParameter(11));

        TF1 *SwapBkg2 = new TF1("SwapBkg2", " ([0]*(1-[5])*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6])) + (1-[4])*(TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6]))))) + ([7] + [8]*x + [9]*x*x + [10]*x*x*x) ", fit_range_low, fit_range_high);
        SwapBkg2->FixParameter(0, F2->GetParameter(0));
        SwapBkg2->FixParameter(1, F2->GetParameter(1));
        SwapBkg2->FixParameter(2, F2->GetParameter(7));
        SwapBkg2->FixParameter(3, F2->GetParameter(15));
        SwapBkg2->FixParameter(4, F2->GetParameter(13));
        SwapBkg2->FixParameter(5, F2->GetParameter(5));
        SwapBkg2->FixParameter(6, F2->GetParameter(6));
        SwapBkg2->FixParameter(7, F2->GetParameter(8));
        SwapBkg2->FixParameter(8, F2->GetParameter(9));
        SwapBkg2->FixParameter(9, F2->GetParameter(10));
        SwapBkg2->FixParameter(10, F2->GetParameter(11));




       // signal1->Draw("same");
        //M1Mass->GetXaxis()->SetRangeUser(1.6,2.1);
       // M1Mass->Add(signal1, -1);
       // M2Mass->Add(signal2, -1);
        //signal2->Draw("same");

        //double bin_width = 1.0e-06;
         double bin_width = M1Mass->GetBinWidth(10);
         

        double count_fg_D, count_fg_Dbar;
        double integralD = signal1->Integral(fit_range_low, fit_range_high);
        cout << "integralD = " << integralD << endl;
        double sb_l2D = signal1->Mean(fit_range_low, fit_range_high);
        double sb_u1D = sb_l2D;
        for (int i = 0; i < 1.0e5; i++)
        {
            count_fg_D = signal1->Integral(sb_l2D, sb_u1D);
            double fraction = count_fg_D / integralD;
            if (fraction > 0.99)
            {
                cout << fixed << setprecision(6) << "determined 99% D0 range, for bin number " << j << endl;
                cout << "mass values are " << sb_l2D << ", " << sb_u1D << endl;
                break;
            }
            else
            {
                sb_l2D -= step_size;
                sb_u1D += step_size;
            }
        }
        double m1D = sb_l2D;
        double m2D = sb_u1D;
        double NB1 = SwapBkg1->Integral(m1D, m2D)/ 0.0005;
        cout << "_______________________ NB1 counts = " << NB1 << endl;

        //M1Mass->GetXaxis()->SetRangeUser(m1D, m2D);
        //double NB1 = M1Mass->Integral();
        for (int i = 0; i < 1.0e9; i++)
        {

            //double count_SB_left = SwapBkg1->Integral(m1D, sb_l2D);
            //double count_SB_left = Mass_nt->GetEntries("M1 > m1D && M1 < sb_l2D");
            double count_SB_left = Mass_nt->GetEntries(Form("M1 >= %f && M1 < %f", m1D, sb_l2D));
            double count_SB_right = Mass_nt->GetEntries(Form("M1 > %f && M1 <= %f", sb_u1D, m2D));
            //M1Mass->GetXaxis()->SetRangeUser(m1D , sb_l2D);
            //double count_SB_left = M1Mass->Integral();
            //double count_SB_right = SwapBkg1->Integral(sb_u1D, m2D);
            //M1Mass->GetXaxis()->SetRangeUser(sb_u1D, m2D);
            //double count_SB_right = M1Mass->Integral();

            double count_SB = count_SB_left + count_SB_right;
            float fraction = count_SB / NB1;

            if (fraction > 1.0)
            {
                cout << fixed << setprecision(6) << "greater than 1.0! D0 limits are " << m1D << " and " << m2D << endl;
                cout << "fraction =" << fraction << " and SB counts = " << count_SB << endl;
                cout << "bin = " << j << endl;
                cout << ".............." << endl;
                break;
            }
            else if(fraction < 0.9) {
                m1D -= 1.0e-3;
                m2D += 1.0e-3;

            }
            else if(fraction < 0.998)  {
                m1D -= 1.0e-4;
                m2D += 1.0e-4;

            }
            else
            {
                cout << "fraction = " << fraction << endl;
                m1D -= step_size;
                m2D += step_size;
            }
        }


        double integralDbar = signal2->Integral(fit_range_low, fit_range_high);
        double sb_l2Dbar = signal2->Mean(fit_range_low, fit_range_high);
        double sb_u1Dbar = sb_l2Dbar;

        for (int i = 0; i < 1.0e5; i++)
        {
            count_fg_Dbar = signal2->Integral(sb_l2Dbar, sb_u1Dbar);
            double fractionbar = count_fg_Dbar / integralDbar;
            if (fractionbar > 0.99)
            {
                cout << fixed << setprecision(6) << "determined 99% Dbar range, for bin number " << j << endl;
                cout << "mass values are " << sb_l2Dbar << ", " << sb_u1Dbar << endl;
                break;
            }
            else
            {
                sb_l2Dbar -= step_size;
                sb_u1Dbar += step_size;
            }
        }
        double m1Dbar = sb_l2Dbar;
        double m2Dbar = sb_u1Dbar;

        double NB2 = SwapBkg2->Integral(m1Dbar, m2Dbar) / 0.0005;
        //M2Mass->GetXaxis()->SetRangeUser(m1Dbar, m2Dbar);
        //double NB2 = M2Mass->Integral();
        for (int i = 0; i < 1.0e6; i++)
        {
            //double count_SB_left = SwapBkg2->Integral(m1Dbar, sb_l2Dbar);
            //double count_SB_right = SwapBkg2->Integral(sb_u1Dbar, m2Dbar);
            double count_SB_left = Mass_nt->GetEntries(Form("M2 >= %f && M2 < %f", m1Dbar, sb_l2Dbar));
            double count_SB_right = Mass_nt->GetEntries(Form("M2 > %f && M2 <= %f", sb_u1Dbar, m2Dbar));

            //M2Mass->GetXaxis()->SetRangeUser(m1Dbar, sb_l2Dbar);
            //double count_SB_left = M2Mass->Integral();
            //double count_SB_right = SwapBkg1->Integral(sb_u1D, m2D);
            //M2Mass->GetXaxis()->SetRangeUser(sb_u1Dbar, m2Dbar);
            //double count_SB_right = M2Mass->Integral();
            double count_SB = count_SB_left + count_SB_right;
            float fraction = count_SB / NB2;

            if (fraction > 1.0)
            {
                cout << fixed << setprecision(6) << "greater than 1.0! Dbar limits are " << m1Dbar << " and " << m2Dbar << endl;
                cout << "fraction =" << fraction << " and SB counts = " << count_SB << " and fg counts" << NB2 << endl;
                cout << "bin = " << j << endl;
                cout << "+++++++++++++++++" << endl;
                cout << "+++++++++++++++++" << endl;
                break;
            }
            else if(fraction < 0.9) {
                m1Dbar -= 1.0e-3;
                m2Dbar += 1.0e-3;

            }
            else if(fraction < 0.998)  {
                m1Dbar -= 1.0e-4;
                m2Dbar += 1.0e-4;

            }
            else
            {
                cout << "fraction = " << fraction << endl;
                m1Dbar -= step_size;
                m2Dbar += step_size;
            }
        }
        c1->cd(j);
        M1Mass->SetMinimum(0.5);
        M1Mass->GetXaxis()->SetRangeUser(1.6,2.1);
        M1Mass->Draw();
        B1Mass->Draw("same");
        SW1Mass->Draw("same");
        //SWB1Mass->SetMinimum(0.5);
        //SWB1Mass->Draw();
        SwapBkg1->SetLineColor(kRed);
        SwapBkg1->Draw("same");
        signal1->SetLineColor(kRed);
        signal1->Draw("same");
        S1Mass->Draw("same");
        background1->Draw("same");
        swap1->Draw("same");
        TLine *l1 = new TLine(sb_l2D, 0, sb_l2D, M1Mass->GetMaximum());
        TLine *l2 = new TLine(sb_u1D, 0, sb_u1D, M1Mass->GetMaximum());
        TLine *l3 = new TLine(m1D, 0, m1D, M1Mass->GetMaximum());
        TLine *l4 = new TLine(m2D, 0, m2D, M1Mass->GetMaximum());
        l1->Draw("same");
        l2->Draw("same");
        l3->Draw("same");
        l4->Draw("same");
        c1->Update();
        c2->cd(j);
        M2Mass->GetXaxis()->SetRangeUser(1.6,2.1);
        M2Mass->SetMinimum(0.5);
        M2Mass->Draw();
        S2Mass->Draw("same");
        B2Mass->Draw("same");
        SW2Mass->Draw("same");
        SwapBkg2->SetLineColor(kRed);
        SwapBkg2->Draw("same");
        signal2->SetLineColor(kRed);
        signal2->Draw("same");
        background2->Draw("same");
        swap2->Draw("same");
        TLine *l5 = new TLine(sb_l2Dbar, 0, sb_l2Dbar, M2Mass->GetMaximum());
        TLine *l6 = new TLine(sb_u1Dbar, 0, sb_u1Dbar, M2Mass->GetMaximum());
        TLine *l7 = new TLine(m1Dbar, 0, m1Dbar, M2Mass->GetMaximum());
        TLine *l8 = new TLine(m2Dbar, 0, m2Dbar, M2Mass->GetMaximum());
        l5->Draw("same");
        l6->Draw("same");
        l7->Draw("same");
        l8->Draw("same");
        c2->Update();
        c3->cd(j);
        //SWB1Mass->Add(M1Mass, -1);
       // SWB1Mass->Divide(M1Mass);
        SWB1Mass->Divide(SwapBkg1);
        SWB1Mass->Draw();
        //l1->Draw("same");
        l2->SetLineColor(kRed);
        l2->SetLineWidth(2);
        l2->Draw("same");
        l3->SetLineColor(kRed);
        l3->SetLineWidth(2);
        l3->Draw("same");
        //l4->Draw("same");
        c4->cd(j);
        //SWB2Mass->Add(M2Mass, -1);
        //SWB2Mass->Divide(M2Mass);
        SWB2Mass->Divide(SwapBkg2);
        SWB2Mass->Draw();
        //l5->Draw("same");
        l6->SetLineColor(kRed);
        l6->SetLineWidth(2);
        l6->Draw("same");
        l7->SetLineColor(kRed);
        l7->SetLineWidth(2);
        l7->Draw("same");
        //l8->Draw("same");
        D_mass_values->Fill(m1D, sb_l2D, sb_u1D, m2D);
        Dbar_mass_values->Fill(m1Dbar, sb_l2Dbar, sb_u1Dbar, m2Dbar);


    }
    outfile->cd();
    D_mass_values->Write();
    Dbar_mass_values->Write();
    outfile->Close();
    /*
    avg_massD.Write("avg_massD");
    avg_massDbar.Write();
    std_massD.Write();
    std_massDbar.Write();
    sb_l1D.Write();
    sb_l1Dbar.Write();
    sb_u2D.Write();
    sb_u2Dbar.Write();
    outfile->Close();
    */

    return 0;
}
