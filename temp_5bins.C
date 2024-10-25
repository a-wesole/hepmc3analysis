#include <iostream>
#include <random>
#include "TLorentzVector.h"
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
#include "TH1F.h"
#include "THStack.h"
#include <TF1.h>

const char *file_path = "data_files/mergedfile29Sept2024_cuts.root";
// const char *file_path = "/scratch/bell/awesole/test_57.root";
const float upper_limit = 3.14159, lower_limit = -3.14159;                 // bounds to transition delta phi from (-2pi, 2pi) to (-0.5pi, 1.5*pi)
const float nbins_mass = 1.0e+04, xmin_mass = 1.6, xmax_mass = 2.1;        // for mass plots
const float nbins = 5.0, xmin = -0.5*TMath::Pi(), xmax = 1.5*TMath::Pi(), ymin = 0.0, ymax = 1.5; // for dphi plots
//const double avg_mass = 1.865, std_mass = 1.902e-02;                       // July 30 updates
//const double avg_mass = 1.86502, std_mass = 0.0147844;                       // July 30 updates
// below are the limits for sidebands...
const float pT_minimum = 0.5; //   GeV/c
/* --99.9%
const float avg_massD[6] = {1.782405, 1.779866, 1.778457, 1.777558, 1.779418, 1.779357};
const float std_massD[6] = {1.947331,  1.949854, 1.951259, 1.952100, 1.950194, 1.950333};
const float sb_l1D[6] = {1.698099, 1.692772, 1.689998, 1.688080, 1.691948, 1.691802}; //99.7 closure   fraction=99.7 data style
const float sb_u2D[6] = {2.031637, 2.036948, 2.039718, 2.041578, 2.037664, 2.037888}; //99.7% closure  fraction=99.7 data style
                                                                                                              
const float avg_massDbar[6] = {1.780648, 1.777830, 1.779124, 1.776793, 1.779139, 1.778525};
const float std_massDbar[6] = {1.948866,  1.951848, 1.950762, 1.952873, 1.950259, 1.951103};
const float sb_l1Dbar[6] = {1.694535, 1.688560, 1.691286, 1.686491, 1.691521, 1.690105}; //99.7 closure   fraction=99.7 data style
const float sb_u2Dbar[6] = {2.034979, 2.041118, 2.038600, 2.043175, 2.037877, 2.039523}; //99.7% closure  fraction=99.7 data style
*/
//-99%
/*
const float avg_massD[6] = {1.806194, 1.804652, 1.804042, 1.802790, 1.804371, 1.804302};
const float std_massD[6] = {1.923542,  1.925068, 1.925674, 1.926868, 1.925241, 1.925388};
const float sb_l1D[6] = {1.746770, 1.743594, 1.742406, 1.739845, 1.743094, 1.742922}; //99.7 closure   fraction=99.7 data style
const float sb_u2D[6] = {1.982966, 1.986126, 1.987310, 1.989813, 1.986518, 1.986768}; //99.7% closure  fraction=99.7 data style
                                                                                                              
const float avg_massDbar[6] = {1.805137, 1.803469, 1.804360, 1.802657, 1.804201, 1.803860};
const float std_massDbar[6] = {1.924377,  1.926209, 1.925526, 1.927009, 1.925197, 1.925768};
const float sb_l1Dbar[6] = {1.744701, 1.741169, 1.742998, 1.739598, 1.742906, 1.742061}; //99.7 closure   fraction=99.7 data style
const float sb_u2Dbar[6] = {1.984813, 1.988509, 1.986888, 1.990068, 1.986492, 1.987567}; //99.7% closure  fraction=99.7 data style
*/

/*--good ones 1k bins 99% 
const float avg_massD[6] = {1.808895, 1.807549, 1.807353, 1.805478, 1.807436, 1.804302};
const float std_massD[6] = {1.920853,  1.922193, 1.922387, 1.924208, 1.922238, 1.925388};
const float sb_l1D[6] = {1.752501, 1.750000, 1.749740, 1.745500, 1.749500, 1.742922}; //99.7 closure   fraction=99.7 data style
const float sb_u2D[6] = {1.977247, 1.979742, 1.980000, 1.984186, 1.980174, 1.986768}; //99.7% closure  fraction=99.7 data style
                                                                                                              
const float avg_massDbar[6] = {1.808221, 1.806867, 1.807167, 1.805911, 1.806857, 1.803860};
const float std_massDbar[6] = {1.921319,  1.922819, 1.922755, 1.923785, 1.922571, 1.925768};
const float sb_l1Dbar[6] = {1.751501, 1.748500, 1.749001, 1.746500, 1.748428, 1.742061}; //99.7 closure   fraction=99.7 data style
const float sb_u2Dbar[6] = {1.978039, 1.981186, 1.980921, 1.983196, 1.981000, 1.987567}; //99.7% closure  fraction=99.7 data style
*/
//4k bins 99%
const float avg_massD[6] = {1.806194, 1.807549, 1.807353, 1.802790, 1.804371, 1.804302};
const float std_massD[6] = {1.923542,  1.922193, 1.922387, 1.926868, 1.925241, 1.925388};
const float sb_l1D[6] = {1.733000, 1.750000, 1.749501, 1.739845, 1.743094, 1.742922}; //99.7 closure   fraction=99.7 data style
const float sb_u2D[6] = {1.996736, 1.979742, 1.980239, 1.989813, 1.986518, 1.986768}; //99.7% closure  fraction=99.7 data style
                                                                                                              
const float avg_massDbar[6] = {1.805137, 1.806867, 1.807167, 1.802657, 1.804201, 1.803860};
const float std_massDbar[6] = {1.924377,  1.922819, 1.922755, 1.927009, 1.925197, 1.925768};
const float sb_l1Dbar[6] = {1.731014, 1.748500, 1.749001, 1.739598, 1.742906, 1.742061}; //99.7 closure   fraction=99.7 data style
const float sb_u2Dbar[6] = {1.998500, 1.981186, 1.980921, 1.990068, 1.986492, 1.987567}; //99.7% closure  fraction=99.7 data style

const float in_val = -0.5 * TMath::Pi(), sc_val = (2 * TMath::Pi() / 5);
const float phi_array[6] = {in_val, in_val + sc_val, in_val + 2 * sc_val, in_val + 3 * sc_val, in_val + 4 * sc_val, in_val + 5 * sc_val};

int bin = 0;
int goodbin = 1;

double step_size = 0.00001;
bool fisB1, fisB2, fisS1, fisS2, fisF1, fisF2, fisSB1, fisSB2;
bool gisB1, gisB2, gisS1, gisS2, gisF1, gisF2, gisSB1, gisSB2;
bool fisD0candidate, fisDbarcandidate, gisD0candidate, gisDbarcandidate;
bool f_foreground, f_sidebands, g_foreground, g_sidebands;

// const float prob_cutoff = 0.1;
const float prob_cutoff = 0.1;

using namespace std;

float transition_phi(float &D0del_phi)
{
    // this function reads in an angle d0del_phi and changes it so that the range is confined to (-pi/2 to 3pi/2)
    if (D0del_phi > upper_limit)
    {
        D0del_phi = D0del_phi - 2 * TMath::Pi();
    }
    if (D0del_phi < lower_limit)
    {
        D0del_phi = 2 * TMath::Pi() + D0del_phi;
    }
    if (D0del_phi < -0.5 * TMath::Pi() && D0del_phi > -1 * TMath::Pi())
    {
        D0del_phi = 2 * TMath::Pi() + D0del_phi;
    }
    return D0del_phi;
}

int temp_5bins()
{
    TString outfile = TString(Form("output/bin%d_99_data_1kbins.root", goodbin));

    std::vector<float> *mass = nullptr;
    std::vector<float> *phi = nullptr;
    std::vector<float> *prob = nullptr;
    std::vector<float> *vertex_no = nullptr;

    std::vector<float> *iparticle_px = nullptr;
    std::vector<float> *iparticle_py = nullptr;
    std::vector<float> *iparticle_pz = nullptr;

    std::vector<float> *jparticle_px = nullptr;
    std::vector<float> *jparticle_py = nullptr;
    std::vector<float> *jparticle_pz = nullptr;
    std::vector<float> *phy_process = nullptr;
    std::vector<float> *TruePdg_i = nullptr;
    std::vector<float> *TruePdg_j = nullptr;
    std::vector<float> *AssignedPdg_i = nullptr;
    std::vector<float> *AssignedPdg_j = nullptr;
    std::vector<float> *fromD0 = nullptr;
    // std::vector<float> F1Mass, F2Mass, SB1Mass, SB2Mass, S1Mass, S2Mass, FullRangeF1, FullRangeF2, FullRangeF1Mass, FullRangeF2Mass;

    /*
    TH1F *inclusiveSignal = new TH1F("inclusiveSignal", "inclusive signal dphi", nbins, xmin, xmax);
    TH1F *correlated_dphi = new TH1F("correlated_dphi", "correlated dphi", nbins, xmin, xmax);
    TH1F *uncorrelated_dphi = new TH1F("uncorrelated_dphi", "uncorrelated dphi", nbins, xmin, xmax);
    TH1F *background_dphi = new TH1F("background_dphi", "background dphi", nbins, xmin, xmax);
    TH1F *Gluon_splitting = new TH1F("Gluon_splitting", "gluon splitting", nbins, xmin, xmax);
    TH1F *Gluon_fusion = new TH1F("Gluon_fusion", "gluon fusion", nbins, xmin, xmax);
    TH1F *flavor_excitation_quark = new TH1F("flavor_excitation_quark", "flavor_excitation_quark ", nbins, xmin, xmax);
    TH1F *flavor_excitation_gluon = new TH1F("flavor_excitation_gluon", "flavor_excitation_gluon ", nbins, xmin, xmax);
    TH1F *Quark_annihlation = new TH1F("Quark_annihlation", "Quark_annihlation ", nbins, xmin, xmax);
    TH1F *nonprompt_phi = new TH1F("nonprompt_dphi", "nonprompt dphi", nbins, xmin, xmax);
    */

    std::vector<TH1F *> inclusiveSignal;
    std::vector<TH1F *> F1Mass;
    std::vector<TH1F *> F2Mass;
    std::vector<TH1F *> F1F2;
    std::vector<TH1F *> F1SB2;
    std::vector<TH1F *> SB1F2;
    std::vector<TH1F *> SB1SB2;
    std::vector<TH1F *> B1B2;
    std::vector<TH1F *> F1B2;
    std::vector<TH1F *> B1F2;
    std::vector<TH1F *> S1S2;
    std::vector<TH1F *> S1B2;
    std::vector<TH1F *> B1S2;
    std::vector<TH2F *> FullRangeM1M2Mass;
    std::vector<TH2F *> S1S2inclusiveSignalMass;

    //for (int i = 1; i < 2; i++)
    //{
    int b=1;

        TString name1 = TString::Format("F1_%d", b);
        TString name2 = TString::Format("F2_%d", b);
        TString name3 = TString::Format("F1F2_%d", b);
        TString name4 = TString::Format("F1SB2_%d", b);
        TString name5 = TString::Format("SB1F2_%d", b);
        TString name6 = TString::Format("SB1SB2_%d", b);
        TString name7 = TString::Format("B1B2_%d", b);
        TString name8 = TString::Format("F1B2_%d", b);
        TString name9 = TString::Format("B1F2_%d", b);
        TString name10 = TString::Format("S1S2_%d", b);
        TString name11 = TString::Format("S1B2_%d", b);
        TString name12 = TString::Format("FB1S2_%d", b);
        TString name13 = TString::Format("FullRangeM1M2_%d", b);
        TString name14 = TString::Format("S1S2InclusiveSignal_%d", b);
        TString name15 = TString::Format("InclusiveSignal_%d", b);

        TH1F *F1 = new TH1F(name1, name1, nbins_mass, xmin_mass, xmax_mass); // Foreground of D0 candidates
        TH1F *F2 = new TH1F(name2, name2, nbins_mass, xmin_mass, xmax_mass); // foreground D0bar candidates
        TH1F *inclusiveSignalH = new TH1F(name15, name15, nbins, xmin, xmax);
        TH1F *B1Mass = new TH1F("B1Mass", "B1 Mass", nbins_mass, xmin_mass, xmax_mass);    // Background of D0 candidates
        TH1F *B2Mass = new TH1F("B2Mass", "B2 Mass", nbins_mass, xmin_mass, xmax_mass);    // Background D0bar candidates
        TH1F *S1Mass = new TH1F("S1Mass", "S1 Mass", nbins_mass, xmin_mass, xmax_mass);    // D0 signal
        TH1F *S2Mass = new TH1F("S2Mass", "S2 Mass", nbins_mass, xmin_mass, xmax_mass);    // D0bar signal
        //TH1F *SW1Mass = new TH1F("SW1Mass", "SW1 Mass", nbins_mass, xmin_mass, xmax_mass); // SW1 mass in 3 sigma
        //TH1F *SW2Mass = new TH1F("SW2Mass", "SW2 Mass", nbins_mass, xmin_mass, xmax_mass); // SW2 mass in 3 sigma
        // TH1F *SBSW1Mass = new TH1F("SBSW1Mass", "SBSW1Mass", nbins_mass, xmin_mass, xmax_mass);
        // TH1F *SBSW2Mass = new TH1F("SBSW2Mass", "SBSW2Mass", nbins_mass, xmin_mass, xmax_mass);
        TH1F *SB1Mass = new TH1F("SB1Mass", "SB1Mass", nbins_mass, xmin_mass, xmax_mass);
        TH1F *SB2Mass = new TH1F("SB2Mass", "SB2Mass", nbins_mass, xmin_mass, xmax_mass);

        TH1F *F1F2H = new TH1F(name3, name3, nbins, xmin, xmax);
        // F1F2->Sumw2();
        TH1F *F1SB2H = new TH1F(name4, name4, nbins, xmin, xmax);
        // F1SB2->Sumw2();
        TH1F *SB1F2H = new TH1F(name5, name5, nbins, xmin, xmax);
        // SB1F2->Sumw2();
        TH1F *SB1SB2H = new TH1F(name6, name6, nbins, xmin, xmax);
        // SB1SB2->Sumw2();
        // TH1F *S1S2 = new TH1F("S1S2", "S1S2", nbins, xmin, xmax);
        TH1F *B1B2H = new TH1F(name7, name7, nbins, xmin, xmax);
        // B1B2->Sumw2();
        TH1F *F1B2H = new TH1F(name8, name8, nbins, xmin, xmax);
        // F1B2->Sumw2();
        TH1F *B1F2H = new TH1F(name9, name9, nbins, xmin, xmax);
        // B1F2->Sumw2();

        TH1F *S1S2H = new TH1F(name10, name10, nbins, xmin, xmax);
        // S1S2->Sumw2();
        TH1F *S1B2H = new TH1F(name11, name11, nbins, xmin, xmax);
        // S1B2->Sumw2();
        TH1F *B1S2H = new TH1F(name12, name12, nbins, xmin, xmax);
        // B1S2->Sumw2();

        TH2F *FullRangeM1M2 = new TH2F(name13, name13, 100.0, 1.6, 2.1, 100.0, 1.6, 2.1);
        FullRangeM1M2->SetXTitle("M1 Mass");
        FullRangeM1M2->SetYTitle("M2 Mass");
        FullRangeM1M2->SetOption("SURF1");

        TH2F *S1S2inclusiveSignal = new TH2F(name14, name14, 16.0, 1.6, 2.1, 16.0, 1.6, 2.1);
        S1S2inclusiveSignal->SetXTitle("M1 Mass");
        S1S2inclusiveSignal->SetYTitle("M2 Mass");
        S1S2inclusiveSignal->SetOption("SURF1");

        F1Mass.push_back(F1);
        F2Mass.push_back(F2);
        inclusiveSignal.push_back(inclusiveSignalH);

        F1F2.push_back(F1F2H);
        F1SB2.push_back(F1SB2H);
        SB1F2.push_back(SB1F2H);
        SB1SB2.push_back(SB1SB2H);
        B1B2.push_back(B1B2H);
        F1B2.push_back(F1B2H);
        B1F2.push_back(B1F2H);
        S1S2.push_back(S1S2H);
        S1B2.push_back(S1B2H);
        B1S2.push_back(B1S2H);
        FullRangeM1M2Mass.push_back(FullRangeM1M2);
        S1S2inclusiveSignalMass.push_back(S1S2inclusiveSignal);
   // }

    bool debug = false;
    float Ifile, event_no, D0_index;
    float pT_i, pT_j;

    TFile *infile = new TFile(file_path);
    TTree *t = (TTree *)infile->Get("event_tree");

    t->SetBranchAddress("mass", &mass);
    t->SetBranchAddress("phi", &phi);
    t->SetBranchAddress("prob", &prob);
    t->SetBranchAddress("vertex_no", &vertex_no);
    t->SetBranchAddress("iparticle_px", &iparticle_px);
    t->SetBranchAddress("iparticle_py", &iparticle_py);
    t->SetBranchAddress("iparticle_pz", &iparticle_pz);
    t->SetBranchAddress("jparticle_px", &jparticle_px);
    t->SetBranchAddress("jparticle_py", &jparticle_py);
    t->SetBranchAddress("jparticle_pz", &jparticle_pz);
    t->SetBranchAddress("phy_process", &phy_process);
    t->SetBranchAddress("TruePdg_i", &TruePdg_i);
    t->SetBranchAddress("TruePdg_j", &TruePdg_j);
    t->SetBranchAddress("AssignedPdg_i", &AssignedPdg_i);
    t->SetBranchAddress("AssignedPdg_j", &AssignedPdg_j);
    t->SetBranchAddress("event_no", &event_no);
    t->SetBranchAddress("Ifile", &Ifile);
    t->SetBranchAddress("fromD0", &fromD0);

    for (int i = 0; i < t->GetEntries(); i++)
    //for (int i = 0; i < 20000; i++)
    {
        t->GetEntry(i);
        // if(Ifile!=42) continue;
        // if(event_no!=8) continue;
        // cout << "Ifile =" << Ifile << " and event_no =" << event_no << endl;
        if (i % 6000 == 0)
            cout << i << " / " << t->GetEntries() << "  " << 100 * i / t->GetEntries() << "%" << endl;

        for (int f = 0; f < phy_process->size() - 1; f++)
        { // for first kpi pair

            fisB1 = fisB2 = fisS1 = fisS2 = fisF1 = fisF2 = fisSB1 = fisSB2 = fisD0candidate = fisDbarcandidate = f_foreground = f_sidebands = false;

            pT_i = 0.0;
            pT_j = 0.0;
            if (prob->at(f) > prob_cutoff && fromD0->at(f) == 0)
                continue; // background cuts
            if (phy_process->at(f) > 5)
                continue; // skip all phy_process that are not 1-5

            pT_i = std::sqrt(iparticle_px->at(f) * iparticle_px->at(f) + iparticle_py->at(f) * iparticle_py->at(f));
            pT_j = std::sqrt(jparticle_px->at(f) * jparticle_px->at(f) + jparticle_py->at(f) * jparticle_py->at(f));
            if (pT_i < pT_minimum || pT_j < pT_minimum)
                continue; // apply pT cuts tp each daughter (iparticle and j particle)

            // if(fromD0->at(f)!=0 && AssignedPdg_i->at(f) == TruePdg_i->at(f) && AssignedPdg_j->at(f) == TruePdg_j->at(f)) continue;//skip all signal
            // if(fromD0->at(f)!=0 && AssignedPdg_i->at(f) == -TruePdg_j->at(f) && AssignedPdg_j->at(f) == -TruePdg_i->at(f)) continue;//skip all swap
            // if(fromD0->at(f) ==0) continue;

            // determine D0 or Dbar candidate
            if ((AssignedPdg_i->at(f) == 211 && AssignedPdg_j->at(f) == -321) || (AssignedPdg_i->at(f) == -321 && AssignedPdg_j->at(f) == 211))
                fisD0candidate = true;
            if ((AssignedPdg_i->at(f) == -211 && AssignedPdg_j->at(f) == 321) || (AssignedPdg_i->at(f) == 321 && AssignedPdg_j->at(f) == -211))
                fisDbarcandidate = true;

            ///////////////////// proceed to second kpi pair //////////////////////////////////////////////////////////////////////////////

            for (int g = f + 1; g < phy_process->size(); g++)
            { // for second kpi pair

                bin = 0;

                gisB1 = gisB2 = gisS1 = gisS2 = gisF1 = gisF2 = gisSB1 = gisSB2 = gisD0candidate = gisDbarcandidate = g_foreground = g_sidebands = false;

                pT_i = 0.0;
                pT_j = 0.0;
                if (prob->at(g) > prob_cutoff && fromD0->at(g) == 0)
                    continue; // background cuts
                if (phy_process->at(g) > 5)
                    continue; // skip all phy_process that are not 1-5

                // if(fromD0->at(g)!=0 && AssignedPdg_i->at(g) == TruePdg_i->at(g) && AssignedPdg_j->at(g) == TruePdg_j->at(g)) continue;//skip all signal
                // if(fromD0->at(g)!=0 && AssignedPdg_i->at(g) == -TruePdg_j->at(g) && AssignedPdg_j->at(g) == -TruePdg_i->at(g)) continue;//skip all swap
                // if(fromD0->at(g) ==0) continue;

                pT_i = std::sqrt(iparticle_px->at(g) * iparticle_px->at(g) + iparticle_py->at(g) * iparticle_py->at(g));
                pT_j = std::sqrt(jparticle_px->at(g) * jparticle_px->at(g) + jparticle_py->at(g) * jparticle_py->at(g));
                if (pT_i < pT_minimum || pT_j < pT_minimum)
                    continue; // apply pT cuts to each daughter (iparticle & jparticle)

                float DDbar_phi = phi->at(f) - phi->at(g); // store the dphi if hte kpi pair are from D0candidate then D0bar candidate
                float DbarD_phi = phi->at(g) - phi->at(f); // store the dphi if hte kpi pair are from D0 bar candidate then D0candidate

                transition_phi(DDbar_phi); // transition the phi to be within (-0.5pi, 1.5pi)
                if (abs(DDbar_phi) < 1.0e-6)
                    continue;
                transition_phi(DbarD_phi);
                if (abs(DbarD_phi) < 1.0e-6)
                    continue;

                // determine D0 and Dbar candidate:
                if ((AssignedPdg_i->at(g) == 211 && AssignedPdg_j->at(g) == -321) || (AssignedPdg_i->at(g) == -321 && AssignedPdg_j->at(g) == 211))
                    gisD0candidate = true;
                if ((AssignedPdg_i->at(g) == -211 && AssignedPdg_j->at(g) == 321) || (AssignedPdg_i->at(g) == 321 && AssignedPdg_j->at(g) == -211))
                    gisDbarcandidate = true;

                /////////////////////////////////////////////////////////////////////////////
                // begin to classify the DDbar phi
                //  fromD0==0 : kpi does not come from D0
                //  fromD0 !=0 : kpi comes from D0
                // vertex_no ==0 : not signal
                // vertex_no !=0 : true signal
                /////////////////////////////////////////////////////////////////////////////

                if (fisD0candidate && gisDbarcandidate)
                { // d0 then d0bar
                  // D0candidate then D0barcandidate


                    bin=0;
                    if (DDbar_phi >= phi_array[0] && DDbar_phi < phi_array[1])
                        bin = 1;
                    if (DDbar_phi >= phi_array[1] && DDbar_phi < phi_array[2])
                        bin = 2;
                    if (DDbar_phi >= phi_array[2] && DDbar_phi < phi_array[3])
                        bin = 3;
                    if (DDbar_phi >= phi_array[3] && DDbar_phi < phi_array[4])
                        bin = 4;
                    if (DDbar_phi >= phi_array[4] && DDbar_phi <= phi_array[5])
                        bin = 5;
                    if (bin != goodbin) continue;
                    if (bin != goodbin) cout << "Error BREAK  fD0 && bin = " << bin << endl;
                    //cout << "f d  bin=" << bin << endl; 
                    /*
                    if (DDbar_phi < phi_array[3] && bin !=0) {
                          cout << "uh oh!  bin should=0 but instead bin = " << bin << endl;
                    }
                    if(bin ==0) continue; 
                    */

                    // determine if f is in foreground or sidebands ( mass regions)
                    //cout << "look here  avg_massD[bin-1] = " << avg_massD[bin-1] << endl;
                    if (mass->at(f) > avg_massD[bin-1] && mass->at(f) < std_massD[bin-1])
                        f_foreground = true;
                    else if ((mass->at(f) <= avg_massD[bin -1] && mass->at(f) >= sb_l1D[bin-1]) || (mass->at(f) >= std_massD[bin-1] && mass->at(f) <= sb_u2D[bin-1]))
                        f_sidebands = true;

                    /////////////////////determine background, swap, signal and sidebands components for FIRST kpi pair////////////////////////////

                    if (f_foreground)
                    { // mass is within foreground region (3 sigma)
                        if (fisD0candidate)
                            fisF1 = true;
                        if (fisDbarcandidate)
                            fisF2 = true;
                        if (fromD0->at(f) == 0) // fromD0==0 means kpi do not go to D0
                        {
                            if (fisD0candidate)
                                fisB1 = true;
                            if (fisDbarcandidate)
                                fisB2 = true;
                        }
                        if (fromD0->at(f) != 0 && AssignedPdg_i->at(f) == TruePdg_i->at(f) && AssignedPdg_j->at(f) == TruePdg_j->at(f)) // True Signal
                        {
                            if (fisD0candidate)
                                fisS1 = true;
                            if (fisDbarcandidate)
                                fisS2 = true;
                        }
                        if (fromD0->at(f) != 0 && AssignedPdg_i->at(f) == -TruePdg_j->at(f) && AssignedPdg_j->at(f) == -TruePdg_i->at(f)) // Swap
                        {
                            if (fisD0candidate)
                                fisB1 = true;
                            if (fisDbarcandidate)
                                fisB2 = true;
                        }
                    }

                    else if (f_sidebands)
                    {                                                                                                                     // mass in within sidebands region
                        if (fisD0candidate)
                            fisSB1 = true;
                        if (fisDbarcandidate)
                            fisSB2 = true;
                        /*
                        if (fromD0->at(f) != 0 && AssignedPdg_i->at(f) == -TruePdg_j->at(f) && AssignedPdg_j->at(f) == -TruePdg_i->at(f)) // Sidebands - swap
                        {
                            if (fisD0candidate)
                                fisSB1 = true;
                            if (fisDbarcandidate)
                                fisSB2 = true;
                        }
                        if (fromD0->at(f) == 0) // sidebands - background
                        {
                            if (fisD0candidate)
                                fisSB1 = true;
                            if (fisDbarcandidate)
                                fisSB2 = true;
                        }
                        */
                    }
                    else continue;
                    // determine if mass is foreground or sidebands region, else continue;
                    if (mass->at(g) > avg_massDbar[bin-1] && mass->at(g) < std_massDbar[bin-1])
                        g_foreground = true;
                    else if ((mass->at(g) <= avg_massDbar[bin - 1] && mass->at(g) >= sb_l1Dbar[bin - 1]) || (mass->at(g) >= std_massDbar[bin - 1] && mass->at(g) <= sb_u2Dbar[bin - 1]))
                            g_sidebands = true;

                    /////////////////////determine background, swap, signal and sidebands components for SECOND kpi pair////////////////////////////

                    if (g_foreground)
                    { // mass is within foreground region (3 sigma)
                        if (gisD0candidate)
                            gisF1 = true;
                        if (gisDbarcandidate)
                            gisF2 = true;
                        if (fromD0->at(g) == 0) // from D0 ==0 means background
                        {
                            if (gisD0candidate)
                                gisB1 = true;
                            if (gisDbarcandidate)
                                gisB2 = true;
                        }
                        if (fromD0->at(g) != 0 && AssignedPdg_i->at(g) == TruePdg_i->at(g) && AssignedPdg_j->at(g) == TruePdg_j->at(g)) // Signal
                        {
                            if (gisD0candidate)
                                gisS1 = true;
                            if (gisDbarcandidate)
                                gisS2 = true;
                        }
                        if (fromD0->at(g) != 0 && AssignedPdg_i->at(g) == -TruePdg_j->at(g) && AssignedPdg_j->at(g) == -TruePdg_i->at(g)) // Swap
                        {
                            if (gisD0candidate)
                                gisB1 = true;
                            if (gisDbarcandidate)
                                gisB2 = true;
                        }
                    }

                    else if (g_sidebands)
                    {                                                                                                                     // mass in within sidebands region
                        if (gisD0candidate)
                            gisSB1 = true;
                        if (gisDbarcandidate)
                            gisSB2 = true;
                        /*
                    if (fromD0->at(g) != 0 && AssignedPdg_i->at(g) == -TruePdg_j->at(g) && AssignedPdg_j->at(g) == -TruePdg_i->at(g)) // SB - swap
                    {
                        if (gisD0candidate)
                            gisSB1 = true;
                        if (gisDbarcandidate)
                            gisSB2 = true;
                    }
                    if (fromD0->at(g) == 0) // SB - background
                    {
                        if (gisD0candidate)
                            gisSB1 = true;
                        if (gisDbarcandidate)
                            gisSB2 = true;
                    }
                    */
                    }
                    else continue;

                    if (fromD0->at(f) != 0 && fromD0->at(g) != 0 && AssignedPdg_i->at(g) == TruePdg_i->at(g) && AssignedPdg_i->at(f) == TruePdg_i->at(f) &&
                        AssignedPdg_j->at(g) == TruePdg_j->at(g) && AssignedPdg_j->at(f) == TruePdg_j->at(f) && f_foreground && g_foreground)
                    { // any kpi pair from DDbar pair

                        inclusiveSignal.at(0)->Fill(DDbar_phi); // any ddbar pair from D0s -- inclusive signal
                        S1S2inclusiveSignalMass.at(0)->Fill(mass->at(f), mass->at(g));

                        /*
                        if (fromD0->at(f) == -101 && fromD0->at(g) == -101)
                            nonprompt_phi->Fill(DDbar_phi);
                        if (((vertex_no->at(f) != vertex_no->at(g)) || (vertex_no->at(f) == 0 && vertex_no->at(g) == 0)) && (fromD0->at(f) != -101 || fromD0->at(g) != -101))
                        {
                            uncorrelated_dphi->Fill(DDbar_phi); // dphi from different origins, unocrrelated
                        }
                        if (vertex_no->at(f) == vertex_no->at(g) && vertex_no->at(f) < 0 && vertex_no->at(g) < 0)
                        {
                            correlated_dphi->Fill(DDbar_phi); // correlated dphi pair, true signal
                            if (phy_process->at(f) == 1)
                                Gluon_splitting->Fill(DDbar_phi);
                            if (phy_process->at(f) == 2)
                                Gluon_fusion->Fill(DDbar_phi);
                            if (phy_process->at(f) == 3)
                                flavor_excitation_gluon->Fill(DDbar_phi);
                            if (phy_process->at(f) == 4)
                                flavor_excitation_quark->Fill(DDbar_phi);
                            if (phy_process->at(f) == 5)
                                Quark_annihlation->Fill(DDbar_phi);
                        }*/
                    }
                    if (mass->at(f) > 1.6 && mass->at(f) < 2.1 && mass->at(g) > 1.6 && mass->at(g) < 2.1)
                    {
                        F1Mass.at(0)->Fill(mass->at(f));
                        F2Mass.at(0)->Fill(mass->at(g));
                        FullRangeM1M2Mass.at(0)->Fill(mass->at(f), mass->at(g));
                    }

                    /*
                    if ((fromD0->at(f) == 0 || fromD0->at(g) == 0))
                        background_dphi.at(bin-1)->Fill(DDbar_phi); // background dphi
                        */

                    // for dphi correlation below!!!!!!!!!!!!!!!!!!!!//

                    if (fisF1 && gisF2)
                        F1F2.at(0)->Fill(DDbar_phi);
                    if (fisF1 && gisSB2)
                        F1SB2.at(0)->Fill(DDbar_phi);
                    if (fisSB1 && gisF2)
                        SB1F2.at(0)->Fill(DDbar_phi);
                    if (fisSB1 && gisSB2)
                    {
                        SB1SB2.at(0)->Fill(DDbar_phi);
                        ////SB1Mass->Fill(mass->at(f));
                        // SB2Mass->Fill(mass->at(g));
                    }
                    if (fisB1 && gisB2)
                    {
                        B1B2.at(0)->Fill(DDbar_phi);
                        // B1Mass->Fill(mass->at(f));
                        // B2Mass->Fill(mass->at(g));
                    }
                    if (fisB1 && gisF2)
                        B1F2.at(0)->Fill(DDbar_phi);
                    if (fisF1 && gisB2)
                        F1B2.at(0)->Fill(DDbar_phi);

                    // if(fisS1 && gisS2) S1S2->Fill(DDbar_phi);
                    if (fisS1 && gisB2)
                        S1B2.at(0)->Fill(DDbar_phi);
                    if (fisB1 && gisS2)
                        B1S2.at(0)->Fill(DDbar_phi);

                } // D0can then D0barcan

                else if (fisDbarcandidate && gisD0candidate)
                { // d0bar then d0
                    // if D0barcand then D0cand


                    bin=0;
                    if (DbarD_phi >= phi_array[0] && DbarD_phi < phi_array[1])
                        bin = 1;
                    if (DbarD_phi >= phi_array[1] && DbarD_phi < phi_array[2])
                        bin = 2;
                    if (DbarD_phi >= phi_array[2] && DbarD_phi < phi_array[3])
                        bin = 3;
                    if (DbarD_phi >= phi_array[3] && DbarD_phi < phi_array[4])
                        bin = 4;
                    if (DbarD_phi >= phi_array[4] && DbarD_phi <= phi_array[5])
                        bin = 5;
                    if (bin != goodbin) continue;
                    if (bin !=goodbin ) cout << "Error BREAK  fDbar && bin = " << bin << endl;
                    //cout << "f dbar  bin=" << bin << endl; 
                    /*
                    if (DbarD_phi < phi_array[3] && bin !=0) {
                          cout << "uh oh!  bin should=0 but instead bin = " << bin << endl;
                    }
                    if(bin ==0) continue; 
                    */

                    // determine if f is in foreground or sidebands ( mass regions)
                    if (mass->at(f) > avg_massDbar[bin-1] && mass->at(f) < std_massDbar[bin-1])
                        f_foreground = true;
                    else if ((mass->at(f) <= avg_massDbar[bin -1] && mass->at(f) >= sb_l1Dbar[bin-1]) || (mass->at(f) >= std_massDbar[bin-1] && mass->at(f) <= sb_u2Dbar[bin-1]))
                        f_sidebands = true;

                    if (f_foreground)
                    { // mass is within foreground region (3 sigma)
                        if (fisD0candidate)
                            fisF1 = true;
                        if (fisDbarcandidate)
                            fisF2 = true;
                        if (fromD0->at(f) == 0) // fromD0==0 means kpi do not go to D0
                        {
                            if (fisD0candidate)
                                fisB1 = true;
                            if (fisDbarcandidate)
                                fisB2 = true;
                        }
                        if (fromD0->at(f) != 0 && AssignedPdg_i->at(f) == TruePdg_i->at(f) && AssignedPdg_j->at(f) == TruePdg_j->at(f)) // True Signal
                        {
                            if (fisD0candidate)
                                fisS1 = true;
                            if (fisDbarcandidate)
                                fisS2 = true;
                        }
                        if (fromD0->at(f) != 0 && AssignedPdg_i->at(f) == -TruePdg_j->at(f) && AssignedPdg_j->at(f) == -TruePdg_i->at(f)) // Swap
                        {
                            if (fisD0candidate)
                                fisB1 = true;
                            if (fisDbarcandidate)
                                fisB2 = true;
                        }
                    }

                    else if (f_sidebands)
                    {                                                                                                                     // mass in within sidebands region
                        if (fisD0candidate)
                            fisSB1 = true;
                        if (fisDbarcandidate)
                            fisSB2 = true;

                        /*
                        if (fromD0->at(f) != 0 && AssignedPdg_i->at(f) == -TruePdg_j->at(f) && AssignedPdg_j->at(f) == -TruePdg_i->at(f)) // Sidebands - swap
                        {
                            if (fisD0candidate)
                                fisSB1 = true;
                            if (fisDbarcandidate)
                                fisSB2 = true;
                        }
                        if (fromD0->at(f) == 0) // sidebands - background
                        {
                            if (fisD0candidate)
                                fisSB1 = true;
                            if (fisDbarcandidate)
                                fisSB2 = true;
                        }
                        */
                    }
                    else continue;
                    // determine if mass is foreground or sidebands region, else continue;
                    if (mass->at(g) > avg_massD[bin-1] && mass->at(g) < std_massD[bin-1])
                        g_foreground = true;
                    else if ((mass->at(g) <= avg_massD[bin -1] && mass->at(g) >= sb_l1D[bin-1]) || (mass->at(g) >= std_massD[bin-1] && mass->at(g) <= sb_u2D[bin-1]))
                        g_sidebands = true;

                    /////////////////////determine background, swap, signal and sidebands components for SECOND kpi pair////////////////////////////

                    if (g_foreground)
                    { // mass is within foreground region (3 sigma)
                        if (gisD0candidate)
                            gisF1 = true;
                        if (gisDbarcandidate)
                            gisF2 = true;
                        if (fromD0->at(g) == 0) // from D0 ==0 means background
                        {
                            if (gisD0candidate)
                                gisB1 = true;
                            if (gisDbarcandidate)
                                gisB2 = true;
                        }
                        if (fromD0->at(g) != 0 && AssignedPdg_i->at(g) == TruePdg_i->at(g) && AssignedPdg_j->at(g) == TruePdg_j->at(g)) // Signal
                        {
                            if (gisD0candidate)
                                gisS1 = true;
                            if (gisDbarcandidate)
                                gisS2 = true;
                        }
                        if (fromD0->at(g) != 0 && AssignedPdg_i->at(g) == -TruePdg_j->at(g) && AssignedPdg_j->at(g) == -TruePdg_i->at(g)) // Swap
                        {
                            if (gisD0candidate)
                                gisB1 = true;
                            if (gisDbarcandidate)
                                gisB2 = true;
                        }
                    }

                    else if (g_sidebands)
                    {                                                                                                                     // mass in within sidebands region
                        if (gisD0candidate)
                            gisSB1 = true;
                        if (gisDbarcandidate)
                            gisSB2 = true;
                        /*
                        if (fromD0->at(g) != 0 && AssignedPdg_i->at(g) == -TruePdg_j->at(g) && AssignedPdg_j->at(g) == -TruePdg_i->at(g)) // SB - swap
                        {
                            if (gisD0candidate)
                                gisSB1 = true;
                            if (gisDbarcandidate)
                                gisSB2 = true;
                        }
                        if (fromD0->at(g) == 0) // SB - background
                        {
                            if (gisD0candidate)
                                gisSB1 = true;
                            if (gisDbarcandidate)
                                gisSB2 = true;
                        }
                        */
                    }
                    else continue;

                    if (fromD0->at(f) != 0 && fromD0->at(g) != 0 && AssignedPdg_i->at(g) == TruePdg_i->at(g) && AssignedPdg_i->at(f) == TruePdg_i->at(f) &&
                        AssignedPdg_j->at(g) == TruePdg_j->at(g) && AssignedPdg_j->at(f) == TruePdg_j->at(f) && f_foreground && g_foreground)
                    {
                        inclusiveSignal.at(0)->Fill(DbarD_phi);
                        S1S2inclusiveSignalMass.at(0)->Fill(mass->at(g), mass->at(f));

                        /*
                        if (fromD0->at(f) == -101 && fromD0->at(g) == -101)
                            nonprompt_phi->Fill(DbarD_phi);
                        if (((vertex_no->at(f) != vertex_no->at(g)) || (vertex_no->at(f) == 0 && vertex_no->at(g) == 0)) && (fromD0->at(f) != -101 || fromD0->at(g) != -101))
                            uncorrelated_dphi->Fill(DbarD_phi);
                        if (vertex_no->at(f) == vertex_no->at(g) && vertex_no->at(f) < 0 && vertex_no->at(g) < 0)
                        {
                            correlated_dphi->Fill(DbarD_phi);
                            if (phy_process->at(f) == 1)
                                Gluon_splitting->Fill(DbarD_phi);
                            if (phy_process->at(f) == 2)
                                Gluon_fusion->Fill(DbarD_phi);
                            if (phy_process->at(f) == 3)
                                flavor_excitation_gluon->Fill(DbarD_phi);
                            if (phy_process->at(f) == 4)
                                flavor_excitation_quark->Fill(DbarD_phi);
                            if (phy_process->at(f) == 5)
                                Quark_annihlation->Fill(DbarD_phi);
                        }
                        */
                    }
                    if (mass->at(f) > 1.6 && mass->at(f) < 2.1 && mass->at(g) > 1.6 && mass->at(g) < 2.1)
                    {
                        F1Mass.at(0)->Fill(mass->at(g));
                        F2Mass.at(0)->Fill(mass->at(f));
                        FullRangeM1M2Mass.at(0)->Fill(mass->at(g), mass->at(f));
                    }

                    /*
                    if ((fromD0->at(f) == 0 || fromD0->at(g) == 0))
                    {
                        background_dphi.at(bin-1)->Fill(DbarD_phi); // background dphi
                    }
                    */

                    //////for dphi correlation below!!!!!!!!!!!!!!!!!!!!//

                    if (fisF2 && gisF1)
                        F1F2.at(0)->Fill(DbarD_phi);
                    if (fisSB2 && gisF1)
                        F1SB2.at(0)->Fill(DbarD_phi);
                    if (fisF2 && gisSB1)
                        SB1F2.at(0)->Fill(DbarD_phi);
                    if (fisSB2 && gisSB1)
                    {
                        SB1SB2.at(0)->Fill(DbarD_phi);
                        // SB1Mass->Fill(mass->at(g));
                        // SB2Mass->Fill(mass->at(f));
                    }
                    if (gisB1 && fisB2)
                    {
                        B1B2.at(0)->Fill(DbarD_phi);
                        // B1Mass->Fill(mass->at(g));
                        // B2Mass->Fill(mass->at(f));
                    }
                    if (gisB1 && fisF2)
                        B1F2.at(0)->Fill(DbarD_phi);
                    if (gisF1 && fisB2)
                        F1B2.at(0)->Fill(DbarD_phi);

                    // if(gisS1 && fisS2) S1S2->Fill(DbarD_phi);
                    if (gisS1 && fisB2)
                        S1B2.at(0)->Fill(DbarD_phi);
                    if (gisB1 && fisS2)
                        B1S2.at(0)->Fill(DbarD_phi);

                } // D0barcand then D0can

            bin=0;
            } // g loop

        } // f loop

    } // for loop

    TCanvas *c = new TCanvas("c", "", 1200, 800);

    /*
    c->Divide(3, 2);
    c->cd(1);
    inclusiveSignal->Draw();
    c->cd(2);
    correlated_dphi->Draw();
    c->cd(3);
    uncorrelated_dphi->Draw();
    c->cd(4);
    background_dphi->Draw();
    c->cd(5);
    nonprompt_phi->Draw();
    c->cd(6);
    auto hs = new THStack("hs", "");
    hs->SetTitle("Dphi by Physics Process");
    Gluon_splitting->SetFillColor(kBlue);
    Gluon_fusion->SetFillColor(kRed);
    flavor_excitation_gluon->SetFillColor(kGreen);
    flavor_excitation_quark->SetFillColor(kCyan);
    Quark_annihlation->SetFillColor(kYellow);
    hs->Add(Quark_annihlation);
    hs->Add(flavor_excitation_quark);
    hs->Add(flavor_excitation_gluon);
    hs->Add(Gluon_fusion);
    hs->Add(Gluon_splitting);
    hs->Draw();
    c->SaveAs("test_histograms_massplots.pdf");

    TCanvas *c1 = new TCanvas("c1", "", 1200, 800);
    c1->Divide(3,2);
    c1->cd(1);
    B1B2->Draw();
    c1->cd(2);
    B1Mass->Draw();
    c1->cd(3);
    B2Mass->Draw();
    c1->cd(4);
    SB1SB2->Draw();
    c1->cd(5);
    SB1Mass->Draw();
    c1->cd(6);
    SB2Mass->Draw();
    c1->SaveAs("sidebands_background_check_after_cuts.pdf");
    */

    TCanvas *c2 = new TCanvas("c2", "", 1200, 1200);
    c2->Divide(2,2);
    c2->cd(1);
    F1F2.at(0)->Draw();
    c2->cd(2);
    F1SB2.at(0)->Draw();
    c2->cd(3);
    SB1F2.at(0)->Draw();
    c2->cd(4);
    SB1SB2.at(0)->Draw();
    //c2->SaveAs("dphi_closure_test_components.pdf");

    TFile *results = new TFile(outfile, "recreate");
    results->cd();
    // SB1Mass->Write();
    // SB2Mass->Write();
    // SW1Mass->Write();
    // SW2Mass->Write();
    // SBSW1Mass->Write();
    // SBSW2Mass->Write();
    // S1Mass->Write();
    // S2Mass->Write();
    F1Mass.at(0)->Write();
    F2Mass.at(0)->Write();
    FullRangeM1M2Mass.at(0)->Write();
    // S1S2inclusiveSignalMass->Write();
    // B1Mass->Write();
    // B2Mass->Write();
    S1B2.at(0)->Write();
    B1S2.at(0)->Write();
    F1F2.at(0)->Write();
    F1SB2.at(0)->Write();
    SB1F2.at(0)->Write();
    SB1SB2.at(0)->Write();
    B1B2.at(0)->Write();
    B1F2.at(0)->Write();
    F1B2.at(0)->Write();
    // S1S2.at(0)->Add(F1F2.at(0));
    // S1S2.at(0)->Add(SB1SB2.at(0));
    // S1S2.at(0)->Add(F1SB2.at(0), -1);
    // S1S2.at(0)->Add(SB1F2.at(0), -1);
    S1S2.at(0)->Add(F1F2.at(0));
    S1S2.at(0)->Add(SB1SB2.at(0));
    S1S2.at(0)->Add(F1SB2.at(0), -1);
    S1S2.at(0)->Add(SB1F2.at(0), -1);
    S1S2.at(0)->Write();
    inclusiveSignal.at(0)->Write();

    TH1F *S1S2e = new TH1F("S1S2e", "S1S2e", nbins, xmin, xmax);
    S1S2e->Sumw2();
    S1S2e->Add(F1F2.at(0));
    S1S2e->Add(F1SB2.at(0), -1);
    S1S2e->Add(SB1F2.at(0), -1);
    S1S2e->Add(SB1SB2.at(0));
    S1S2e->Write();

    TCanvas *c3 = new TCanvas("c3", "", 1200, 1200);
    c3->Divide(2,2);
    c3->cd(1);
    inclusiveSignal.at(0)->SetMinimum(0);
    inclusiveSignal.at(0)->SetMaximum(30000);
    inclusiveSignal.at(0)->SetLineColor(kBlack);
    inclusiveSignal.at(0)->Draw();
    c3->cd(2);
    S1S2e->SetMinimum(0);
    S1S2e->SetTitle("S1S2 via SB Method");
    S1S2e->SetLineColor(kRed);
    S1S2e->SetMaximum(30000);
    inclusiveSignal.at(0)->SetMaximum(30000);
    S1S2e->Draw();
    c3->cd(3);

    TH1F *hist = new TH1F("hist", "S1S2 Template Fit", nbins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
    double values[5] = {21780.5, 23632.3, 22999.0, 26579.6, 23930.4};
    double errors[5] = {465.707, 440.069, 431.497, 536.060, 560.55};

    for (int i = 1; i < 6; i++)
    {
        hist->SetBinContent(i, values[i - 1]);
        hist->SetBinError(i, errors[i - 1]);
    }
    hist->SetMinimum(0);
    hist->SetMaximum(30000);
    hist->SetLineColor(kBlue);
    hist->Draw("E");
    c3->cd(4);
    inclusiveSignal.at(0)->Draw();
    S1S2e->Draw("same");
    hist->Draw("same");



    cout << "---" << endl;
    for (int i=1; i<6; i++){
        cout << "for bin " << i << " fraction of S1S2 / inclusive signal = " << S1S2.at(0)->GetBinContent(i) / inclusiveSignal.at(0)->GetBinContent(i) << endl;
        cout << "---" << endl;
    }
    cout << "overall closure = " << S1S2.at(0)->GetEntries() / inclusiveSignal.at(0)->GetEntries() << endl;
    cout << "---" << endl;
    /*
    correlated_dphi->Write();
    uncorrelated_dphi->Write();
    background_dphi->Write();
    nonprompt_phi->Write();
    Gluon_splitting->Write();
    Gluon_fusion->Write();
    flavor_excitation_gluon->Write();
    flavor_excitation_quark->Write();
    Quark_annihlation->Write();
    */

    return 0;
}
