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

// const char *file_path = "/scratch/bell/awesole/analysis_hepmc_07_18_updatedBkg/ROOT/july29mergedfile.root";
//const char *file_path = "data_files/mergedfile_sept24_2024.root";
const char *file_path = "data_files/mergedfile29Sept2024_cuts.root";
//const char *file_path = "data_files/bin4test.root";
// const char *file_path = "/scratch/bell/awesole/test_57.root";
const float upper_limit = 3.14159, lower_limit = -3.14159;          // bounds to transition delta phi from (-2pi, 2pi) to (-0.5pi, 1.5*pi)
const float nbins_mass = 1.0e+04, xmin_mass = 1.6, xmax_mass = 2.1; // for mass plots
const double nbins = 1000.0, xmin = -2.0, xmax = 5.0;
float ymin = 0.0, ymax = 1.5;                        // for dphi plots
const double avg_mass = 1.865, std_mass = 1.902e-02; // July 30 updates
const float pT_minimum = 0.5;                        //   GeV/c
const float in_val = -0.5 * TMath::Pi(), sc_val = (2 * TMath::Pi() / 5);
const float phi_array[6] = {in_val, in_val + sc_val, in_val + 2 * sc_val, in_val + 3 * sc_val, in_val + 4 * sc_val, in_val + 5 * sc_val};

bool fisB, fisS, fisF, fisSW;
bool gisB, gisS, gisF, gisSW;
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

void MakeSingleParticleMassPlots_5bins()
{
    //TString outfile = TString("TH2F_output_10k.root");
    TString outfile = TString("data_files/TH2F_skippingswap.root");
    TFile *results = new TFile(outfile, "recreate");

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

    std::vector<TH2F *> M1M2Mass;
    std::vector<TH2F *> S1S2Mass;
    std::vector<TH2F *> S1SW2Mass;
    std::vector<TH2F *> SW1S2Mass;
    std::vector<TH2F *> B1B2Mass;
    std::vector<TH2F *> SW1SW2Mass;
    std::vector<TH2F *> SignalSwapMass;
    std::vector<TH2F *> S1B2Mass;
    std::vector<TH2F *> B1S2Mass;
    std::vector<TH2F *> SW1B2Mass;
    std::vector<TH2F *> B1SW2Mass;
    std::vector<TNtuple *> MassTups;

    for (int i = 1; i < 7; i++)
    {
        TString name1 = TString::Format("M1M2Mass_%d", i);
        TString name2 = TString::Format("S1S2Mass_%d", i);
        TString name3 = TString::Format("S1SW2Mass_%d", i);
        TString name4 = TString::Format("SW1S2Mass_%d", i);
        TString name5 = TString::Format("B1B2Mass_%d", i);
        TString name6 = TString::Format("SW1SW2Mass_%d", i);
        TString name7 = TString::Format("SignalSwapMass_%d", i);
        TString name8 = TString::Format("S1B2Mass_%d", i);
        TString name9 = TString::Format("B1S2Mass_%d", i);
        TString name10 = TString::Format("SW1B2Mass_%d", i);
        TString name11 = TString::Format("B1SW2Mass_%d", i);
        TString name12 = TString::Format("MassTuple_%d", i);

        TNtuple *nt = new TNtuple(name12, name12 , "M1:M2");

        TH2F *M1M2 = new TH2F(name1, name1, nbins, 1.6, 2.1, nbins, 1.6, 2.1);
        M1M2->SetXTitle("M1 Mass");
        M1M2->SetYTitle("M2 Mass");
        M1M2->SetOption("SURF1");
        TH2F *S1S2 = new TH2F(name2, name2, nbins, 1.6, 2.1, nbins, 1.6, 2.1);
        S1S2->SetXTitle("S1 Mass");
        S1S2->SetYTitle("S2 Mass");
        S1S2->SetOption("SURF1");
        TH2F *S1SW2 = new TH2F(name3, name3, nbins, 1.6, 2.1, nbins, 1.6, 2.1);
        S1SW2->SetXTitle("S1 Mass");
        S1SW2->SetYTitle("SW2 Mass");
        S1SW2->SetOption("SURF1");
        TH2F *SW1S2 = new TH2F(name4, name4, nbins, 1.6, 2.1, nbins, 1.6, 2.1);
        SW1S2->SetXTitle("SW1 Mass");
        SW1S2->SetYTitle("S2 Mass");
        SW1S2->SetOption("SURF1");
        TH2F *B1B2 = new TH2F(name5, name5, nbins, 1.6, 2.1, nbins, 1.6, 2.1);
        B1B2->SetXTitle("B1 Mass");
        B1B2->SetYTitle("B2 Mass");
        B1B2->SetOption("SURF1");
        TH2F *SW1SW2 = new TH2F(name6, name6, nbins, 1.6, 2.1, nbins, 1.6, 2.1);
        SW1SW2->SetXTitle("SW1 Mass");
        SW1SW2->SetYTitle("SW2 Mass");
        SW1SW2->SetOption("SURF1");
        TH2F *SignalSwap = new TH2F(name7, name7, nbins, 1.6, 2.1, nbins, 1.6, 2.1);
        SignalSwap->SetXTitle("SSW1 Mass");
        SignalSwap->SetYTitle("SSW2 Mass");
        SignalSwap->SetOption("SURF1");
        TH2F *S1B2 = new TH2F(name8, name8, nbins, 1.6, 2.1, nbins, 1.6, 2.1);
        S1B2->SetXTitle("S1 Mass");
        S1B2->SetYTitle("B2 Mass");
        S1B2->SetOption("SURF1");
        TH2F *B1S2 = new TH2F(name9, name9, nbins, 1.6, 2.1, nbins, 1.6, 2.1);
        B1S2->SetXTitle("B1 Mass");
        B1S2->SetYTitle("S2 Mass");
        B1S2->SetOption("SURF1");
        TH2F *SW1B2 = new TH2F(name10, name10, nbins, 1.6, 2.1, nbins, 1.6, 2.1);
        SW1B2->SetXTitle("SW1 Mass");
        SW1B2->SetYTitle("B2 Mass");
        SW1B2->SetOption("SURF1");
        TH2F *B1SW2 = new TH2F(name11, name11, nbins, 1.6, 2.1, nbins, 1.6, 2.1);
        B1SW2->SetXTitle("B1 Mass");
        B1SW2->SetYTitle("SW2 Mass");
        B1SW2->SetOption("SURF1");

        M1M2Mass.push_back(M1M2);
        S1S2Mass.push_back(S1S2);
        S1SW2Mass.push_back(S1SW2);
        SW1S2Mass.push_back(SW1S2);
        B1B2Mass.push_back(B1B2);
        SW1SW2Mass.push_back(SW1SW2);
        SignalSwapMass.push_back(SignalSwap);
        S1B2Mass.push_back(S1B2);
        B1S2Mass.push_back(B1S2);
        SW1B2Mass.push_back(SW1B2);
        B1SW2Mass.push_back(B1SW2);
        MassTups.push_back(nt);
    }

    bool debug = false;
    float Ifile, event_no, D0_index;
    float pT_i, pT_j;
    int bin = 0, tally=0;

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
    // for (int i = 0; i < 20; i++)
    {
        t->GetEntry(i);
        if (phy_process->size() <= 1)
            continue;
         /*
         if(Ifile!=2) continue;
         if(event_no!=14) continue;
         */
        // cout << "Ifile =" << Ifile << " and event_no =" << event_no << endl;
        if (i % 6000 == 0)
            cout << i << " / " << t->GetEntries() << "  " << 100 * i / t->GetEntries() << "%" << endl;

        for (int f = 0; f < phy_process->size() - 1; f++)
        { // for first kpi pair

            fisB = fisS = fisF = fisSW = fisD0candidate = fisDbarcandidate = f_foreground = f_sidebands = false;

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
            if(AssignedPdg_i->at(f) == -TruePdg_j->at(f) && AssignedPdg_j->at(f) == -TruePdg_i->at(f)) continue;//skip all swap
            // if (fromD0->at(f) == 0) continue;

            // determine D0 or Dbar candidate
            if ((AssignedPdg_i->at(f) == 211 && AssignedPdg_j->at(f) == -321) || (AssignedPdg_i->at(f) == -321 && AssignedPdg_j->at(f) == 211))
                fisD0candidate = true;
            if ((AssignedPdg_i->at(f) == -211 && AssignedPdg_j->at(f) == 321) || (AssignedPdg_i->at(f) == 321 && AssignedPdg_j->at(f) == -211))
                fisDbarcandidate = true;
            if (AssignedPdg_i->at(f) == TruePdg_i->at(f) && AssignedPdg_j->at(f) == TruePdg_j->at(f) && fromD0->at(f) != 0)
                fisS = true;
            if (AssignedPdg_i->at(f) == -TruePdg_j->at(f) && AssignedPdg_j->at(f) == -TruePdg_i->at(f) && fromD0->at(f) != 0)
                fisSW = true;
            if (fromD0->at(f) == 0)
                fisB = true;

            if (mass->at(f) > 1.6 && mass->at(f) < 2.1)
            {
                for (int g = f + 1; g < phy_process->size(); g++)
                {
                    gisB = gisS = gisF = gisSW = gisD0candidate = gisDbarcandidate = g_foreground = g_sidebands = false;
                    bin = 0;

                    pT_i = 0.0;
                    pT_j = 0.0;
                    if (prob->at(g) > prob_cutoff && fromD0->at(g) == 0)
                        continue; // background cuts
                    if (phy_process->at(g) > 5)
                        continue; // skip all phy_process that are not 1-5
                    if (iparticle_px->at(f) == iparticle_px->at(g) && iparticle_py->at(f) == iparticle_py->at(g) && iparticle_pz->at(f) == iparticle_pz->at(g) 
                        && jparticle_px->at(f) == jparticle_px->at(g) && jparticle_py->at(f) == jparticle_py->at(g) && jparticle_pz->at(f) == jparticle_pz->at(g)) {
                         tally +=1;
                         continue;
                        }
                        //cout  << "skipping! this is particle && its own swap mass at f = " << mass->at(f) << " and g = " << mass->at(g) << " f and g = " << f << "&" << g << endl;

                    // if(fromD0->at(g)!=0 && AssignedPdg_i->at(g) == TruePdg_i->at(g) && AssignedPdg_j->at(g) == TruePdg_j->at(g)) continue;//skip all signal
                    if(AssignedPdg_i->at(g) == -TruePdg_j->at(g) && AssignedPdg_j->at(g) == -TruePdg_i->at(g)) continue;//skip all swap
                    // if(fromD0->at(g) ==0) continue;

                    pT_i = std::sqrt(iparticle_px->at(g) * iparticle_px->at(g) + iparticle_py->at(g) * iparticle_py->at(g));
                    pT_j = std::sqrt(jparticle_px->at(g) * jparticle_px->at(g) + jparticle_py->at(g) * jparticle_py->at(g));
                    if (pT_i < pT_minimum || pT_j < pT_minimum)
                        continue; // apply pT cuts to each daughter (iparticle & jparticle)

                    float DDbar_phi = phi->at(f) - phi->at(g); // store the dphi if hte kpi pair are from D0candidate then D0bar candidate
                    float DbarD_phi = phi->at(g) - phi->at(f); // store the dphi if hte kpi pair are from D0 bar candidate then D0candidate

                    transition_phi(DDbar_phi); // transition the phi to be within (-0.5pi, 1.5pi)
                    /*
                    if (abs(DDbar_phi) < 1.0e-6)
                        continue;
                        */
                    transition_phi(DbarD_phi);
                    /*
                    if (abs(DbarD_phi) < 1.0e-6)
                        continue;
                        */

                    // determine D0 and Dbar candidate:
                    if ((AssignedPdg_i->at(g) == 211 && AssignedPdg_j->at(g) == -321) || (AssignedPdg_i->at(g) == -321 && AssignedPdg_j->at(g) == 211))
                        gisD0candidate = true;
                    if ((AssignedPdg_i->at(g) == -211 && AssignedPdg_j->at(g) == 321) || (AssignedPdg_i->at(g) == 321 && AssignedPdg_j->at(g) == -211))
                        gisDbarcandidate = true;
                    if (AssignedPdg_i->at(g) == TruePdg_i->at(g) && AssignedPdg_j->at(g) == TruePdg_j->at(g) && fromD0->at(g) != 0)
                        gisS = true;
                    if (AssignedPdg_i->at(g) == -TruePdg_j->at(g) && AssignedPdg_j->at(g) == -TruePdg_i->at(g) && fromD0->at(g) != 0)
                        gisSW = true;
                    if (fromD0->at(g) == 0)
                        gisB = true;

                    if ((fisD0candidate && gisD0candidate) || (fisDbarcandidate && gisDbarcandidate))
                        continue;
                    if (mass->at(g) > 1.6 && mass->at(g) < 2.1)
                    {
                        if (fisD0candidate && gisDbarcandidate)
                        {
                            if (bin != 0)
                                cout << "Error!!!!!!!  phi = " << DDbar_phi << " and bin =" << bin << endl;
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
                                //if (bin != 4) continue;
                            M1M2Mass.at(bin - 1)->Fill(mass->at(f), mass->at(g)); // signal+swap+bkg
                            MassTups.at(bin - 1)->Fill(mass->at(f), mass->at(g)); 
                            MassTups.at(5)->Fill(mass->at(f), mass->at(g)); 

                            if ((fisS || fisSW) && (gisS || gisSW))
                                SignalSwapMass.at(bin - 1)->Fill(mass->at(f), mass->at(g));
                            if (fisS && gisS)
                                S1S2Mass.at(bin - 1)->Fill(mass->at(f), mass->at(g)); // only signal
                            if (fisSW && gisS)
                                SW1S2Mass.at(bin - 1)->Fill(mass->at(f), mass->at(g)); // only signal
                            if (fisS && gisSW)
                                S1SW2Mass.at(bin - 1)->Fill(mass->at(f), mass->at(g)); // only signal
                            if (fisSW && gisSW)
                                SW1SW2Mass.at(bin - 1)->Fill(mass->at(f), mass->at(g)); // only signal
                            if (fisB && gisB)
                                B1B2Mass.at(bin - 1)->Fill(mass->at(f), mass->at(g));
                            if (fisS && gisB)
                                S1B2Mass.at(bin - 1)->Fill(mass->at(f), mass->at(g)); // s1b2
                            if (fisSW && gisB)
                                SW1B2Mass.at(bin - 1)->Fill(mass->at(f), mass->at(g)); // sw1b2
                            if (fisB && gisS)
                                B1S2Mass.at(bin - 1)->Fill(mass->at(f), mass->at(g)); // s1b2
                            if (fisB && gisSW)
                                B1SW2Mass.at(bin - 1)->Fill(mass->at(f), mass->at(g)); // sw1b2
                            /*
                            else{
                                cout << "Error dphi not filled!! bin = " << bin << " phi=" << DDbar_phi  << " fromD0 at f = " << fromD0->at(f) << " and g= " << fromD0->at(g) << endl;
                                cout << "gisS = " << gisS << endl;
                            }
                            */ 
                        }
                        else if (fisDbarcandidate && gisD0candidate)
                        {
                            if (bin != 0)
                                cout << "Error!!!!!!!  phi = " << DbarD_phi << " and bin =" << bin << endl;
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
                               // if (bin != 4) continue;
                            M1M2Mass.at(bin - 1)->Fill(mass->at(g), mass->at(f)); // signal+swap+bkg
                            MassTups.at(bin - 1)->Fill(mass->at(g), mass->at(f)); 
                            MassTups.at(5)->Fill(mass->at(g), mass->at(f)); 

                            if ((fisS || fisSW) && (gisS || gisSW))
                                SignalSwapMass.at(bin - 1)->Fill(mass->at(g), mass->at(f));
                            if (gisS && fisS)
                                S1S2Mass.at(bin - 1)->Fill(mass->at(g), mass->at(f)); // only signal
                            if (gisSW && fisSW)
                                SW1SW2Mass.at(bin - 1)->Fill(mass->at(g), mass->at(f)); // only signal
                            if (gisS && fisSW)
                                S1SW2Mass.at(bin - 1)->Fill(mass->at(g), mass->at(f)); // only signal
                            if (gisSW && fisS)
                                SW1S2Mass.at(bin - 1)->Fill(mass->at(g), mass->at(f)); // only signal
                            if (gisB && fisB)
                                B1B2Mass.at(bin - 1)->Fill(mass->at(g), mass->at(f));
                            if (gisS && fisB)
                                S1B2Mass.at(bin - 1)->Fill(mass->at(g), mass->at(f)); // s1b2
                            if (gisSW && fisB)
                                SW1B2Mass.at(bin - 1)->Fill(mass->at(g), mass->at(f)); // sw1b2
                            if (gisB && fisS)
                                B1S2Mass.at(bin - 1)->Fill(mass->at(g), mass->at(f)); // s1b2
                            if (gisB && fisSW)
                                B1SW2Mass.at(bin - 1)->Fill(mass->at(g), mass->at(f)); // sw1b2
                            /*
                            else
                                cout << "Error dphi not filled!! bin = " << bin << " phi=" << DbarD_phi << endl;
                                */
                        }
                    }
                }
            }
        }
    }
    results->cd();

    for (int i = 0; i < 5; i++)
    {
        M1M2Mass.at(5)->Add(M1M2Mass.at(i));
        S1S2Mass.at(5)->Add(S1S2Mass.at(i));
        SW1SW2Mass.at(5)->Add(SW1SW2Mass.at(i));
        B1B2Mass.at(5)->Add(B1B2Mass.at(i));
        SignalSwapMass.at(5)->Add(SignalSwapMass.at(i));
        SW1S2Mass.at(5)->Add(SW1S2Mass.at(i));
        S1SW2Mass.at(5)->Add(S1SW2Mass.at(i));
        S1B2Mass.at(5)->Add(S1B2Mass.at(i));
        B1S2Mass.at(5)->Add(B1S2Mass.at(i));
        SW1B2Mass.at(5)->Add(SW1B2Mass.at(i));
        B1SW2Mass.at(5)->Add(B1SW2Mass.at(i));
    }

    for (int i = 0; i < 6; i++)
    {
        M1M2Mass.at(i)->Write();
        S1S2Mass.at(i)->Write();
        SW1SW2Mass.at(i)->Write();
        B1B2Mass.at(i)->Write();
        SignalSwapMass.at(i)->Write();
        SW1S2Mass.at(i)->Write();
        S1SW2Mass.at(i)->Write();
        S1B2Mass.at(i)->Write();
        B1S2Mass.at(i)->Write();
        SW1B2Mass.at(i)->Write();
        B1SW2Mass.at(i)->Write();
        MassTups.at(i)->Write();
    }
    for (int i = 0; i < 6; i++)
    {
        M1M2Mass.at(i)->Add(S1S2Mass.at(i), -1);
        M1M2Mass.at(i)->Add(S1B2Mass.at(i), -1);
        M1M2Mass.at(i)->Add(B1S2Mass.at(i), -1);
        M1M2Mass.at(i)->Add(S1SW2Mass.at(i), -1);
        M1M2Mass.at(i)->Add(SW1S2Mass.at(i), -1);
        M1M2Mass.at(i)->Add(SW1SW2Mass.at(i), -1);
        M1M2Mass.at(i)->Add(B1B2Mass.at(i), -1);
        M1M2Mass.at(i)->Add(B1SW2Mass.at(i), -1);
        M1M2Mass.at(i)->Add(SW1B2Mass.at(i), -1);
        if (M1M2Mass.at(i)->GetEntries() != 0)
            cout << "Error m1m2mass not correct!   bin=" << i + 1 << "entries =" << M1M2Mass.at(i)->GetEntries() << endl;
    }
}
