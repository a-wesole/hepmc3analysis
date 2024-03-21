#include <iostream>
#include <Pythia8/Pythia.h>
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>
#include <HepMC3/WriterAscii.h>
#include <HepMC3/ReaderAscii.h>
#include "Math/Vector4D.h"
#include "TLorentzVector.h"
#include <vector>
#include <TTree.h>
#include <TBranch.h>
#include <cmath>
#include <TMath.h>
#include <TRandom3.h>


#include "TFile.h"
#include "TNtuple.h"


using namespace std;

void tracecdaughters(const std::shared_ptr<const HepMC3::GenParticle>& particle1, bool& D0, float& D0px, float& D0py, float& D0pz, float& D0e, float& Kpx, float& Kpy, float& Kpz, float& Ke, float& Pipx, float& Pipy, float& Pipz, float& Pie);
void tracecbardaughters(const std::shared_ptr<const HepMC3::GenParticle>& particle2, bool& D0bar, float& D0barpx, float& D0barpy, float& D0barpz);
void tracecmothers(const std::shared_ptr<const HepMC3::GenParticle>& InParticle2, int c_index, bool& flavor_exc);
void tracecbarmothers(const std::shared_ptr<const HepMC3::GenParticle>& InParticle2, int cbar_index, bool& flavor_exc);



//int main(int istart, int iend) {
int main(int argc, char *argv[]) {

    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " istart iend" << std::endl;
        return 1;  // Indicate an error
    }

    int istart = std::stoi(argv[1]);
    int iend = std::stoi(argv[2]);
    //TString outfile = TString("garbage.root");
   TString outfile = TString::Format("/scratch/bell/awesole/analysis_hepmc_0318/ROOT/hepmc_tree_%d_%d.root",istart, iend);
        TNtuple *T = new TNtuple("T", "T","event_no:vtx_no:D0_count:Dbar0_count:phys_process:pT_c:pT_cbar:delta_phi:delta_eta:color:delta_p:index:pT_D0:pT_D0bar:D0_phi:D0bar_phi:D0_delta_phi:D0_eta:D0bar_eta:c_eta:cbar_eta:Kmass:PiMass:D0mass:smearD0mass:deltapxK:deltapyK:deltapzK:deltapxPi:deltapyPi:deltapzPi");
    TFile *results = new TFile(outfile, "recreate");
    //ifstream file_stream("pythia_output_02_19.txt");
    ifstream file_stream("D0_updates_0312.txt");

    string filename;

    int ifile = 0;

    while (true) {

        file_stream >> filename;
        if (file_stream.eof()) {
            break;
        }
        if (ifile < istart) {
            ifile++;
            continue;
        }
        if (ifile >= iend) break;

        cout << "ifile=" << ifile << endl;
        cout << "outfile=" << outfile << endl;
        //TFile *fin = TFile::Open(filename.c_str());
        //TTree *tree = (TTree *) fin->Get("T");

        //std::string hepmcFilePath;
        //std::getline(file_stream, hepmcFilePath);
        cout << "file path!= " << filename.c_str() << endl;
        HepMC3::ReaderAscii reader(filename.c_str());

        /*
           if (fin->IsZombie()) {
           ifile++;
           continue;
           }
           */

        int debug = 1;

        int c_index = -1;
        int cbar_index = -1;
        bool D0, D0bar;
        float c_phi, cbar_phi, c_eta, cbar_eta, c_theta, cbar_theta, cpx, cpy, cpz, cbarpx, cbarpy, cbarpz, gluon_momentum;
        float D0_phi, D0bar_phi, D0_eta, D0bar_eta, D0_theta, D0bar_theta, D0px, D0py, D0pz, D0barpx, D0barpy, D0barpz, mass, D0e;
        float K_phi, K_eta, K_theta, Kpx, Kpy, Kpz, Kmass, Ke;
        float Pi_phi, Pi_eta, Pi_theta, Pipx, Pipy, Pipz, Pimass, Pie;
        float upper_limit = 3.14159, lower_limit = -3.14159;


        //std::ofstream textfile("HepD0info.txt");
        //TFile outFile("outfile.hepmc", "RECREATE");
        //TFile results("analysis_hepmc.root", "RECREATE");
        //const std::string hepmcFilePath = "/scratch/bell/awesole/D0out_02_19_24_1/job0_1.hepmc"; //define file path
        // const std::string hepmcFilePath = "/home/awesole/pythia8/new_hepmc_200_event.hepmc"; //define file path
        //const std::string hepmcFilePath = "/home/awesole/pythia8/small.hepmc"; //define file path
        /////HepMC3::ReaderAscii reader(hepmcFilePath); //open the HepMC file

        float array[31];
        TLorentzVector particleVector;
        TLorentzVector KVector;
        TLorentzVector PiVector;
        TLorentzVector D0Vector;
        TLorentzVector smearKVector;
        TLorentzVector smearPiVector;
        TLorentzVector smearD0Vector;


        int phys_process;

        while (!reader.failed()) { //loop over all events
            HepMC3::GenEvent event; //load the next event
            reader.read_event(event); //read the event
            int c_count = 0, cbar_count = 0;
            int vc_count = 0, vcbar_count = 0;
            int color;
            array[11] = iend;

            //if (event.event_number() != 9832) continue;

            int d0c = 0;
            int d0bc = 0;

            for (const auto &particle1: event.particles()) { //preliminary for loop that only selects events with the same number of c and c-bar
                if (particle1->pdg_id() == 4) c_count++;
                if (particle1->pdg_id() == -4) cbar_count++;
                if (particle1->pdg_id() == 421) d0c++;
                if (particle1->pdg_id() == -421) d0bc++;
            }

            if (c_count > 0 && cbar_count == c_count) {//begin analysis only for events with same number of c and cbar
                //count/cbarcount=" << c_count << "/" << cbar_count << endl;
                D0 = false;
                D0bar = false;

                for (const auto &vertex: event.vertices()) {//for all vertices
                    bool flavor_exc = false;


                    const std::vector<HepMC3::GenParticlePtr> &particlesOut = vertex->particles_out();//create a vector to contain all the incoming particles from a vertex
                    const std::vector<HepMC3::GenParticlePtr> &particlesIn = vertex->particles_in();//create a vector to contain all the incoming particles from a vertex
                    int count4 = 0;
                    int count_4 = 0;


                    //Loop over all outgoing particles
                    for (const auto &particle: particlesOut) {

                        if (particle->pdg_id() == 4) {
                            count4++;  //keep a tally of all c quarks in the vertex
                        }
                        if (particle->pdg_id() == -4) {
                            count_4++; //keep a tally of all cbar quarks
                        }

                    }//loop to count number of c and cbar coming out of each vertex
                    if (count4 > 0 && count4 == count_4) {//only for vertices containing at least 1 c,cbar pair
                                    if (debug > 1) {
                                        cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
                                        cout << "Event number=" << event.event_number() << endl;
                                        cout << "Number of particles=" << event.particles().size() << endl;
                                        cout << "*****New Vertex****" << endl;
                                        cout << "Vertex no. =" << vertex->id() << endl;
                                    }//debug
                        //cout << "LOOK HERE! VERTEX SIZE= " << particlesOut.size() << endl;
                        //D0=false;
                        //D0bar=false;
                        D0 = false;
                        D0bar = false;

                        D0px = 0;
                        D0py = 0;
                        D0pz = 0;
                        D0barpx = 0;
                        D0barpy = 0;
                        D0barpz = 0;

                        for (const auto &particle3: particlesOut) {//loop over all outoging particles to determine physics process, main loop

                            //if (particle3->status() < 60 && particle3->status() > 50) {//only for particles from final state showers, gluon splitting  only
                                if (particle3->pdg_id() == 4) {//trace cquark daughters
                                    // cout << "here we go to start tracing charm" << endl;
                                    tracecdaughters(particle3, D0, D0px, D0py, D0pz, D0e, Kpx, Kpy, Kpz, Ke, Pipx, Pipy,
                                                    Pipz, Pie);
                                }

                                //if (particle3->pdg_id()==-4 && particle3->status()>49){
                                if (particle3->pdg_id() == -4) {//tracecbar daughters
                                    //cout << "here we go to start tracing cbar" << endl;
                                    tracecbardaughters(particle3, D0bar, D0barpx, D0barpy, D0barpz);
                                }
                            //}//only for particles from final state showers, gluon splitting  only
                            if (D0 and D0bar) {//if the vertex leads to a D0 and D0bar
                                if (debug > 0) {
                                    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
                                    cout << "Event number=" << event.event_number() << endl;
                                    cout << "Number of particles=" << event.particles().size() << endl;
                                    cout << "*****New Vertex****" << endl;
                                    cout << "Vertex no. =" << vertex->id() << endl;
                                    cout << "D0 momentum/energy= " << D0px << "/" << D0py << "/" << D0pz << "/" << D0e << endl;
                                }

                                D0 = false;
                                D0bar = false;
                                KVector.SetPxPyPzE(Kpx,Kpy,Kpz,Ke);
                                array[21]=KVector.M();
                                PiVector.SetPxPyPzE(Pipx, Pipy, Pipz, Pie);
                                array[22]=PiVector.M();
                                D0Vector=KVector + PiVector;
                                array[23]=D0Vector.M();


                                array[0] = event.event_number();
                                array[1] = vertex->id();
                                array[2] = d0c;
                                array[3] = d0bc;
                                float D0_pT = std::sqrt((D0px * D0px) + (D0py * D0py));//D0pT
                                float D0bar_pT = std::sqrt((D0barpx * D0barpx) + (D0barpy * D0barpy));//D0pT
                                array[12] = D0_pT;
                                array[13] = D0bar_pT;//D0barpT
                                float D0_phi = std::atan2(D0py, D0px);
                                float D0bar_phi = std::atan2(D0barpy, D0barpx);
                                array[14] = D0_phi;//D0_phi
                                array[15] = D0bar_phi;//D0_phi
                                D0_theta = std::atan2(D0_pT, D0pz);
                                D0_eta = -1 * std::log(std::tan(D0_theta / 2.0));
                                array[17] = D0_eta;
                                D0bar_theta = std::atan2(D0bar_pT, D0barpz);
                                D0bar_eta = -1 * std::log(std::tan(D0bar_theta / 2.0));
                                array[18] = D0bar_eta;

                                float D0del_phi = D0_phi - D0bar_phi;
                                color = 0;
                                if (D0del_phi > upper_limit) {
                                    D0del_phi = D0del_phi - 2 * TMath::Pi();
                                    color = 1;
                                }

                                if (D0del_phi < lower_limit) {
                                    D0del_phi = 2 * TMath::Pi() + D0del_phi;
                                    color = -1;
                                }
                                if (D0del_phi < -0.5 * TMath::Pi() && D0del_phi > -1 * TMath::Pi()) {
                                    D0del_phi = 2 * TMath::Pi() + D0del_phi;
                                    color = 2;
                                }
                                array[9] = color;
                                array[16] = D0del_phi;

                                if (debug > 1) {
                                    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
                                    cout << "Event number=" << event.event_number() << endl;
                                    cout << "Number of particles=" << event.particles().size() << endl;
                                    cout << "*****New Vertex****" << endl;
                                    cout << "Vertex no. =" << vertex->id() << endl;
                                }//debug

                                //if (particle3->status() < 60 && particle3->status() > 50) {//only for particles from final state showers, gluon splitting  only

                                //cout << "Particles in include:" << endl;
                                cpx = 0;
                                cpy = 0;
                                cpz = 0;
                                cbarpx = 0;
                                cbarpy = 0;
                                cbarpz = 0;
                                for (const auto &ccbar1: particlesOut) {
                                    if (ccbar1->pdg_id() == 4) {
                                        float cpT = std::sqrt((ccbar1->momentum().x() * (ccbar1->momentum().x())) +
                                                              (ccbar1->momentum().y() *
                                                               (ccbar1->momentum().y())));//pt_c
                                        array[5] = cpT;//pt_c
                                        c_phi = std::atan2(ccbar1->momentum().y(), ccbar1->momentum().x());
                                        c_theta = std::atan2(cpT,ccbar1->momentum().z());
                                        c_eta = -1 * std::log(std::tan(c_theta / 2.0));
                                        array[19] = c_eta;
                                        cpx = ccbar1->momentum().x();
                                        cpy = ccbar1->momentum().y();
                                        cpz = ccbar1->momentum().z();
                                    }
                                    if (ccbar1->pdg_id() == -4) {
                                        float cbarpT = std::sqrt((ccbar1->momentum().x() * (ccbar1->momentum().x())) +
                                                                 (ccbar1->momentum().y() *
                                                                  (ccbar1->momentum().y())));//pTcbar
                                        array[6] = cbarpT;//pTcbar
                                        cbar_phi = std::atan2(ccbar1->momentum().y(), ccbar1->momentum().x());
                                        cbar_theta = std::atan2(cbarpT,ccbar1->momentum().z());
                                        cbar_eta = -1 * std::log(std::tan(cbar_theta / 2.0));
                                        array[20] = cbar_eta;
                                        cbarpx = ccbar1->momentum().x();
                                        cbarpy = ccbar1->momentum().y();
                                        cbarpz = ccbar1->momentum().z();
                                    }

                                    //cout << "pdg/inex" << ccbar->pdg_id() << "/" << ccbar->id() << endl;
                                }
                                float del_phi = c_phi - cbar_phi;
                                color = 0;
                                if (del_phi > upper_limit) {
                                    del_phi = del_phi - 2 * TMath::Pi();
                                    //color =1;
                                }

                                if (del_phi < lower_limit) {
                                    del_phi = 2 * TMath::Pi() + del_phi;
                                    //color =-1;
                                }
                                if (del_phi < -0.5 * TMath::Pi() && del_phi > -1 * TMath::Pi()) {
                                    del_phi = 2 * TMath::Pi() + del_phi;
                                    //color=2;
                                }
                                //if (del_phi  > -1.5*TMath::Pi() && del_phi < -0.5*TMath::Pi()) del_phi=abs(del_phi);
                                array[7] = del_phi;
                                array[8] = c_eta - cbar_eta;
                                gluon_momentum = 0;

                                if (particle3->status() < 60 && particle3->status() > 50) {//only for particles from final state showers, gluon splitting  only
                                    cout << "BOTH D0 AND D0BAR -- PROCEED WITH ANALYSIS OF FINAL STATE SHOWERS" << endl;
                                    phys_process=-1;

                                    for (const auto &particle4: particlesIn) {//if both D0 and D0bar, begin analysis of incoming particles of vertex

                                        if (particlesIn.size() == 1 && particle4->pdg_id() == 21) {
                                            phys_process = 1;//gluon splitting
                                            gluon_momentum = sqrt(
                                                    particle4->momentum().x() * particle4->momentum().x() +
                                                    particle4->momentum().y() * particle4->momentum().y() +
                                                    particle4->momentum().z() * particle4->momentum().z());
                                        }
                                        if (phys_process == 1) {
                                            cout << "-*-*-*-*-Physics process is:";
                                            cout << "gluon splitting-*-*-*-*-" << endl;

                                            float cp = sqrt(
                                                    (cpx + cbarpx) * (cpx + cbarpx) + (cpy + cbarpy) * (cpy + cbarpy) +
                                                    (cpz + cbarpz) * (cpz + cbarpz));

                                            array[10] = gluon_momentum - cp;

                                        } else {
                                            cout << "ERROR: Final State showers not belonging to a category! "
                                                 << endl;
                                            phys_process = -1;
                                        }
                                        array[4] = phys_process;

                                    }//for all incoming particles
                                }//final state

                                else if (particle3->status() < 50) {//only for particles from initial state showers, gluon fusion and flavor excitation
                                phys_process = -4;
                                    cout << "BOTH D0 AND D0BAR -- PROCEED WITH ANALYSIS OF INITIAL STATE SHOWERS" << endl;


                                    for (const auto &particle4: particlesIn) {//if both d0 and d0bar, begin analysis of incoming particles of vertex

                                        //cout << "particle pdg=" << particle4->pdg_id() << ", indexno." << particle4->id() << endl;
                                        //cout << "particle status=" << particle3->status() << endl;

                                        if (particlesIn.size() == 2 && particle4->pdg_id() == 21) {
                                            phys_process = 2;//gluon fusion
                                            array[10] = 0;
                                        }
                                        if (particlesIn.size() == 2 && particle4->pdg_id() < 7) {
                                            phys_process = 4; //quark annihlation
                                            array[10] = 0;
                                        }

                                        if (particlesIn.size() == 1 && particle4->pdg_id() == 21) {//begin analysis to determine flavor excitation
                                            cout << "&&&&&& beginning analysis for flavor excitation &&&&&&" << endl;
                                            int gluon_in = 0, gluon_out = 0, c_in = 0, c_out = 0, cbar_in = 0, cbar_out = 0;
                                            for (const auto &vertex: event.vertices()) {//loop over all vertices to only select the ones with g+c comming in and g+c coming out
                                                for (const auto &inparticle: vertex->particles_in()) {//counts incoming gluons and c/cbar
                                                    if (inparticle->pdg_id() == 21) gluon_in++;
                                                    if (inparticle->pdg_id() == 4) c_in++;
                                                    if (inparticle->pdg_id() == -4) cbar_in++;
                                                }//counts incoming gluons and c/cbar
                                                for (const auto &outparticle: vertex->particles_out()) {//counts outgoing gluons and c/cbar
                                                    if (outparticle->pdg_id() == 21) gluon_out++;
                                                    if (outparticle->pdg_id() == 4) c_out++;
                                                    if (outparticle->pdg_id() == -4) cbar_out++;
                                                }//counts outgoing gluons and c/cbar


                                                for (const auto &inparticle2: vertex->particles_in()) {
                                                    if (inparticle2->pdg_id() == 4 && gluon_in == 1 and c_in == 1 and c_out == 1) {
                                                        int c_index = inparticle2->id();
                                                        //cout << "gin/gout/cin/cout/cbarin/cbarout=" << gluon_in << "/" << gluon_out << "/" << c_in << "/" << c_out << "/" << cbar_in << "/" << cbar_out << endl;
                                                        tracecmothers(inparticle2, c_index, flavor_exc);
                                                    }
                                                    if (inparticle2->pdg_id() == -4 && gluon_in == 1 and cbar_in == 1 and cbar_out == 1) {
                                                        int c_index = inparticle2->id();
                                                        //cout << "gin/gout/cin/cout/cbarin/cbarout=" << gluon_in << "/" << gluon_out << "/" << c_in << "/" << c_out << "/" << cbar_in << "/" << cbar_out << endl;
                                                        tracecbarmothers(inparticle2, cbar_index, flavor_exc);
                                                    }
                                                }
                                                if (flavor_exc) {
                                                    phys_process = 3; //flavor excitation
                                                    array[10] = 0;
                                                }
                                                gluon_in = 0;
                                                gluon_out = 0;
                                                c_in = 0;
                                                c_out = 0;
                                                cbar_in = 0;
                                                cbar_out = 0;

                                            }//for all vertices
                                        }//flavor excitation
                                    }//incoming particles
                                    if (phys_process == 2) {
                                        cout << "-*-*-*-*-physics process is:";
                                        cout << "gluon fusion-*-*-*-*-" << endl;
                                        array[10] = 0;
                                    }
                                    if (phys_process == 3) {
                                        cout << "-*-*-*-*-physics process is:";
                                        cout << "flavor excitation-*-*-*-*-" << endl;
                                        array[10] = 0;
                                    }
                                    if (phys_process == 4) {
                                        cout << "-*-*-*-*-physics process is:";
                                        cout << "quark annihlation -*-*-*-*-" << endl;
                                        array[10] = 0;
                                    }
                                    if (phys_process == -4) {
                                        cout << "error: initial state showers not belonging to a category! " << endl;
                                        array[10] = 0;
                                    } else {//d0 only nothing to be done
                                        //cout <<"only d0, continue" << endl;
                                    }
                                    array[4] = phys_process;
                                }//inital state
                                else continue;
                                //post-analysis
                                //create bit about smearing px py and pz
                                float mean = 0, width = 0.01;
                                gRandom->SetSeed(istart + event.event_number());
                                double smearKpx = Kpx*(1+gRandom->Gaus(mean,width));
                                double smearKpy = Kpy*(1+gRandom->Gaus(mean,width));
                                double smearKpz = Kpz*(1+gRandom->Gaus(mean,width));
                                double smearPipx = Pipx*(1+gRandom->Gaus(mean,width));
                                double smearPipy = Pipy*(1+gRandom->Gaus(mean,width));
                                double smearPipz = Pipz*(1+gRandom->Gaus(mean,width));
                                double Kmass = 0.493677;
                                double Pimass = 0.139570;
                                double  Kenergy = std::sqrt(smearKpx*smearKpx + smearKpy*smearKpy + smearKpz*smearKpz + Kmass*Kmass);
                                double  Pienergy = std::sqrt(smearPipx*smearPipx + smearPipy*smearPipy + smearPipz*smearPipz + Pimass*Pimass);
                                smearKVector.SetPxPyPzE(smearKpx, smearKpy, smearKpz, Kenergy);
                                smearPiVector.SetPxPyPzE(smearPipx, smearPipy, smearPipz,Pienergy);
                                smearD0Vector = smearKVector + smearPiVector;
                                array[24]=smearD0Vector.M();
                                array[25]=(Kpx-smearKpx)/Kpx;//deltapxK
                                array[26]=(Kpy-smearKpy)/Kpy;//deltapxK
                                array[27]=(Kpz-smearKpz)/Kpz;//deltapxK
                                array[28]=(Pipx-smearPipx)/Pipx;
                                array[29]=(Pipy-smearPipy)/Pipy;
                                array[30]=(Pipz-smearPipz)/Pipz;
                                T->Fill(array);
                                cout << "K=" << Kpx << "," << Kpy << "," << Kpz << endl;
                                cout << "smearK=" << smearKpx << "," << smearKpy << "," << smearKpz << endl;
                                cout << "Pi=" << Pipx << "," << Pipy << "," << Pipz << endl;
                                cout << "smearPi=" << smearPipx << "," << smearPipy << "," << smearPipz << endl;
                                cout << "Kenergy/Pienergy =" << Kenergy << "/" << Pienergy << endl;
                            }//if D0 and D0bar
                        }//loop over all outoging particles to determine physics process, main loop
                    }//only for events containing at least 1c/cbar
                }//for all vertices
            }//only for same number of c/cbar
        }//all events loop
            results->cd();
            T->Write();
        //textfile.close();
        reader.close();
        ifile++;
    }//loop ove all files
            cout << "~~~~~~~done~~~~~~~" << endl;
    results->Close();
    return 0;
}//main

            void tracecdaughters(const std::shared_ptr<const HepMC3::GenParticle>& particle6, bool& d0, float& d0px, float& d0py, float& d0pz, float& d0e, float& kpx, float& kpy, float& kpz, float& ke, float& pipx, float& pipy, float& pipz, float& pie){
                if (particle6->pdg_id()==421) {

                    const auto &daughters = particle6->children();
                    if (daughters.size() == 2 && (daughters[0]->pdg_id() == -321 && daughters[1]->pdg_id() == 211) || (daughters[0]->pdg_id() == 211 && daughters[1]->pdg_id() == -321)) {
                        d0 = true;
                        d0px = particle6->momentum().x();
                        d0py = particle6->momentum().y();
                        d0pz = particle6->momentum().z();
                        d0e = particle6->momentum().e();
                        //cout << "great news! this d0 decays into exactly 2 particles!" << endl;
                        //cout << daughters[0]->pdg_id() << "= daughter 1 pdg_id" << endl;
                        //cout << daughters[1]->pdg_id() << "= daughter 2 pdg_id" << endl;
                        //cout << "***************" << endl;
                            for (const auto& daughterparticle:daughters){
                                if (daughterparticle->pdg_id() == -321){
                                    kpx = daughterparticle->momentum().x();
                                    kpy = daughterparticle->momentum().y();
                                    kpz = daughterparticle->momentum().z();
                                    ke = daughterparticle->momentum().e();
                                }
                                if (daughterparticle->pdg_id() == 211){
                                    pipx = daughterparticle->momentum().x();
                                    pipy = daughterparticle->momentum().y();
                                    pipz= daughterparticle->momentum().z();
                                    pie = daughterparticle->momentum().e();
                                }
                            }
                        return;
                    }
                }
                if (d0) return; //stop process if d0bar is already true


                // cout << "locating the daughters for particle with pdg=" << particle6->pdg_id() << "and index no. " << particle6->id() << "status =" << particle6->status() << endl;
                //cout << "locating the daughters for particle with pdg=" << particle6->pdg_id() << " and index no." << particle6->id() << "         ";
                const std::vector<std::shared_ptr<const HepMC3::GenParticle>> &daughters = particle6->children();
                for (const auto &daughter: daughters) {

                    //cout << "daughter index=" << daughter->id() << "daughter pdg_id =" << daughter->pdg_id() << endl;

                    if (daughter->pdg_id() ==421) {
                        const auto& daughters = daughter->children();
                        if(daughters.size()==2 && (daughters[0]->pdg_id() == -321 && daughters[1]->pdg_id() == 211) || (daughters[0]->pdg_id() == 211 && daughters[1]->pdg_id() == -321)) {
                            d0 = true;
                            d0px = daughter->momentum().x();
                            d0py = daughter->momentum().y();
                            d0pz = daughter->momentum().z();
                            d0e = daughter->momentum().e();
                            //cout << "great news! this d0 decays into exactly 2 particles!" << endl;
                        //cout << daughters[0]->pdg_id() << "= daughter 1 pdg_id" << endl;
                        //cout << daughters[1]->pdg_id() << "= daughter 2 pdg_id" << endl;
                        //cout << "***************" << endl;
                        //cout << "located a d0" << "index no. = " << daughter->id() << endl;
                        //cout << "located a d0, " << "momentum =(" << daughter->momentum().x() << "," << daughter->momentum().y() << "," << daughter->momentum().z() << ")" << endl;

                            for (const auto& daughterparticle:daughters){
                                if (daughterparticle->pdg_id() == -321){
                                    kpx = daughterparticle->momentum().x();
                                    kpy = daughterparticle->momentum().y();
                                    kpz = daughterparticle->momentum().z();
                                    ke = daughterparticle->momentum().e();
                                }
                                if (daughterparticle->pdg_id() == 211){
                                    pipx = daughterparticle->momentum().x();
                                    pipy = daughterparticle->momentum().y();
                                    pipz= daughterparticle->momentum().z();
                                    pie = daughterparticle->momentum().e();
                                }
                            }
                        }
                        return;
                    }//if daughter is d0

                    //else tracecdaughters(daughter, d0);
                    else if (daughter->pdg_id()==4 || daughter->status()==74|| (daughter->status()<100 && daughter->status()>80 && daughter->pdg_id()!=21) || (daughter->status()==2 && daughter->pdg_id()!=21)) tracecdaughters(daughter, d0, d0px, d0py, d0pz,d0e,kpx, kpy, kpz, ke, pipx, pipy, pipz, pie );
                }//loop over daughters
            }

            void tracecbardaughters(const std::shared_ptr<const HepMC3::GenParticle>& particle6, bool& d0bar, float& d0barpx, float& d0barpy, float& d0barpz) {

                if (particle6->pdg_id()==-421) {
                    const auto& daughters = particle6->children();
                    if(daughters.size()==2 && (daughters[0]->pdg_id() == 321 && daughters[1]->pdg_id() == -211) || (daughters[0]->pdg_id() == -211 && daughters[1]->pdg_id() == 321)) {
                        d0bar = true;
                        d0barpx = particle6->momentum().x();
                        d0barpy = particle6->momentum().y();
                        d0barpz = particle6->momentum().z();
                        //cout << "great news! this dbar0 decays into exactly 2 particles!" << endl;
                        //cout << daughters[0]->pdg_id() << "= daughter 1 pdg_id" << endl;
                        //cout << daughters[1]->pdg_id() << "= daughter 2 pdg_id" << endl;
                        //cout << "***************" << endl;
                        return;
                    }
                }
                if (d0bar) return; //stop process if d0bar

                // cout << "locating the daughters for particle with pdg=" << particle6->pdg_id() << "and index no. " << particle6->id() << "status =" << particle6->status() << endl;
                const std::vector<std::shared_ptr<const HepMC3::GenParticle>> &daughters = particle6->children();
                for (const auto &daughter: daughters) {
                    // cout << "daughter index=" << daughter->id() << "daughter pdg_id =" << daughter->pdg_id() << endl;

                    if (daughter->pdg_id() ==-421) {
                        const auto& daughters = daughter->children();
                        if(daughters.size()==2 && (daughters[0]->pdg_id() == 321 && daughters[1]->pdg_id() == -211) || (daughters[0]->pdg_id() == -211 && daughters[1]->pdg_id() == 321)) {
                            d0bar = true;
                            d0barpx = daughter->momentum().x();
                            d0barpy = daughter->momentum().y();
                            d0barpz = daughter->momentum().z();
                            //   cout << "located a d0bar" << "index no. = " << daughter->id() << " status " << daughter->status() << endl;
                        //cout << "great news! this dbar0 decays into exactly 2 particles!" << endl;
                        //cout << daughters[0]->pdg_id() << "= daughter 1 pdg_id" << endl;
                        //cout << daughters[1]->pdg_id() << "= daughter 2 pdg_id" << endl;
                        //cout << "***************" << endl;

                            return;
                        }
                    }//if daughter is d0
                    else if (daughter->pdg_id()==-4  || daughter->status()==74|| (daughter->status()<100 && daughter->status()>80 && daughter->pdg_id()!=21) || (daughter->status()==2 && daughter->pdg_id()!=21)) tracecbardaughters(daughter, d0bar, d0barpx, d0barpy, d0barpz);
                    //else tracecbardaughters(daughter, d0bar);
                }//loop over daughters
            }

            void tracecmothers(const std::shared_ptr<const HepMC3::GenParticle>& particle6, int c_index, bool& flavor_exc) {


                if (particle6->id()==c_index) {
                    flavor_exc=true;
                    return;
                }

                if (flavor_exc) return;
                //cout << "we desire particle index no." << c_index << endl;
                //cout << "locating the mothers for particle with pdg=" << particle6->pdg_id() << "and index no." << particle6->id()<< "         ";
                const std::vector<std::shared_ptr<const HepMC3::GenParticle>> &mothers = particle6->parents();
                for (const auto &mother: mothers) {
                    //cout << "mother index=" << mother->id() << "mother pdg_id =" << mother->pdg_id() << endl;

                    if (mother->id() ==c_index) {
                        cout << "mother of c is orginal c" << endl;
                        flavor_exc=true;
                        return;
                    }//if daughter is d0
                    tracecmothers(mother, c_index, flavor_exc);
                }//loop over daughters
            }
            void tracecbarmothers(const std::shared_ptr<const HepMC3::GenParticle>& particle6, int cbar_index, bool& flavor_exc) {

                if (particle6->id()==cbar_index) {
                    flavor_exc=true;
                    return;
                }
                if(flavor_exc) return;
                //cout << "we desire particle index no." << cbar_index << endl;
                // cout << "locating the mothers for particle with pdg=" << particle6->pdg_id() << "and index no." << particle6->id()<< "         ";
                const std::vector<std::shared_ptr<const HepMC3::GenParticle>> &mothers = particle6->parents();
                for (const auto &mother: mothers) {
                    // cout << "mother index=" << mother->id() << ", mother pdg_id =" << mother->pdg_id() << endl;

                    if (mother->id() ==cbar_index) {
                        cout << "mother of cbar is orginal cbar" << endl;
                        flavor_exc=true;
                        return;
                    }//if daughter is d0
                    tracecbarmothers(mother, cbar_index, flavor_exc);
                }//loop over daughters
            }
