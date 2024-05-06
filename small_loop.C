#include <iostream>
#include <random>
#include <Pythia8/Pythia.h>
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>
#include <HepMC3/WriterAscii.h>
#include <HepMC3/ReaderAscii.h>
#include "TLorentzVector.h"
#include "TFile.h"
#include "TNtuple.h"
#include <TRandom3.h>

using namespace std;

void smearing(double& x,
              double& y,
              double& z);

void find_signal_and_process(const HepMC3::GenEvent& Event,
                             std::vector<std::vector<int>> &all_vectors,
                             std::vector<int> &SignalValues,
                             std::vector<int> &newd0,
                             std::vector<int> &newd0bar,
                             bool& D0,
                             bool& D0bar,
                             std::vector<float> &d0_p_phi_eta,
                             std::vector<float> &d0bar_p_phi_eta);

void tracedaughters(const std::shared_ptr<const HepMC3::GenParticle>& particle1,
                    bool& D0, bool& D0bar,
                    std::vector<int> &newd0,
                    std::vector<int> &newd0bar,
                    std::vector<float> &d0_p_phi_eta,
                    std::vector<float> &d0bar_p_phi_eta);

void tracemothers(const std::shared_ptr<const HepMC3::GenParticle>& InParticle2,
                  int index,
                  bool& flavor_exc);

float transition_phi(float D0del_phi);

const float upper_limit = 3.14159, lower_limit = -3.14159;//bounds to transition delta phi from (-2pi, 2pi) to (-0.5pi, 1.5*pi)
const float mass_upper_limit = 2.1, mass_lower_limit = 1.7;//bounds on mass
const double cut_off = 64.0/5.6e6; //scaling factor to reduce background, probability to be selected
const double Kmass = 0.493677, Pimass = 0.139570;//mass of kaon on pion
const double mean = 0.0; //used for smearing momentum
const double width = 0.01; //used for smearing momentum
const int d0_pid = 421, kaon_pid = 321, pion_pid = 211, gluon_pid = 21, charm_pid = 4, quarks_pid=7;//pdg_id of particles
const int gluon_splitting =1, gluon_fusion = 2, flavor_excitation_gluon = 3, flavor_excitation_quark = 4, quark_annihlation = 5;//physics processes

string get_nt_name()
//define the list of variables to be used in Ntuple
{
    string name = "a";
    name += "ifile:";
    name += "event_no:";
    name += "vertex_no:";
    name += "mass:";
    name += "phys_process:";
    name += "d0_index:";
    name += "d0bar_index:";
    name += "d0_delta_phi:";
    name += "d0_px:";
    name += "d0_py:";
    name += "d0_pz:";
    name += "d0_phi:";
    name += "d0_eta:";
    name += "d0bar_px:";
    name += "d0bar_py:";
    name += "d0bar_pz:";
    name += "d0bar_phi:";
    name += "d0bar_eta:";
    name += "TruePdg_i:";
    name += "TruePdg_j:";
    name += "iparticle_eta:";
    name += "jparticle_eta:";
    name += "AssignedPdg_i:";
    name += "AssignedPdg_j";
    return name;
}

int main(int argc, char *argv[]) {

    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " istart iend" << std::endl;
        return 1;  // Indicate an error
    }

    int istart = std::stoi(argv[1]);
    int iend = std::stoi(argv[2]);
    string name = get_nt_name();


    TString outfile = TString("btest.root");
    //TString outfile = TString::Format("/scratch/bell/awesole/analysis_hepmc_bg_cuts_05_01/ROOT/1hepmc_tree_%d_%d.root",istart, iend);

    TNtuple *T = new TNtuple("T","", name.c_str());

    TFile *results = new TFile(outfile, "recreate");
    //ifstream file_stream("pythia_output_02_19.txt");
    ifstream file_stream("D0_updates_0312.txt");

    vector<float> nt_val;

    string filename;
    bool debug = false, D0 = false, D0bar = false;
    int phys_process, fromD = 0, Phys_Process = 0;
    int k_pi_count = 0;
    float reconstructed_mass, D0_phi, D0bar_phi, delta_phi;
    float i_pT, j_pT, i_theta, j_theta, i_eta, j_eta;
    double prob;

    std::vector<std::vector<int>> all_vectors;
    all_vectors.push_back(std::vector<int>(7,0));
    std::vector<int> signalvalues = {0, 0, 0, 0, 0, 0, 0,0};
    std::vector<int> newd0 = {0, 0};
    std::vector<int> newd0bar = {0, 0};
    std::vector<float> d0_p_phi_eta = {0,0,0,0,0};
    std::vector<float> d0bar_p_phi_eta = {0,0,0,0,0};
    std::vector<float> c_p_phi_eta = {0,0,0,0,0};
    std::vector<float> cbar_p_phi_eta = {0,0,0,0,0};

    int ifile = 0;
    while (true) {
        //loop over all files in the list
        gRandom->SetSeed(ifile);

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
        cout << "file path!= " << filename.c_str() << endl;
        HepMC3::ReaderAscii reader(filename.c_str());



        while (!reader.failed()) { //loop over all events
            HepMC3::GenEvent event; //load the next event
            reader.read_event(event); //read the event
            int particlesSize = event.particles().size();
            int AssignedPid_i, AssignedPid_j, factor1, factor2;
            TLorentzVector KVector;
            TLorentzVector PiVector;
            TLorentzVector D0Vector;

            //set all vectors to contain all 0s to begin with
            all_vectors.clear();
            fill(signalvalues.begin(), signalvalues.end(), 0);
            fill(newd0.begin(), newd0.end(), 0);
            fill(newd0bar.begin(), newd0bar.end(), 0);
            all_vectors.push_back(signalvalues);

            if (event.event_number() % 1000 == 0) cout << event.event_number() << "  " << event.event_number() / 100 << "%" << endl;

            //for debugging purposes
            //if (event.event_number() != 1800 ) continue;
            //if (event.event_number() < 2949 ) continue;
            //if (event.event_number() > 2950 ) break;

            //Print information regarding each event
            if (debug) {
                cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
                cout << "Event number=" << event.event_number() << endl;
                cout << "Number of particles=" << event.particles().size() << endl;
            }

            //call the program that traces each c/cbar pair to determine if they fo to D0 and D0bar.  If they do determine the physics process:
            find_signal_and_process(event, all_vectors, signalvalues, newd0, newd0bar, D0, D0bar, d0_p_phi_eta, d0bar_p_phi_eta);

            //define delta phi and send it to the function that will transition phi from (-2pi,2pi) to (-pi/2, 1.5pi) as desired
            float D0del_phi = d0_p_phi_eta[3] - d0bar_p_phi_eta[3];
            transition_phi(D0del_phi);

            //debugging purposes
            if (event.event_number() == 29490) {
                cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
                cout << "Event number=" << event.event_number() << endl;
                cout << "Number of particles=" << event.particles().size() << endl;
                for (auto value: signalvalues) {
                    cout << "Value in signalvalues vector: " << value << endl;
                }
                for (auto value: newd0) {
                    cout << "Value in newd0 vector: " << value << endl;
                }
                for (auto value: newd0bar) {
                    cout << "Value in newd0bar vector: " << value << endl;
                }
            }

            //this begings the main data processing loop, matching each kaon or pion with the kaons and pions that follow
            for (int iParticle = 0; iParticle < particlesSize - 1; ++iParticle) { //for all particles
                const HepMC3::GenParticle &particle = *event.particles().at(iParticle); // Get the particle object
                int iParticle_pdg = particle.pdg_id();//store the pdg_id
                if (abs(iParticle_pdg) != kaon_pid && abs(iParticle_pdg) != pion_pid) continue; //only analyze kaons and pions
                for (int jParticle = iParticle + 1; jParticle < particlesSize; ++jParticle) {
                    const HepMC3::GenParticle &j_particle = *event.particles().at(jParticle); // Get the particle object
                    int jParticle_pdg = j_particle.pdg_id();//store the pdg_id

                    if (abs(jParticle_pdg) != kaon_pid && abs(jParticle_pdg) != pion_pid) continue; //only analyze kaons and pions
                    if (iParticle_pdg > 0 && jParticle_pdg > 0) continue;//ensure the particles have opposite charge
                    if (iParticle_pdg < 0 && jParticle_pdg < 0) continue;//ensure the particles have opposite charge

                    //store the momentum and energy of both particles
                    double k_x = particle.momentum().x();
                    double k_y = particle.momentum().y();
                    double k_z = particle.momentum().z();
                    double k_e = particle.momentum().e();
                    double pi_x = j_particle.momentum().x();
                    double pi_y = j_particle.momentum().y();
                    double pi_z = j_particle.momentum().z();
                    double pi_e = j_particle.momentum().e();

                    //send the momentum of the particles to be smeared by a factor of about 1%.  See details in the smearing function.
                    // we need to set a different seed for each event to ensure the smearing factor changes
                    //gRandom->SetSeed(jParticle + event.event_number());
                    smearing(k_x, k_y, k_z);
                    smearing(pi_x, pi_y, pi_z);

                    //define the energy based on momentum and smeared mass for kaon and pion
                    // define the kaon and pion 4vetors then add them to create d0 4vector
                    double Kenergy = std::sqrt(k_x * k_x + k_y * k_y + k_z * k_z + Kmass * Kmass);
                    double Pienergy = std::sqrt(pi_x * pi_x + pi_y * pi_y + pi_z * pi_z + Pimass * Pimass);
                    KVector.SetPxPyPzE(k_x, k_y, k_z, Kenergy);
                    PiVector.SetPxPyPzE(pi_x, pi_y, pi_z, Pienergy);
                    D0Vector = KVector + PiVector;
                    reconstructed_mass = D0Vector.M();

                    //assign factor 1 and 2 to match the charge of the particles, will be used later
                    if (iParticle_pdg > 0) factor1 = 1;
                    if (iParticle_pdg < 0) factor1 = -1;
                    if (jParticle_pdg > 0) factor2 = 1;
                    if (jParticle_pdg < 0) factor2 = -1;

                    //define the pT, theta and eta for each particle
                    i_pT = std::sqrt((particle.momentum().x() *particle.momentum().x()) + (particle.momentum().y()* particle.momentum().y()));//D0pT
                    j_pT = std::sqrt((j_particle.momentum().x() *j_particle.momentum().x()) + (j_particle.momentum().y()* j_particle.momentum().y()));//D0pT
                    i_theta = std::atan2(i_pT,particle.momentum().z());
                    j_theta = std::atan2(j_pT,j_particle.momentum().z());
                    i_eta = (-1 * std::log(std::tan(i_theta / 2.0)));
                    j_eta = (-1 * std::log(std::tan(j_theta / 2.0)));


                    //by default set physics process to 0 - means assume the pair is background
                    Phys_Process = 0;

                    //begin analysis of matching the particle pair to known d0 and d0bar signal pairs
                    for (const auto& sig_values : all_vectors) {

                        if (((particle.id() == sig_values[0] || particle.id() == sig_values[1]) &&
                             (j_particle.id() == sig_values[0] || j_particle.id() == sig_values[1])) ||
                            ((particle.id() == sig_values[3] || particle.id() == sig_values[4]) &&
                             (j_particle.id() == sig_values[3] || j_particle.id() == sig_values[4]))) {
                            //if the i and j particle correspond with the k-&pi+ of d0 or they correspond with k+ and pi- of D0bar

                            Phys_Process = sig_values[7]; //assign the actual physics process
                            prob = 0.00;//ensure the signal will be included in the ntuple
                            cout << "yes " << endl;
                        }
                        else prob = gRandom->Rndm(); //for each background assign each a random number
                        //cout << "prob = " << prob << endl;

                        //fill the nt_val vector that will be filled in ntuple
                        nt_val.clear();
                        nt_val.push_back(ifile);
                        nt_val.push_back(event.event_number());
                        nt_val.push_back(sig_values[6]); //vertex no
                        nt_val.push_back(reconstructed_mass);//mass
                        nt_val.push_back(Phys_Process);//physics process
                        nt_val.push_back(sig_values[2]);//D0 index
                        nt_val.push_back(sig_values[5]);//D0bar index
                        nt_val.push_back(D0del_phi);//delta_phi
                        nt_val.push_back(d0_p_phi_eta[0]);//d0 px
                        nt_val.push_back(d0_p_phi_eta[1]);//d0 py
                        nt_val.push_back(d0_p_phi_eta[2]);// d0 pz
                        nt_val.push_back(d0_p_phi_eta[3]);// d0 phi
                        nt_val.push_back(d0_p_phi_eta[4]);//d0 eta
                        nt_val.push_back(d0bar_p_phi_eta[0]);//d0b px
                        nt_val.push_back(d0bar_p_phi_eta[1]);//d0b py
                        nt_val.push_back(d0bar_p_phi_eta[2]);//d0b pz
                        nt_val.push_back(d0bar_p_phi_eta[3]);//d0b phi
                        nt_val.push_back(d0bar_p_phi_eta[4]);//d0b eta
                        nt_val.push_back(iParticle_pdg);
                        nt_val.push_back(jParticle_pdg);
                        nt_val.push_back(i_eta);
                        nt_val.push_back(j_eta);

//                        if(prob<cut_off) cout << "cut_off =" << cut_off << " prob = " << prob << endl;
                        if (reconstructed_mass <= mass_upper_limit && reconstructed_mass >= mass_lower_limit && prob < cut_off ) {//if mass is within desired range and the probability is selected
                            nt_val.push_back(factor1*kaon_pid);//assigned pdg_i
                            nt_val.push_back(factor2*pion_pid);//assigned pdg j
                            T->Fill(nt_val.data());
                            //cout << "yes 1" << endl;
                        }

                        //refine the particles with assumed opposite mass this time --generates the "swap" portion
                        Kenergy = std::sqrt(k_x * k_x + k_y * k_y + k_z * k_z + Pimass * Pimass);
                        Pienergy = std::sqrt(pi_x * pi_x + pi_y * pi_y + pi_z * pi_z + Kmass * Kmass);
                        KVector.SetPxPyPzE(k_x, k_y, k_z, Kenergy);
                        PiVector.SetPxPyPzE(pi_x, pi_y, pi_z, Pienergy);
                        D0Vector = KVector + PiVector;
                        reconstructed_mass = D0Vector.M();

                        //fill the ntuple again
                        if (reconstructed_mass <= mass_upper_limit && reconstructed_mass >= mass_lower_limit && prob < cut_off ) {
                            nt_val[3] = reconstructed_mass;
                            nt_val[nt_val.size()-2] = factor1*pion_pid;
                            nt_val[nt_val.size()-1] = factor2*kaon_pid;

                            //cout << "yes 2" << endl;
                            T->Fill(nt_val.data());
                        }
                        //nt_val.clear();
                    }
                }//for all secondary particles
            }//for all K+ pi-
        }//events - while loop

        results->cd();
        T->Write();
        //textfile.close();
        reader.close();
        ifile++;
    }//loop over all files
    cout << "~~~~~~~done~~~~~~~" << endl;
    results->Close();
    return 0;
}//main



void smearing(double& x, double& y, double& z) {
    //this function reads in the x,y,z momentum for a particle, generates a a random number according to the gaussian function and smears the value by that percentage
    //this helps model what the data will look like better
    x =  x*(1+gRandom->Gaus(mean,width));
    y =  y*(1+gRandom->Gaus(mean,width));
    z =  z*(1+gRandom->Gaus(mean,width));
}
void find_signal_and_process(const HepMC3::GenEvent& event,
                             std::vector<std::vector<int>>& all_vectors,
                             std::vector<int>& signalvalues,
                             std::vector<int> &newd0,
                             std::vector<int> &newd0bar,
                             bool& D0, bool& D0bar,
                             std::vector<float>& d0_p_phi_eta,
                             std::vector<float>& d0bar_p_phi_eta) {
    //this function is the most important one.  for each event it traces all c/cbar pairs and if they hadronize into d0 and d0bar then k/pi pairs it determines the physics process
    int c_count = 0, cbar_count = 0, phys_process;

    for (const auto &particle1: event.particles()) { //keep a count of c/cbar in event
        if (particle1->pdg_id() == charm_pid) c_count++;
        if (particle1->pdg_id() == -charm_pid) cbar_count++;
    }

    if (c_count > 0 && cbar_count == c_count) {//begin analysis only for events with at least 1 c/cbar pair
        D0 = false;
        D0bar = false;

        all_vectors.clear();
        d0bar_p_phi_eta.clear();
        d0_p_phi_eta.clear();
        for (const auto &vertex: event.vertices()) {//for all vertices in the event
            bool flavor_exc = false;


            //define 2 vectors - one for the incoming and one for the outgoing particles of each event
            const std::vector<std::shared_ptr<const HepMC3::GenParticle>>& particlesOut = vertex->particles_out();
            const std::vector<std::shared_ptr<const HepMC3::GenParticle>>& particlesIn = vertex->particles_in();

            int count4 = 0;
            int count_4 = 0;

            //Loop over all outgoing particles and count the number of outgoing c and cbar
            for (const auto &particle: particlesOut) {
                if (particle->pdg_id() == charm_pid) count4++;
                if (particle->pdg_id() == -charm_pid) count_4++;
            }

            //begin analysis for vertices with at least 1 outgoing c/cbar pairs
            if (count4 > 0 && count4 == count_4) {

                //variables initialization
                D0 = false;
                D0bar = false;
                newd0.clear();
                newd0bar.clear();

                for (const auto &particle3: particlesOut) {//loop over the outgoing c and cbar
                    //send the c and cbar particles to the trace daughters function.  see function for details
                    //if outgoing particle is not c/cbar skip
                    if (particle3->pdg_id() == charm_pid) tracedaughters(particle3, D0, D0bar, newd0, newd0bar, d0_p_phi_eta, d0bar_p_phi_eta);
                    if (particle3->pdg_id() == -charm_pid) tracedaughters(particle3, D0, D0bar, newd0, newd0bar, d0_p_phi_eta, d0bar_p_phi_eta);
                    else continue;
                }

                if (D0 and D0bar) {//if the vertex leads to a D0 and D0bar that both decay into K/Pi pairs
                    signalvalues.clear();
                    cout << "*^*^*^*^*^*^*" << endl;
                    cout << "signal!!!" << endl;

                    for (int id : newd0){
                        //fill the signval values with the d0, k-, pi+ indices
                        signalvalues.push_back(id);
                    }
                    for (int id : newd0bar){
                        //fill the signval values with the d0bar, k+, pi- indices
                        signalvalues.push_back(id);
                    }
                    //and fill the vertex_id
                    signalvalues.push_back(vertex->id());

                    //begin analysis to determine the physics process of creating c/cbar
                    for (const auto &particle3: particlesOut) {
                        if (particle3->status() < 60 && particle3->status() > 50) {//if the c/cbar have pythia code inidacting created during final state showers
                            phys_process = -1; //set an arbitrary physics process to detect errors
                            for (const auto &particle4: particlesIn) {
                                //if the incoming particles of the c/cbar vertex contain 1 gluon, declare gluon splitting
                                if (particlesIn.size() == 1 && particle4->pdg_id() == gluon_pid) phys_process = gluon_splitting;//gluon splitting
                            }
                        }
                        else if (particle3->status() <50) {//only for particles from initial state showers, gluon fusion and flavor excitation
                            phys_process = -1;//to note errors
                            for (const auto &particle4: particlesIn) {//if both d0 and d0bar, begin analysis of incoming particles of vertex

                                //if there are 2 incoming particles and they are gluons, gluon fusion, they are quarks, quark annihlation
                                if (particlesIn.size() == 2 && particle4->pdg_id() == gluon_pid) phys_process = gluon_fusion;//gluon fusion
                                if (particlesIn.size() == 2 && particle4->pdg_id() < quarks_pid) phys_process = quark_annihlation; //quark annihlation

                                //if one gluon comes in and the particles are initial state showers begin analysis for flavor excitation
                                if (particlesIn.size() == 1 && particle4->pdg_id() == gluon_pid) {//begin analysis to determine flavor excitation
                                    int gluon_in = 0,  c_in = 0, c_out = 0, cbar_in = 0, cbar_out = 0, q_in=0;

                                    for (const auto &vertex: event.vertices()) {//loop over all vertices to only select the ones with g+c comming in and g+c coming out
                                        for (const auto &inparticle: vertex->particles_in()) {//counts incoming gluons and c/cbar
                                            if (inparticle->pdg_id() == gluon_pid) gluon_in++;
                                            if (inparticle->pdg_id() == charm_pid) c_in++;
                                            if (inparticle->pdg_id() == -charm_pid) cbar_in++;
                                            if (inparticle->pdg_id()< quarks_pid && abs(inparticle->pdg_id())!=charm_pid) q_in++;
                                        }//counts incoming gluons and c/cbar
                                        for (const auto &outparticle: vertex->particles_out()) {//counts outgoing gluons and c/cbar
                                            if (outparticle->pdg_id() == charm_pid) c_out++;
                                            if (outparticle->pdg_id() == -charm_pid) cbar_out++;
                                        }//counts outgoing gluons and c/cbar

                                        //cout <<"particles in =" << endl;
                                        for (const auto &inparticle2: vertex->particles_in()) { //analyze vertices with 1 charm & 1 gluon or quark coming in and one charm going out
                                            //cout << "particle pdg = " << inparticle2->pdg_id() << endl;
                                            if (abs(inparticle2->pdg_id()) == charm_pid && gluon_in == 1 and c_in == 1 and c_out == 1) {
                                                int c_index = inparticle2->id();
                                                tracemothers(inparticle2, c_index, flavor_exc); //see trace mothers for more details
                                            }
                                            else if (abs(inparticle2->pdg_id()) == charm_pid && q_in == 1 and c_in == 1 and c_out == 1) {
                                                int c_index = inparticle2->id();
                                                tracemothers(inparticle2, c_index, flavor_exc); //see tracemothers for more details
                                            }
                                        }
                                        if (flavor_exc && gluon_in ==1) phys_process = flavor_excitation_gluon; //flavor excitation
                                        //cout << "The end " << endl;
                                        if (flavor_exc && q_in ==1) phys_process = flavor_excitation_quark; //flavor excitation
                                        flavor_exc = false;
                                        gluon_in = 0;
                                        q_in = 0;
                                        c_in = 0;
                                        c_out = 0;
                                        cbar_in = 0;
                                        cbar_out = 0;
                                    }//for all vertices
                                }//flavor excitation
                            }//incoming particles
                        }//inital state
                    }
                    cout << "phys process = " << phys_process << endl;
                    signalvalues.push_back(phys_process);//add physics process
                    all_vectors.push_back(signalvalues); //add signalvalues to all_vectors
                    fill(signalvalues.begin(), signalvalues.end(), 0);
                }//if D0 and D0bar ** add analysis here
            }//only for vertices containing at least 1 c/cbar pair
        }//for all vertices
    }//only for events with same number of c/cbar
    if (all_vectors.size() == 0) {
        fill(signalvalues.begin(), signalvalues.end(), 0);
        all_vectors.push_back(signalvalues);
    }
    if (d0_p_phi_eta.size() == 0) fill(d0_p_phi_eta.begin(), d0_p_phi_eta.end(), 0);
    if (d0bar_p_phi_eta.size() == 0) fill(d0bar_p_phi_eta.begin(), d0bar_p_phi_eta.end(), 0);
    return;
}//trace for signal

void tracedaughters(const std::shared_ptr<const HepMC3::GenParticle>& particle6, bool& d0, bool& d0bar,  std::vector<int> &newd0, std::vector<int> &newd0bar,  std::vector<float> &d0_p_phi_eta,  std::vector<float> &d0bar_p_phi_eta){
    //this is a recursive function that traces the daughters of each c/cbar
    const std::vector<std::shared_ptr<const HepMC3::GenParticle>> &daughters = particle6->children();//create a vector of all the daughers of the particle

    for (const auto &daughter: daughters) {//analyze each daughter
        if (abs(daughter->pdg_id()) ==d0_pid) {//if the daughter is d0/d0bar
            const auto& daughters = daughter->children();//anaylze the d0 daughters
            if (daughters.size() == 2 && (abs(daughters[0]->pdg_id()) == kaon_pid && abs(daughters[1]->pdg_id()) == pion_pid)) {//if k then pi
                if (daughter->pdg_id()==d0_pid) {
                    d0 = true;
                    newd0.push_back(daughters[0]->id());//pushback index of k-
                    newd0.push_back(daughters[1]->id());//pushback index of pi+
                    newd0.push_back(daughter->id());//push back d0 index
                    d0_p_phi_eta.push_back(daughter->momentum().x());
                    d0_p_phi_eta.push_back(daughter->momentum().y());
                    d0_p_phi_eta.push_back(daughter->momentum().z());

                    d0_p_phi_eta.push_back(std::atan2(daughter->momentum().y(), daughter->momentum().x()));
                    float D0_pT = std::sqrt((daughter->momentum().x() *daughter->momentum().x()) + (daughter->momentum().y()* daughter->momentum().y()));//D0pT
                    float D0_theta = std::atan2(D0_pT,daughter->momentum().z());
                    d0_p_phi_eta.push_back(-1 * std::log(std::tan(D0_theta / 2.0)));

                    return;
                }
                if (daughter->pdg_id() == -d0_pid) {//same for d0bar
                    d0bar = true;
                    newd0bar.push_back(daughters[0]->id());//pushback id of k+
                    newd0bar.push_back(daughters[1]->id());//pushback id of pi-
                    newd0bar.push_back(daughter->id());
                    d0bar_p_phi_eta.push_back(daughter->momentum().x());
                    d0bar_p_phi_eta.push_back(daughter->momentum().y());
                    d0bar_p_phi_eta.push_back(daughter->momentum().z());

                    d0bar_p_phi_eta.push_back(std::atan2(daughter->momentum().y(), daughter->momentum().x()));
                    float D0_pT = std::sqrt((daughter->momentum().x() *daughter->momentum().x() ) + ( daughter->momentum().y()* daughter->momentum().y()));//D0pT
                    float D0_theta = std::atan2(D0_pT,daughter->momentum().z() );
                    d0bar_p_phi_eta.push_back(-1 * std::log(std::tan(D0_theta / 2.0)));
                    return;
                }
            }//if k then pi
            else if (daughters.size() == 2 && (abs(daughters[0]->pdg_id()) == pion_pid && abs(daughters[1]->pdg_id()) == kaon_pid)){//if pi then k
                if (daughter->pdg_id()==d0_pid) {
                    d0 = true;
                    newd0.push_back(daughters[1]->id());//pushback id of k-
                    newd0.push_back(daughters[0]->id());//pushback id of pi+
                    newd0.push_back(daughter->id());
                    d0_p_phi_eta.push_back(daughter->momentum().x());
                    d0_p_phi_eta.push_back(daughter->momentum().y());
                    d0_p_phi_eta.push_back(daughter->momentum().z());

                    d0_p_phi_eta.push_back(std::atan2(daughter->momentum().y(), daughter->momentum().x()));
                    float D0_pT = std::sqrt((daughter->momentum().x() *daughter->momentum().x() ) + ( daughter->momentum().y()* daughter->momentum().y()));//D0pT
                    float D0_theta = std::atan2(D0_pT,daughter->momentum().z() );
                    d0_p_phi_eta.push_back(-1 * std::log(std::tan(D0_theta / 2.0)));
                }
                if (daughter->pdg_id() == -d0_pid) {
                    d0bar = true;
                    newd0bar.push_back(daughters[1]->id());//pushback id of k+
                    newd0bar.push_back(daughters[0]->id());//pushback id of pi-
                    d0bar_p_phi_eta.push_back(daughter->id());
                    d0bar_p_phi_eta.push_back(daughter->momentum().x());
                    d0bar_p_phi_eta.push_back(daughter->momentum().y());
                    d0bar_p_phi_eta.push_back(daughter->momentum().z());

                    d0bar_p_phi_eta.push_back(std::atan2(daughter->momentum().y(), daughter->momentum().x()));
                    float D0_pT = std::sqrt((daughter->momentum().x() *daughter->momentum().x() ) + ( daughter->momentum().y()* daughter->momentum().y()));//D0pT
                    float D0_theta = std::atan2(D0_pT,daughter->momentum().z() );
                    d0bar_p_phi_eta.push_back(-1 * std::log(std::tan(D0_theta / 2.0)));
                }
                return;
            }
        }
        //if the daughter is not d0 trace the next daughter with a few exceptions mainly we just want to trace the charm quarks and any diqaurks
        else if (abs(daughter->pdg_id())==charm_pid || daughter->status()==74|| (daughter->status()<100 && daughter->status()>80 && daughter->pdg_id()!=gluon_pid) || (daughter->status()==2 && daughter->pdg_id()!=gluon_pid)) tracedaughters(daughter, d0, d0bar, newd0, newd0bar, d0_p_phi_eta, d0bar_p_phi_eta);
    }//loop over daughters

}
void tracemothers(const std::shared_ptr<const HepMC3::GenParticle>& particle6, int c_index, bool& flavor_exc) {
    //this function takes an incoming charm quark and traces its mothers.  if the desired mother charm is reached we denote it as flavor excitation
    if (particle6->id()==c_index) {
        flavor_exc=true;
        return;
    }

    if (flavor_exc) return;
    const std::vector<std::shared_ptr<const HepMC3::GenParticle>> &mothers = particle6->parents();
    for (const auto &mother: mothers) {
        if (mother->id() ==c_index) {
            flavor_exc=true;
            return;
        }//if daughter is d0
        tracemothers(mother, c_index, flavor_exc);
    }//loop over daughters
}

float transition_phi(float D0del_phi){
//this function reads in an angle d0del_phi and changes it so that the range is confined to (-pi/2 to 3pi/2)
    if (D0del_phi > upper_limit) {
        D0del_phi = D0del_phi - 2 * TMath::Pi();
    }
    if (D0del_phi < lower_limit) {
        D0del_phi = 2 * TMath::Pi() + D0del_phi;
    }
    if (D0del_phi < -0.5 * TMath::Pi() && D0del_phi > -1 * TMath::Pi()) {
        D0del_phi = 2 * TMath::Pi() + D0del_phi;
    }
    return D0del_phi;
}
