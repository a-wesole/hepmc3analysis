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
#include <TTree.h>
#include <TBranch.h>


using namespace std;

void smearing(double& x,
              double& y,
              double& z);

void find_signal_and_process(const HepMC3::GenEvent& Event,
                             std::vector<std::array<int, 8>>& all_vectors,
                             std::array<int,8> SignalValues,
                             std::vector<int>& newd0,
                             std::vector<int>& newd0bar,
                             bool& D0,
                             bool& D0bar,
                             std::vector<float>& d0_p_phi_eta,
                             std::vector<float>& d0bar_p_phi_eta,
                             bool& D0_notsignal, bool& D0bar_notsignal);

void tracedaughters(const std::shared_ptr<const HepMC3::GenParticle>& particle1,
                    bool& D0, bool& D0bar,
                    std::vector<int> &newd0,
                    std::vector<int> &newd0bar,
                    std::vector<float> &d0_p_phi_eta,
                    std::vector<float> &d0bar_p_phi_eta,
                    bool& D0_notsignal, bool& D0bar_notsignal);

void tracemothers(const std::shared_ptr<const HepMC3::GenParticle>& InParticle2,
                  int index,
                  bool& flavor_exc);

float transition_phi(float& D0del_phi);

void check_for_duplicates(std::vector<float> vertex_no,
                           std::vector<float> d0_index,
                           std::vector<float> d0bar_index,
                           std::vector<float> &phy_process);

const float upper_limit = 3.14159, lower_limit = -3.14159;//bounds to transition delta phi from (-2pi, 2pi) to (-0.5pi, 1.5*pi)
const float mass_upper_limit = 2.1, mass_lower_limit = 1.6;//bounds on mass
//const double cut_off = 1000.0/5.6e6; //scaling factor to reduce background, probability to be selected
const double cut_off = 1.00; //scaling factor to reduce background, probability to be selected
const double Kmass = 0.493677, Pimass = 0.139570;//mass of kaon on pion
const double mean = 0.0; //used for smearing momentum
const double width = 0.01; //used for smearing momentum
const int d0_pid = 421, kaon_pid = 321, pion_pid = 211, gluon_pid = 21, charm_pid = 4, quarks_pid=7;//pdg_id of particles
const int gluon_splitting =1, gluon_fusion = 2, flavor_excitation_gluon = 3, flavor_excitation_quark = 4, quark_annihlation = 5;//physics processes

float Ifile, event_no;

std::vector<float> vertex_no;
std::vector<float> mass;
std::vector<float> phy_process;
std::vector<float> d0_index;
std::vector<float> d0bar_index;
std::vector<float> phi1;
std::vector<float> prob1;
std::vector<float> d0_px;
std::vector<float> d0_py;
std::vector<float> d0_pz;
std::vector<float> d0_phi;
std::vector<float> d0_eta;
std::vector<float> d0bar_px;
std::vector<float> d0bar_py;
std::vector<float> d0bar_pz;
std::vector<float> d0bar_phi;
std::vector<float> d0bar_eta;
std::vector<float> TruePdg_i;
std::vector<float> TruePdg_j;
std::vector<float> iparticle_eta;
std::vector<float> jparticle_eta;
std::vector<float> iparticle_px;
std::vector<float> iparticle_py;
std::vector<float> iparticle_pz;
std::vector<float> jparticle_px;
std::vector<float> jparticle_py;
std::vector<float> jparticle_pz;
std::vector<float> AssignedPdg_i;
std::vector<float> AssignedPdg_j;
std::vector<float> fromD0;

void set_branches(TTree *tree)
// define the list of variables to be set as branches of tree
{

    tree->Branch("Ifile", &Ifile, "Ifile/F");
    tree->Branch("event_no", &event_no, "event_no/F");

    tree->Branch("vertex_no", &vertex_no);
    tree->Branch("mass", &mass);
    tree->Branch("phy_process", &phy_process);
    tree->Branch("d0_index", &d0_index);
    tree->Branch("d0bar_index", &d0bar_index);
    tree->Branch("phi", &phi1);
    tree->Branch("prob", &prob1);

    tree->Branch("d0_px", &d0_px);
    tree->Branch("d0_py", &d0_py);
    tree->Branch("d0_pz", &d0_pz);
    tree->Branch("d0_phi", &d0_phi);
    tree->Branch("d0_eta", &d0_eta);

    tree->Branch("d0bar_px", &d0bar_px);
    tree->Branch("d0bar_py", &d0bar_py);
    tree->Branch("d0bar_pz", &d0bar_pz);
    tree->Branch("d0bar_phi", &d0bar_phi);
    tree->Branch("d0bar_eta", &d0bar_eta);
    
    tree->Branch("TruePdg_i", &TruePdg_i);
    tree->Branch("TruePdg_j", &TruePdg_j);
    tree->Branch("iparticle_eta", &iparticle_eta);
    tree->Branch("jparticle_eta", &jparticle_eta);
    tree->Branch("iparticle_px", &iparticle_px);
    tree->Branch("iparticle_py", &iparticle_py);
    tree->Branch("iparticle_pz", &iparticle_pz);
    tree->Branch("jparticle_py", &jparticle_py);
    tree->Branch("jparticle_px", &jparticle_px);
    tree->Branch("jparticle_pz", &jparticle_pz);
    
    tree->Branch("AssignedPdg_i", &AssignedPdg_i);
    tree->Branch("AssignedPdg_j", &AssignedPdg_j);
    tree->Branch("fromD0", &fromD0);
}

int main(int argc, char *argv[]) {

    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " istart iend" << std::endl;
        return 1;  // Indicate an error
    }

    int istart = std::stoi(argv[1]);
    int iend = std::stoi(argv[2]);
    gRandom->SetSeed(istart);


    //TString outfile = TString("/scratch/bell/awesole/atest.root");
    //TString outfile = TString::Format("/scratch/bell/awesole/analysis_hepmc_bg_cuts_06_18/ROOT/small_list_hepmc_tree_%d_%d.root",istart, iend);
    //TString outfile = TString::Format("/scratch/bell/awesole/analysis_hepmc_bg_cuts_06_18/ROOT/small_list_hepmc_tree_%d_%d.root",istart, iend);
    //TString outfile = TString::Format("/scratch/bell/awesole/analysis_hepmc_07_15_orig_prod/ROOT/hepmc_tree_%d_%d.root",istart, iend);
    //TString outfile = TString::Format("/scratch/bell/awesole/analysis_hepmc_07_18_updatedBkg/ROOT/July29updated_hepmc_tree_%d_%d.root",istart, iend);
    TString outfile = TString("test_1file.root");


    //TNtuple *T = new TNtuple("T","", name.c_str());
    TTree *event_tree = new TTree("event_tree","");
    set_branches(event_tree);

    TFile *results = new TFile(outfile, "recreate");
    //ifstream file_stream("D0_updates_0312.txt");
    ifstream file_stream("new_production.list");
    //ifstream file_stream("sorted.list");

    //ifstream file_stream("test.list");

    vector<float> nt_val;

    string filename;
    bool debug = true, D0 = false, D0bar = false, D0_notsignal=false, D0bar_notsignal = false;
    int phys_process, Phys_Process = 0;
    int k_pi_count = 0;
    float reconstructed_mass, D0_phi, D0bar_phi, delta_phi, phi;
    float i_pT, j_pT, i_theta, j_theta, i_eta, j_eta, fromD, fromD0_index;
    double prob;

    std::vector<std::array<int, 8>> all_vectors;

    //    std::vector<std::vector<int>> all_vectors;
    std::array<int, 8> signalvalues = {0,0,0,0,0,0,0,0} ;
    all_vectors.push_back(signalvalues);
    std::vector<int> newd0 = {0, 0};
    std::vector<int> newd0bar = {0, 0};
    std::vector<float> d0_p_phi_eta = {0,0,0,0,0};
    std::vector<float> d0bar_p_phi_eta = {0,0,0,0,0};
    std::vector<float> c_p_phi_eta = {0,0,0,0,0};
    std::vector<float> cbar_p_phi_eta = {0,0,0,0,0};

    int ifile = 0;
    while (true) {
        //loop over all files in the list

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

        //TFile fin(filename.c_str());

        /*
        if (fin.IsZombie())
        {
            ifile++; // you need this here to keep track of the file
            continue;
        }
        */

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
            for (int i = 0; i < 8; ++i)
            {
                signalvalues[i] = 0;
            }
            fill(newd0.begin(), newd0.end(), 0);
            fill(newd0bar.begin(), newd0bar.end(), 0);
            all_vectors.push_back(signalvalues);

            if (event.event_number() % 1000 == 0) cout << event.event_number() << "  " << event.event_number() / 100 << "%" << endl;

            //for debugging purposes
            //if (event.event_number() != 13 ) continue;
            //if (event.event_number() < 2949 ) continue;
            //if (event.event_number() < 9000) continue;

            //Print information regarding each event
            if (debug) {
                cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
                cout << "Event number=" << event.event_number() << endl;
                cout << "Number of particles=" << event.particles().size() << endl;
            }

            //call the program that traces each c/cbar pair to determine if they fo to D0 and D0bar.  If they do determine the physics process:
            find_signal_and_process(event, all_vectors, signalvalues, newd0, newd0bar, D0, D0bar, d0_p_phi_eta, d0bar_p_phi_eta, D0_notsignal, D0bar_notsignal);

            //define delta phi and send it to the function that will transition phi from (-2pi,2pi) to (-pi/2, 1.5pi) as desired

            //debugging purposes
            if (event.event_number() == 156471532476) {
                cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
                cout << "Event number=" << event.event_number() << endl;
                cout << "Number of particles=" << event.particles().size() << endl;
            }
            Ifile = ifile;
            event_no = event.event_number();
            //if (event_no > 0) break;

            if(vertex_no.size()>0) vertex_no.clear();
            if(mass.size()>0) mass.clear();
            if(phy_process.size()>0) phy_process.clear();
            if(d0_index.size()>0) d0_index.clear();
            if(d0bar_index.size()>0) d0bar_index.clear();
            if(phi1.size()>0) phi1.clear();
            if(prob1.size()>0) prob1.clear();
            if(d0_px.size()>0) d0_px.clear();
            if(d0_py.size()>0) d0_py.clear();
            if(d0_pz.size()>0) d0_pz.clear();
            if(d0_phi.size()>0) d0_phi.clear();
            if(d0_eta.size()>0) d0_eta.clear();
            if(d0bar_px.size()>0) d0bar_px.clear();
            if(d0bar_py.size()>0) d0bar_py.clear();
            if(d0bar_pz.size()>0) d0bar_pz.clear();
            if(d0bar_phi.size()>0) d0bar_phi.clear();
            if(d0bar_eta.size()>0) d0bar_eta.clear();
            if(TruePdg_i.size()>0) TruePdg_i.clear();
            if(TruePdg_j.size()>0) TruePdg_j.clear();
            if(iparticle_eta.size()>0) iparticle_eta.clear();
            if(jparticle_eta.size()>0) jparticle_eta.clear();
            if(iparticle_px.size()>0) iparticle_px.clear();
            if(iparticle_py.size()>0) iparticle_py.clear();
            if(iparticle_pz.size()>0) iparticle_pz.clear();
            if(jparticle_px.size()>0) jparticle_px.clear();
            if(jparticle_py.size()>0) jparticle_py.clear();
            if(jparticle_pz.size()>0) jparticle_pz.clear();
            if(AssignedPdg_i.size()>0) AssignedPdg_i.clear();
            if(AssignedPdg_j.size()>0) AssignedPdg_j.clear();
            if(fromD0.size()>0) fromD0.clear();

            // this begings the main data processing loop, matching each kaon or pion with the kaons and pions that follow
            for (int iParticle = 0; iParticle < particlesSize - 1; ++iParticle)
            {                                                                           // for all particles
                const HepMC3::GenParticle &particle = *event.particles().at(iParticle); // Get the particle object
                int iParticle_pdg = particle.pdg_id();//store the pdg_id
                if (abs(iParticle_pdg) != kaon_pid && abs(iParticle_pdg) != pion_pid) continue; //only analyze kaons and pions
                for (int jParticle = iParticle + 1; jParticle < particlesSize; ++jParticle) {
                    const HepMC3::GenParticle &j_particle = *event.particles().at(jParticle); // Get the particle object
                    int jParticle_pdg = j_particle.pdg_id();//store the pdg_id

                    if (abs(jParticle_pdg) != kaon_pid && abs(jParticle_pdg) != pion_pid) continue; //only analyze kaons and pions
                    if (iParticle_pdg > 0 && jParticle_pdg > 0) continue;//ensure the particles have opposite charge
                    if (iParticle_pdg < 0 && jParticle_pdg < 0) continue;//ensure the particles have opposite charge


                    fromD = 0.0;
                    fromD0_index = -999;
                    const std::vector<std::shared_ptr<const HepMC3::GenParticle>> &iparticle_mothers = particle.parents();
                    const std::vector<std::shared_ptr<const HepMC3::GenParticle>> &jparticle_mothers = j_particle.parents();
                    for (const auto &iparticle_mother: iparticle_mothers) {
                        if (abs(iparticle_mother->pdg_id()) == d0_pid) {
                            fromD0_index = iparticle_mother->id();
                           }//if daughter is d0
                    }
                    for (const auto &jparticle_mother: jparticle_mothers){
                        if(jparticle_mother->id() == fromD0_index) {
                            fromD = fromD0_index;
                            const std::vector<std::shared_ptr<const HepMC3::GenParticle>> &D0_mothers =jparticle_mother->parents();
                            for (const auto &D0particle_mother : D0_mothers) {
                                if (abs(int(D0particle_mother->pdg_id() / 100) %10) == 5 || abs(int(D0particle_mother->pdg_id()/1000)%10)==5) fromD = -101; //denoting nonprompt D0
                            }
                            //cout << "this should be a D0....pid=" << jparticle_mother->pdg_id() << endl;
                            break;
                        }
                    }
                    if(abs(iParticle_pdg)==abs(jParticle_pdg)) fromD=0.0;
                    //if (fromD!=0.0) cout << "located a k/pi that goes to D0 ---- d0 id is " << fromD0_index << endl;


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
                    smearing(k_x, k_y, k_z);
                    smearing(pi_x, pi_y, pi_z);

                    //define the energy based on momentum and smeared mass for kaon and pion
                    // define the kaon and pion 4vetors then add them to create d0 4vector
                    double Kenergy = std::sqrt(k_x * k_x + k_y * k_y + k_z * k_z + Kmass * Kmass);
                    double Pienergy = std::sqrt(pi_x * pi_x + pi_y * pi_y + pi_z * pi_z + Pimass * Pimass);
                    KVector.SetPxPyPzE(k_x, k_y, k_z, Kenergy);
                    PiVector.SetPxPyPzE(pi_x, pi_y, pi_z, Pienergy);
                    D0Vector = KVector + PiVector;
                    phi = std::atan2(D0Vector.Py(),D0Vector.Px());
                    transition_phi(phi);
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



                    //begin analysis of matching the particle pair to known d0 and d0bar signal pairs
                    for (const auto& sig_values : all_vectors) {
                    //by default set physics process to 0 - means assume the pair is background
                        Phys_Process = 0;
                        //if (sig_values[7] == 101) cout << "event=" << event_no << " ifile=" << Ifile << " vertex = " << sig_values[6] << endl;
                        //if(D0_no and !D0bar) Phys_Process=101;//Update as of July 18.  Previous categorization did not in linear backgournd, pileup at D0mass becuase all D0 candidates nt directly signal were listed as background. necessary cha
                        //if(!D0 and D0bar_no) Phys_Process=101;//if there is a nonsignal D0 and there is the other is not a signal D0bar
                        //if (reconstructed_mass > 1.82 && reconstructed_mass < 1.90 && D0 && !D0bar) cout << "D0 out only phy_process =" << Phys_Process << endl; 
                        //if (reconstructed_mass > 1.82 && reconstructed_mass < 1.90 && !D0 && D0bar) cout << "D0BAR out only phy_process =" << Phys_Process << endl; 

                        if((particle.id() == sig_values[0] && j_particle.id() == sig_values[1]) || 
                           (particle.id() == sig_values[1] && j_particle.id() == sig_values[0]) ||
                           (particle.id() == sig_values[3] && j_particle.id() == sig_values[4]) ||
                           (particle.id() == sig_values[4] && j_particle.id() == sig_values[3])) {

                       /*
                       if (((particle.id() == sig_values[0] || particle.id() == sig_values[1]) &&
                             (j_particle.id() == sig_values[0] || j_particle.id() == sig_values[1])) ||
                            ((particle.id() == sig_values[3] || particle.id() == sig_values[4]) &&
                             (j_particle.id() == sig_values[3] || j_particle.id() == sig_values[4]))) {
                            //if the i and j particle correspond with the k-&pi+ of d0 or they correspond with k+ and pi- of D0bar
                            */

                            Phys_Process = sig_values[7]; //assign the actual physics process
                            prob = 0.00;//ensure the signal will be included in the ntuple
                            cout << "yes " << endl;
                            //cout << "updated phy_process==" << Phys_Process << endl;
                        }
                        else prob = gRandom->Rndm(); //for each background assign each a random number
                        //cout << "prob = " << prob << endl;

                        //fill the nt_val vector that will be filled in ntuple
                        
                        if (reconstructed_mass <= mass_upper_limit && reconstructed_mass >= mass_lower_limit && prob < cut_off ) {//if mass is within desired range and the probability is selected
                        vertex_no.push_back(sig_values[6]); //vertex no
                        mass.push_back(reconstructed_mass);//mass
                        phy_process.push_back(Phys_Process);//physics process
                        d0_index.push_back(sig_values[2]);//D0 index
                        d0bar_index.push_back(sig_values[5]);//D0bar index
                        phi1.push_back(phi);//delta_phi
                        prob1.push_back(prob);
                        //cout << "phi =" << D0del_phi << endl;
                        d0_px.push_back(d0_p_phi_eta[0]);//d0 px
                        d0_py.push_back(d0_p_phi_eta[1]);//d0 py
                        d0_pz.push_back(d0_p_phi_eta[2]);// d0 pz
                        d0_phi.push_back(d0_p_phi_eta[3]);// d0 phi
                        d0_eta.push_back(d0_p_phi_eta[4]);//d0 eta
                        d0bar_px.push_back(d0bar_p_phi_eta[0]);//d0b px
                        d0bar_py.push_back(d0bar_p_phi_eta[1]);//d0b py
                        d0bar_pz.push_back(d0bar_p_phi_eta[2]);//d0b pz
                        d0bar_phi.push_back(d0bar_p_phi_eta[3]);//d0b phi
                        d0bar_eta.push_back(d0bar_p_phi_eta[4]);//d0b eta
                        TruePdg_i.push_back(iParticle_pdg);
                        TruePdg_j.push_back(jParticle_pdg);
                        iparticle_eta.push_back(i_eta);
                        jparticle_eta.push_back(j_eta);
                        iparticle_px.push_back(k_x);
                        iparticle_py.push_back(k_y);
                        iparticle_pz.push_back(k_z);
                        jparticle_px.push_back(pi_x);
                        jparticle_py.push_back(pi_y);
                        jparticle_pz.push_back(pi_z);

                        AssignedPdg_i.push_back(factor1 * kaon_pid); // assigned pdg_i
                        AssignedPdg_j.push_back(factor2 * pion_pid); // assigned pdg j
                        fromD0.push_back(fromD);


                        }

                        //refine the particles with assumed opposite mass this time --generates the "swap" portion
                        Kenergy = std::sqrt(k_x * k_x + k_y * k_y + k_z * k_z + Pimass * Pimass);
                        Pienergy = std::sqrt(pi_x * pi_x + pi_y * pi_y + pi_z * pi_z + Kmass * Kmass);
                        KVector.SetPxPyPzE(k_x, k_y, k_z, Kenergy);
                        PiVector.SetPxPyPzE(pi_x, pi_y, pi_z, Pienergy);
                        D0Vector = KVector + PiVector;
                        phi = std::atan2(D0Vector.Py(), D0Vector.Px());
                        transition_phi(phi);
                        reconstructed_mass = D0Vector.M();

                        if (reconstructed_mass <= mass_upper_limit && reconstructed_mass >= mass_lower_limit && prob < cut_off ) {
                        vertex_no.push_back(sig_values[6]); //vertex no
                        //cout << "vertex id = " << sig_values[6] << endl;
                        mass.push_back(reconstructed_mass);//mass
                        phy_process.push_back(Phys_Process);//physics process
                        d0_index.push_back(sig_values[2]);//D0 index
                        d0bar_index.push_back(sig_values[5]);//D0bar index
                        phi1.push_back(phi);//delta_phi
                        prob1.push_back(prob);
                        d0_px.push_back(d0_p_phi_eta[0]);//d0 px
                        d0_py.push_back(d0_p_phi_eta[1]);//d0 py
                        d0_pz.push_back(d0_p_phi_eta[2]);// d0 pz
                        d0_phi.push_back(d0_p_phi_eta[3]);// d0 phi
                        d0_eta.push_back(d0_p_phi_eta[4]);//d0 eta
                        d0bar_px.push_back(d0bar_p_phi_eta[0]);//d0b px
                        d0bar_py.push_back(d0bar_p_phi_eta[1]);//d0b py
                        d0bar_pz.push_back(d0bar_p_phi_eta[2]);//d0b pz
                        d0bar_phi.push_back(d0bar_p_phi_eta[3]);//d0b phi
                        d0bar_eta.push_back(d0bar_p_phi_eta[4]);//d0b eta
                        TruePdg_i.push_back(iParticle_pdg);
                        TruePdg_j.push_back(jParticle_pdg);
                        iparticle_eta.push_back(i_eta);
                        jparticle_eta.push_back(j_eta);
                        iparticle_px.push_back(k_x);
                        iparticle_py.push_back(k_y);
                        iparticle_pz.push_back(k_z);
                        jparticle_px.push_back(pi_x);
                        jparticle_py.push_back(pi_y);
                        jparticle_pz.push_back(pi_z);
                        AssignedPdg_i.push_back(factor1 * pion_pid);
                        AssignedPdg_j.push_back(factor2 * kaon_pid);
                        fromD0.push_back(fromD);

                        }

                        if (vertex_no.size()>=2) check_for_duplicates(vertex_no, d0_index, d0bar_index, phy_process);

                        //nt_val.clear();
                    }
                }//for all secondary particles
            } // for all K+ pi-
            if(mass.size()>0) event_tree->Fill();
        }//events - while loop

        results->cd();
        event_tree->Write();
        //textfile.close();
        reader.close();
        ifile++;
    }//loop over all files
    cout << "~~~~~~~done~~~~~~~" << endl;
    cout << " " << endl;
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
                             std::vector<std::array<int, 8>>& all_vectors,
                             std::array<int,8> signalvalues,
                             std::vector<int> &newd0,
                             std::vector<int> &newd0bar,
                             bool& D0, bool& D0bar,
                             std::vector<float>& d0_p_phi_eta,
                             std::vector<float>& d0bar_p_phi_eta,
                             bool& D0_notsignal, bool& D0bar_notsignal) {
    //this function is the most important one.  for each event it traces all c/cbar pairs and if they hadronize into d0 and d0bar then k/pi pairs it determines the physics process
    int c_count = 0, cbar_count = 0, phys_process;

    for (const auto &particle1: event.particles()) { //keep a count of c/cbar in event
        if (particle1->pdg_id() == 4) c_count++;
        if (particle1->pdg_id() == -4) cbar_count++;
    }

    if (c_count > 0 && cbar_count > 0) {//begin analysis only for events with at least 1 c/cbar pair
        D0_notsignal = false;
        D0bar_notsignal = false;
        D0 = false;
        D0bar = false;

        all_vectors.clear();
        d0bar_p_phi_eta.clear();
        d0_p_phi_eta.clear();
        for (const auto &vertex: event.vertices()) {//for all vertices in the event
            bool flavor_exc = false;
            //if(vertex->id()!=-108) continue;


            //define 2 vectors - one for the incoming and one for the outgoing particles of each event
            const std::vector<std::shared_ptr<const HepMC3::GenParticle>>& particlesOut = vertex->particles_out();
            const std::vector<std::shared_ptr<const HepMC3::GenParticle>>& particlesIn = vertex->particles_in();

            int count4 = 0;
            int count_4 = 0;

            //Loop over all outgoing particles and count the number of outgoing c and cbar
            for (const auto &particle: particlesOut) {
                if (particle->pdg_id() == 4) count4++;
                if (particle->pdg_id() == -4) count_4++;
            }

            //begin analysis for vertices with at least 1 outgoing c/cbar pairs
            if (count4 > 0 && count_4 > 0) {

                //variables initialization
                D0_notsignal=false;
                D0bar_notsignal=false;
                D0 = false;
                D0bar = false;
                newd0.clear();
                newd0bar.clear();

                for (const auto &particle3: particlesOut) {//loop over the outgoing c and cbar
                    //send the c and cbar particles to the trace daughters function.  see function for details
                    //if outgoing particle is not c/cbar skip
                    if (particle3->pdg_id() == 4) tracedaughters(particle3, D0, D0bar, newd0, newd0bar, d0_p_phi_eta, d0bar_p_phi_eta, D0_notsignal, D0bar_notsignal);
                    if (particle3->pdg_id() == -4) tracedaughters(particle3, D0, D0bar, newd0, newd0bar, d0_p_phi_eta, d0bar_p_phi_eta, D0_notsignal, D0bar_notsignal);
                    else continue;
                    cout << "D0=" << D0 << endl;
                }
                cout << "new vertex" << endl;

                if(D0 && !D0bar) {
                    cout << "D0 only, break for checking. ifile=" << ifile << " event =" << event_no << "vertex" << vertex_no << endl;
                    break;
                }
                if(!D0 && D0bar) {
                    cout << "Dbar only, break for checking. ifile=" << ifile << " event =" << event_no << "vertex" << vertex_no << endl;
                    break;
                }

                if (D0 and D0bar) {//if the vertex leads to a D0 and D0bar that both decay into K/Pi pairs
                    //signalvalues.clear();
                    cout << "*^*^*^*^*^*^*" << endl;
                    cout << "signal!!!" << endl;


                    signalvalues[0] = newd0[0]; //k- index
                    signalvalues[1] = newd0[1]; //pi- index
                    signalvalues[2] = newd0[2]; //d0 index
                    signalvalues[3] = newd0bar[0]; //k+ index
                    signalvalues[4] = newd0bar[1]; //pi- index
                    signalvalues[5] = newd0bar[2]; //d0bar index
                    signalvalues[6] = vertex->id();//vertex id


                    //begin analysis to determine the physics process of creating c/cbar
                    for (const auto &particle3: particlesOut) {
                        if (particle3->status() < 60 && particle3->status() > 50) {//if the c/cbar have pythia code inidacting created during final state showers
                            phys_process = 100; //set an arbitrary physics process to detect errors
                            for (const auto &particle4: particlesIn) {
                                //if the incoming particles of the c/cbar vertex contain 1 gluon, declare gluon splitting
                                //if (particlesIn.size() == 1 && particle4->pdg_id() == gluon_pid) phys_process = gluon_splitting;//gluon splitting
                                if (particlesIn.size() == 1 && particle4->pdg_id() == 21) phys_process = 1;//gluon splitting
                            }
                        }
                        else if (particle3->status() <50) {//only for particles from initial state showers, gluon fusion and flavor excitation
                            phys_process = 100;//to note errors in initial state showers
                            for (const auto &particle4: particlesIn) {//if both d0 and d0bar, begin analysis of incoming particles of vertex
                            //cout << "particles in size =" << particlesIn.size() << endl;
                            //cout << "particle pdg =" << particle4->pdg_id() << endl;

                                //if there are 2 incoming particles and they are gluons, gluon fusion, they are quarks, quark annihlation
                                //if (particlesIn.size() == 2 && particle4->pdg_id() == 21) phys_process = gluon_fusion;//gluon fusion
                                if (particlesIn.size() == 2 && particle4->pdg_id() == 21) phys_process = 2;//gluon fusion
                                //if (particlesIn.size() == 2 && particle4->pdg_id() < 7) phys_process = quark_annihlation; //quark annihlation
                                if (particlesIn.size() == 2 && particle4->pdg_id() < 7) phys_process = 5; //quark annihlation

                                //if one gluon comes in and the particles are initial state showers begin analysis for flavor excitation
                                if (particlesIn.size() == 1 && particle4->pdg_id() == 21) {//begin analysis to determine flavor excitation
                                    int gluon_in = 0,  c_in = 0, c_out = 0, cbar_in = 0, cbar_out = 0, q_in=0;

                                    for (const auto &vertex: event.vertices()) {//loop over all vertices to only select the ones with g+c comming in and g+c coming out
                                        for (const auto &inparticle: vertex->particles_in()) {//counts incoming gluons and c/cbar
                                            if (inparticle->pdg_id() == 21) gluon_in++;
                                            if (inparticle->pdg_id() == 4) c_in++;
                                            if (inparticle->pdg_id() == -4) cbar_in++;
                                            if (inparticle->pdg_id()< 7 && abs(inparticle->pdg_id())!=4) q_in++;
                                        }//counts incoming gluons and c/cbar
                                        for (const auto &outparticle: vertex->particles_out()) {//counts outgoing gluons and c/cbar
                                            if (outparticle->pdg_id() == 4) c_out++;
                                            if (outparticle->pdg_id() == -4) cbar_out++;
                                        }//counts outgoing gluons and c/cbar

                                        //cout << "gluon in =" << gluon_in << " q_in=" << q_in << " c_in=" << c_in << " c_out" << c_out << endl;

                                        //cout <<"particles in =" << endl;
                                        for (const auto &inparticle2: vertex->particles_in()) { //analyze vertices with 1 charm & 1 gluon or quark coming in and one charm going out
                                            //cout << "particle pdg = " << inparticle2->pdg_id() << endl;
                                            if (abs(inparticle2->pdg_id()) == 4 && gluon_in == 1 and c_in == 1 and c_out == 1) {
                                                int c_index = inparticle2->id();
                                                tracemothers(inparticle2, c_index, flavor_exc); //see trace mothers for more details
                                                if(flavor_exc) phys_process=3;
                                            }
                                            if (abs(inparticle2->pdg_id()) == 4 && q_in == 1 and c_in == 1 and c_out == 1) {
                                                int c_index = inparticle2->id();
                                                tracemothers(inparticle2, c_index, flavor_exc); //see trace mothers for more details
                                                if(flavor_exc) phys_process=4;
                                            }
                                        }
                                        if(!flavor_exc) phys_process=100;

                                        //cout << "The end " << endl;
                                        //flavor_exc = false;
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
                        else if(particle3->status() == 61 || particle3->status()==63) phys_process=100; 
                        else phys_process=100;
                    }
                    cout << "phys process = " << phys_process << endl;
                    signalvalues[7] =phys_process;//add physics process
                    all_vectors.push_back(signalvalues); //add signalvalues to all_vectors
                    for (int i = 0; i < 8; ++i)
                    {
                        signalvalues[i] = 0;
                    }
                }//if D0 and D0bar ** add analysis here
                /*
                if(D0 && D0bar_notsignal) {
                    for (int i = 0; i < 7; ++i)
                    {
                        signalvalues[i] = 0;
                    }
                    cout << "Hey!11111111111111111111" << endl;
                    signalvalues[7]=101;
                    all_vectors.push_back(signalvalues);
                    signalvalues[7] = 0;
                }
                if(D0bar && D0_notsignal) {
                    for (int i = 0; i < 7; ++i)
                    {
                        signalvalues[i] = 0;
                    }
                    cout << "Hey!111111111111111111" << endl;
                    signalvalues[7]=101;
                    all_vectors.push_back(signalvalues);
                    signalvalues[7] = 0;
                }
                */

            }//only for vertices containing at least 1 c/cbar pair
        }//for all vertices
    }//only for events with same number of c/cbar
    if (all_vectors.size() == 0)
    {
        for (int i = 0; i < 8; ++i)
        {
            signalvalues[i] = 0;
        }
        all_vectors.push_back(signalvalues);
    }
    if (d0_p_phi_eta.size() == 0) fill(d0_p_phi_eta.begin(), d0_p_phi_eta.end(), 0);
    if (d0bar_p_phi_eta.size() == 0) fill(d0bar_p_phi_eta.begin(), d0bar_p_phi_eta.end(), 0);
    return;
}//trace for signal

void tracedaughters(const std::shared_ptr<const HepMC3::GenParticle>& particle6, bool& d0, bool& d0bar,  std::vector<int> &newd0, std::vector<int> &newd0bar,  std::vector<float> &d0_p_phi_eta,  std::vector<float> &d0bar_p_phi_eta, bool& d0_notsignal, bool& d0bar_notsignal){
    //this is a recursive function that traces the daughters of each c/cbar
    const std::vector<std::shared_ptr<const HepMC3::GenParticle>> &daughters = particle6->children();//create a vector of all the daughers of the particle

   // cout << "made it this far" << endl;
    for (const auto &daughter: daughters) {//analyze each daughter
        //cout << "daughter pid = " << daughter->pdg_id() << " and status!" << daughter->status() << endl; 
        if (abs(daughter->pdg_id()) ==d0_pid) {//if the daughter is d0/d0bar
            
            //cout << "daught pid =" << daughter->pdg_id() << endl; 
            const auto& daughters = daughter->children();//anaylze the d0 daughters
            if (daughters.size() == 2 && (abs(daughters[0]->pdg_id()) == kaon_pid && abs(daughters[1]->pdg_id()) == pion_pid)) {//if k then pi
                if (daughter->pdg_id()==d0_pid && !d0) {
                    d0 = true;
                    //cout << "do! v1" << endl;
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
                if (daughter->pdg_id() == -d0_pid && !d0bar) {//same for d0bar
                    d0bar = true;
                    //cout << "dob! v1" << endl;
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
                if (daughter->pdg_id()==d0_pid && !d0 ) {
                    d0 = true;
                    //cout << "do! v2" << endl;
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
                if (daughter->pdg_id() == -d0_pid && !d0bar) {
                    d0bar = true;
                    //cout << "dob! v2" << endl;
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
        else if (abs(daughter->pdg_id())==4 || daughter->status()==74|| (abs(daughter->status())<100 && abs(daughter->status())>80 && daughter->pdg_id()!=21) || (daughter->status()==2 && daughter->pdg_id()!=21)) tracedaughters(daughter, d0, d0bar, newd0, newd0bar, d0_p_phi_eta, d0bar_p_phi_eta, d0_notsignal, d0bar_notsignal);
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

float transition_phi(float& D0del_phi){
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

void check_for_duplicates(std::vector<float> vertex_no,
                          std::vector<float> d0_index,
                          std::vector<float> d0bar_index,
                          std::vector<float>& phy_process)
{
    for (int j = 0; j < vertex_no.size() - 1; j++)
    {
        for (int k = j + 1; k < vertex_no.size(); k++)
        {
            // Check if d0_index and d0bar_index are the same for both vertex numbers
            if ((vertex_no)[j] != (vertex_no)[k])
            {
                if ((d0_index)[j] == (d0_index)[k] || (d0bar_index)[j] == (d0bar_index)[k])
                {
                    /*
                    cout << "VERTEX_NO SIZE=" << vertex_no.size() << endl;
                    cout << "phys processj =" << (phy_process)[j] << endl;
                    cout << "d0_index=" << (d0_index)[j] << " and " << (d0_index)[k] << endl;
                    cout << "d0bar_index=" << (d0bar_index)[j] << " and " << (d0bar_index)[k] << endl;
                    cout << "new phys processj =" << (phy_process)[j] << endl;
                    cout << endl;
                    */
                    (phy_process)[j] = 100;
                    (phy_process)[k] = 100;
                    
                }
            }
        }
    return;
    }
}
