#include "Pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include <iostream>
#include <fstream>

using namespace Pythia8;
using namespace std;
using namespace fastjet;

int main() {
  Pythia pythia;
  
  pythia.readString("Beams:frameType = 4");
  pythia.readString("Beams:LHEF = ttbar/Events/run_01_decayed_1/events.lhe.gz");
  
  // Import the MC@NLO command file 
  pythia.readFile("main-nlo.cmd");

  // Some extra settings
  pythia.readString("PartonLevel:Remnants          = off");
  pythia.readString("HadronLevel:All               = off");
  pythia.readString("TimeShower:nGluonToQuark      = 4");
  pythia.readString("SpaceShower:nQuarkIn          = 4");
  pythia.readString("Check:Event                   = off");
  pythia.readString("ProcessLevel:resonanceDecays  = off");

  pythia.init();

  Event &event = pythia.event;

  ofstream outfSamples("negative_weight_samples.csv");
  ofstream outfWeights("negative_weight_weights.csv");

  // Fastjet analysis - select algorithm and parameters
  double Rparam = 0.4;
  fastjet::Strategy               strategy = fastjet::Best;
  fastjet::RecombinationScheme    recombScheme = fastjet::E_scheme;
  fastjet::JetDefinition         *jetDef = NULL;
  jetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm, Rparam,
                                      recombScheme, strategy);

  // Fastjet input
  vector <fastjet::PseudoJet> fjInputs;

  int nPos = 0, nNeg = 0;

  double Emax = sqrt(pythia.info.s())/2;
  double Emin = 0.1;
  int nSkippedPt = 0;
  int nSkippedMass = 0;

  double mt = pythia.particleData.m0(6);
  double mW = pythia.particleData.m0(24);

  for (int iEvent = 0; ; ++iEvent) {
    if (pythia.info.atEndOfFile()) break;

    int ie = -1, imu = -1;
    int inue = -1, inumu = -1;
    if (pythia.next()) {       
      // Reset Fastjet input
      fjInputs.resize(0);

      for (int i=event.size()-1; i>=0; i--) {
        if (!event[i].isFinal()) continue;

        if      (abs(event[i].id()) == 11) ie = i;
        else if (abs(event[i].id()) == 12) inue = i;
        else if (abs(event[i].id()) == 13) imu = i;
        else if (abs(event[i].id()) == 14) inumu = i;
        else {
          // It's a jet constituent
          fjInputs.push_back( fastjet::PseudoJet( pythia.event[i].px(),
            pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e() ) );

          // Tag the b quarks
          if (event[i].id() == 5) {
            fjInputs.back().set_user_index(5);
          }
          else if (event[i].id() == -5) {
            fjInputs.back().set_user_index(-5);
          }
        }
      }

      if (ie == -1 || inue == -1 || imu == -1 || inumu == -1) continue;

      // Run Fastjet algorithm
      if (fjInputs.size() == 0) {continue;}

      fastjet::ClusterSequence clustSeq(fjInputs, *jetDef);
      vector<PseudoJet> sortedJets = sorted_by_pt(clustSeq.inclusive_jets(20.0));

      // Look for the hardest b jet
      int ibJet = -1; 
      for (unsigned int i=0; i<sortedJets.size(); i++) {
        vector<PseudoJet> constituents = sortedJets[i].constituents();
        for (unsigned int j=0; j<constituents.size(); j++) {
          if (constituents[j].user_index() == 5) {
            ibJet = i;
            break;
          }
        }
        if (ibJet != -1) break;
      }      
      
      // Look for the hardest bbar jet
      int ibbarJet = -1;
      for (unsigned int i=0; i<sortedJets.size(); i++) {
        vector<PseudoJet> constituents = sortedJets[i].constituents();
        for (unsigned int j=0; j<constituents.size(); j++) {
          if (constituents[j].user_index() == -5) {
            ibbarJet = i;
            break;
          }
        }
        if (ibbarJet != -1) break;
      } 
      
      if (ibJet == -1 || ibbarJet == -1) continue;

      if (ibJet == ibbarJet) continue;

      // Find the W and top momenta
      Vec4 pWpls = event[imu].p() + event[inumu].p();
      Vec4 pWmin = event[ie].p()  + event[inue].p();
      Vec4 pt(sortedJets[ibJet].px() + pWpls.px(), 
              sortedJets[ibJet].py() + pWpls.py(), 
              sortedJets[ibJet].pz() + pWpls.pz(), 
              sortedJets[ibJet].E() + pWpls.e());
      Vec4 ptbar(sortedJets[ibbarJet].px() + pWmin.px(), 
              sortedJets[ibbarJet].py() + pWmin.py(), 
              sortedJets[ibbarJet].pz() + pWmin.pz(), 
              sortedJets[ibbarJet].E() + pWmin.e());

      // Skip events with too small top or W pts
      if (pWpls.pT() < Emin || pWmin.pT() < Emin || pt.pT() < Emin || ptbar.pT() < Emin) {
        nSkippedPt++;
        continue;
      }

      // Skip events with ridiculous resonance masses
      if (pt.mCalc() > 2*mt || ptbar.mCalc() > 2*mt || pWpls.mCalc() > 2*mW || pWmin.mCalc() > 2*mW){
        nSkippedMass++;
        continue;
      }

      vector<double> xOut; 

      // Top
      xOut.push_back(pt.mCalc()/mt/2);
      xOut.push_back(log(pt.pT()/Emin)/log(Emax/Emin));
      xOut.push_back(pt.theta()/M_PI);
      xOut.push_back(pt.phi()/2/M_PI + 0.5);
      // antiTop
      xOut.push_back(ptbar.mCalc()/mt/2);
      xOut.push_back(log(ptbar.pT()/Emin)/log(Emax/Emin));
      xOut.push_back(ptbar.theta()/M_PI);
      xOut.push_back(ptbar.phi()/2/M_PI + 0.5);
      // W+
      xOut.push_back(pWpls.mCalc()/mW/2);
      xOut.push_back(log(pWpls.pT()/Emin)/log(Emax/Emin));
      xOut.push_back(pWpls.theta()/M_PI);
      xOut.push_back(pWpls.phi()/2/M_PI + 0.5);
      // W-
      xOut.push_back(pWmin.mCalc()/mW/2);
      xOut.push_back(log(pWmin.pT()/Emin)/log(Emax/Emin));
      xOut.push_back(pWmin.theta()/M_PI);
      xOut.push_back(pWmin.phi()/2/M_PI + 0.5);
      // W+ decay - use the muon as reference
      xOut.push_back(event[imu].theta()/M_PI);
      xOut.push_back(event[imu].phi()/2/M_PI + 0.5);
      // W- decay - use the electron as reference
      xOut.push_back(event[ie].theta()/M_PI);
      xOut.push_back(event[ie].phi()/2/M_PI + 0.5);

      for (unsigned int i=0; i<xOut.size(); i++) {
        if (xOut[i] < 0 || xOut[i] > 1) {
          cout << "Warning: x outside (0,1) in feature " << i << endl;
          continue;
        }
      }
      
      // Momentum reconstruction for validation
      /*
      // Top momentum
      double top_pT_Reco = exp(xOut[1]*log(Emax/Emin))*Emin;
      double top_eta_Reco = -log(tan(xOut[2]*M_PI/2));
      double top_phi_Reco = (xOut[3] + 0.5)*M_PI*2;
      double top_m_Reco = xOut[0]*mt*2;
      Vec4 ptReco(top_pT_Reco*cos(top_phi_Reco), 
                  top_pT_Reco*sin(top_phi_Reco), 
                  top_pT_Reco*sinh(top_eta_Reco),
                  sqrt(pow2(top_pT_Reco*cosh(top_eta_Reco)) + pow2(top_m_Reco)));

      double antitop_pT_Reco = exp(xOut[5]*log(Emax/Emin))*Emin;
      double antitop_eta_Reco = -log(tan(xOut[6]*M_PI/2));
      double antitop_phi_Reco = (xOut[7] + 0.5)*M_PI*2;
      double antitop_m_Reco = xOut[4]*mt*2;
      Vec4 ptbarReco( antitop_pT_Reco*cos(antitop_phi_Reco), 
                      antitop_pT_Reco*sin(antitop_phi_Reco), 
                      antitop_pT_Reco*sinh(antitop_eta_Reco),
                      sqrt(pow2(antitop_pT_Reco*cosh(antitop_eta_Reco)) + pow2(antitop_m_Reco)));

      // W momentum
      double Wpls_pT_Reco = exp(xOut[9]*log(Emax/Emin))*Emin;
      double Wpls_eta_Reco = -log(tan(xOut[10]*M_PI/2));
      double Wpls_phi_Reco = (xOut[11] + 0.5)*M_PI*2;
      double Wpls_m_Reco = xOut[8]*mW*2;
      Vec4 pWplsReco( Wpls_pT_Reco*cos(Wpls_phi_Reco), 
                      Wpls_pT_Reco*sin(Wpls_phi_Reco), 
                      Wpls_pT_Reco*sinh(Wpls_eta_Reco),
                      sqrt(pow2(Wpls_pT_Reco*cosh(Wpls_eta_Reco)) + pow2(Wpls_m_Reco)));

      double Wmin_pT_Reco = exp(xOut[13]*log(Emax/Emin))*Emin;
      double Wmin_eta_Reco = -log(tan(xOut[14]*M_PI/2));
      double Wmin_phi_Reco = (xOut[15] + 0.5)*M_PI*2;
      double Wmin_m_Reco = xOut[12]*mW*2;
      Vec4 pWminReco( Wmin_pT_Reco*cos(Wmin_phi_Reco), 
                      Wmin_pT_Reco*sin(Wmin_phi_Reco), 
                      Wmin_pT_Reco*sinh(Wmin_eta_Reco),
                      sqrt(pow2(Wmin_pT_Reco*cosh(Wmin_eta_Reco)) + pow2(Wmin_m_Reco)));

      // Wpls decay
      double mu_eta_Reco = -log(tan(xOut[16]*M_PI/2));
      double mu_phi_Reco = (xOut[17] + 0.5)*M_PI*2;
      double mu_pT_Reco = pWplsReco.m2Calc()/2/(
            cosh(mu_eta_Reco)*pWplsReco.e() - 
            cos(mu_phi_Reco)*pWplsReco.px() - 
            sin(mu_phi_Reco)*pWplsReco.py() - 
            sinh(mu_eta_Reco)*pWplsReco.pz());
      Vec4 pmuReco(mu_pT_Reco*cos(mu_phi_Reco),
                  mu_pT_Reco*sin(mu_phi_Reco),
                  mu_pT_Reco*sinh(mu_eta_Reco),
                  mu_pT_Reco*cosh(mu_eta_Reco));
      Vec4 pnumuReco = pWplsReco - pmuReco;


      // Wmin decay
      double e_eta_Reco = -log(tan(xOut[18]*M_PI/2));
      double e_phi_Reco = (xOut[19] + 0.5)*M_PI*2;
      double e_pT_Reco = pWminReco.m2Calc()/2/(
            cosh(e_eta_Reco)*pWminReco.e() -
            cos(e_phi_Reco)*pWminReco.px() - 
            sin(e_phi_Reco)*pWminReco.py() - 
            sinh(e_eta_Reco)*pWminReco.pz()
      );

      Vec4 peReco(e_pT_Reco*cos(e_phi_Reco),
                  e_pT_Reco*sin(e_phi_Reco),
                  e_pT_Reco*sinh(e_eta_Reco),
                  e_pT_Reco*cosh(e_eta_Reco));
      Vec4 pnueReco = pWminReco - peReco;

      cout << "pt " << pt << ptReco << endl;
      cout << "ptbar " << ptbar << ptbarReco << endl;
      cout << "pWpls " << pWpls << pWplsReco << endl;
      cout << "pWmin " << pWmin << pWminReco << endl;
      cout << "pmu " << event[imu].p() << pmuReco << endl;
      cout << "pnumu " << event[inumu].p() << pnumuReco << endl;
      cout << "pe " << event[ie].p() << peReco << endl;
      cout << "pnue " << event[inue].p() << pnueReco << endl;
      */

      for (unsigned int i=0; i<xOut.size(); i++) {
        if (i == xOut.size() - 1) {outfSamples << setprecision(8) << xOut[i] << endl;}
        else                      {outfSamples << setprecision(8) << xOut[i] << ",";}
      }

      if (pythia.info.weights_detailed_vector[0] > 0) {nPos++;}
      else {nNeg++;}
      
      outfWeights << setprecision(8) << pythia.info.weights_detailed_vector[0] << endl;
    }
  }
}