#include <iostream>
#include <iomanip> 
#include <random>
#include <cmath>
#include <string>
#include <fstream>

#include "CPPProcess.h"

// ------------------------------------------------------------------------
// ----------------------- Importance sampling metrics --------------------
// ------------------------------------------------------------------------
void compute_metrics(vector<double> &w) {
  double w_average = 0;
  double w2 = 0;

  for (int i=0; i<w.size(); i++) {
    w_average += w[i]/w.size();
    w2 += w[i]*w[i];
  }

  cout << "Cross section " << w_average << " +- " << sqrt((w2 - w.size()*w_average*w_average))/w.size() << endl;
  cout << "ESS           " << w.size()*w_average*w_average/w2 << endl;

  // Compute the unweighting efficiencies and corresponding coverage
  sort(w.begin(), w.end());

  // Compute maximum quantile
  int max_quantile = 0;
  int size_w = w.size();
  while(size_w > 1) {
    size_w /= 10;
    max_quantile += 1;
  }

  for (int i = 2; i<=max_quantile; i++) {
    int quantile_index = w.size()*(1. - pow(10, -i));

    // Clipped average
    double w_average_clipped = 0;
    for (int i=0; i<w.size(); i++) {
      w_average_clipped += w[i] < w[quantile_index] ? w[i]/w.size() : w[quantile_index]/w.size();
    }
    cout << "UE (10^-" << i << ")    " << w_average_clipped/w[quantile_index] << " Coverage " << w_average_clipped/w_average << endl;
  }
}

// Kallen function
double sqrt_kallen(double x, double y, double z) {
  return sqrt(max(0., x*x + y*y + z*z - 2*x*y - 2*x*z - 2*y*z));
}

// ------------------------------------------------------------------------
// ------------------------------- Vec4 class -----------------------------
// ------------------------------------------------------------------------

struct Vec4 {
  Vec4(double E, double px, double py, double pz): E(E), px(px), py(py), pz(pz) {}

  double m2Calc() {
    return E*E - px*px - py*py - pz*pz;
  }
  double mCalc() {
    return sqrt(max(0., m2Calc()));
  }

  Vec4 operator+(Vec4& pOther) {
    return Vec4(E + pOther.E, px + pOther.px, py + pOther.py, pz + pOther.pz);
  }

  // Boost to rest frame of pRef
  void toRest(Vec4 &pRef) {
    double M = pRef.mCalc();
    double Eout = (E*pRef.E - px*pRef.px - py*pRef.py - pz*pRef.pz)/M;
    double pFac = (E + Eout)/(pRef.E + M);

    E = Eout;
    px -= pFac*pRef.px;
    py -= pFac*pRef.py;
    pz -= pFac*pRef.pz;
  }

  // Boost from rest frame of pRef
  void fromRest(Vec4 &pRef) {
    double M = pRef.mCalc();
    double Eout = (E*pRef.E + px*pRef.px + py*pRef.py + pz*pRef.pz)/M;
    double pFac = (E + Eout)/(pRef.E + M);

    E = Eout;
    px += pFac*pRef.px;
    py += pFac*pRef.py;
    pz += pFac*pRef.pz;
  }

  void print() {
    cout << "(" << E << ", " << px << ", " << py << ", " << pz << ") m2 = " << m2Calc() << endl;
  }
 
  double E, px, py, pz;
};

// ------------------------------------------------------------------------
// ---------------------- Matrix element evaluator class-------------------
// ------------------------------------------------------------------------

struct ttbar_base {

  ttbar_base() {
    // Hard-code the CoM energy for the moment
    s = 1E6;
    process.initProc("/Users/rob/Desktop/PhaseSpaceAutoregressiveFlow/ee_to_ttbar/ME_VEGAS/param_card.dat");

    mt = process.pars->mdl_MT;
    mb = process.pars->mdl_MB;
    mw = process.pars->mdl_MW;

    mt2 = mt*mt;
    mb2 = mb*mb;
    mw2 = mw*mw;

    mub2 = mb2/s;

    gt = process.pars->mdl_WW;
    gw = process.pars->mdl_WT;
  }

  // Convert unit box to phase space
  vector<Vec4> convert_to_phase_space(vector<double> &sample) {
    // This functions computes the cross section weight as a function of a hypercube point
    // As such, it also includes the relevant phase space/flux factors

    vector<Vec4> result;

    // Check the phase space 
    for (int i=0; i<(int)sample.size(); i++) {
      if (sample[i] >= 1 || sample[i] <= 0) {return result;}
    }

    // u-d momenta
    double cos_ud = sample[4]*2. - 1.;
    double sin_ud = sqrt(max(0., 1. - cos_ud*cos_ud));
    double phi_ud = sample[5]*2*pi;
    double E_ud = sqrt(sample[2]*s)/2;
    Vec4 pu(E_ud,  E_ud*sin_ud*sin(phi_ud),  E_ud*sin_ud*cos(phi_ud),  E_ud*cos_ud);
    Vec4 pd(E_ud, -E_ud*sin_ud*sin(phi_ud), -E_ud*sin_ud*cos(phi_ud), -E_ud*cos_ud);

    // e-v momenta
    double cos_ev = sample[6]*2. - 1.;
    double sin_ev = sqrt(max(0., 1. - cos_ev*cos_ev));
    double phi_ev = sample[7]*2*pi;
    double E_ev = sqrt(sample[3]*s)/2;
    Vec4 pe(E_ev,  E_ev*sin_ev*sin(phi_ev),  E_ev*sin_ev*cos(phi_ev),  E_ev*cos_ev);
    Vec4 pv(E_ev, -E_ev*sin_ev*sin(phi_ev), -E_ev*sin_ev*cos(phi_ev), -E_ev*cos_ev);

    // w+ b momenta
    double cos_wplsb = sample[8]*2. - 1.;
    double sin_wplsb = sqrt(max(0., 1. - cos_wplsb*cos_wplsb));
    double phi_wplsb = sample[9]*2*pi;
    double Ewpls = sqrt(s/sample[0])*(sample[0] + sample[2] - mub2)/2;
    double Ebpls = sqrt(s/sample[0])*(sample[0] - sample[2] + mub2)/2;
    double ppls =  sqrt(s/sample[0])*sqrt_kallen(sample[0], sample[2], mub2)/2;

    // Check phase space 
    if (Ewpls < 0 || Ebpls < 0 || ppls == 0.) {return result;}

    Vec4 pwpls(Ewpls,  ppls*sin_wplsb*sin(phi_wplsb),  ppls*sin_wplsb*cos(phi_wplsb),  ppls*cos_wplsb);
    Vec4 pbpls(Ebpls, -ppls*sin_wplsb*sin(phi_wplsb), -ppls*sin_wplsb*cos(phi_wplsb), -ppls*cos_wplsb);

    // w- b momenta
    double cos_wminb = sample[10]*2. - 1.;
    double sin_wminb = sqrt(max(0., 1. - cos_wminb*cos_wminb));
    double phi_wminb = sample[11]*2*pi;
    double Ewmin = sqrt(s/sample[1])*(sample[1] + sample[3] - mub2)/2;
    double Ebmin = sqrt(s/sample[1])*(sample[1] - sample[3] + mub2)/2;
    double pmin =  sqrt(s/sample[1])*sqrt_kallen(sample[1], sample[3], mub2)/2;

    // Check phase space 
    if (Ewmin < 0 || Ebmin < 0 || pmin == 0.) {return result;}
    
    Vec4 pwmin(Ewmin,  pmin*sin_wminb*sin(phi_wminb),  pmin*sin_wminb*cos(phi_wminb),  pmin*cos_wminb);
    Vec4 pbmin(Ebmin, -pmin*sin_wminb*sin(phi_wminb), -pmin*sin_wminb*cos(phi_wminb), -pmin*cos_wminb);

    // Top momenta
    double cos_tt = sample[12]*2. - 1.;
    double sin_tt = sqrt(max(0., 1. - cos_tt*cos_tt));
    double phi_tt = sample[13]*2*pi;
    double Et1 = sqrt(s)*(1 + sample[0] - sample[1])/2.;
    double Et2 = sqrt(s)*(1 - sample[0] + sample[1])/2.;
    double pt = sqrt(s)*sqrt_kallen(1, sample[0], sample[1])/2;

    // Check phase space 
    if (Et1 < 0 || Et2 < 0 || pt == 0.) {return result;}

    Vec4 pt1(Et1,  pt*sin_tt*sin(phi_tt),  pt*sin_tt*cos(phi_tt),  pt*cos_tt);
    Vec4 pt2(Et2, -pt*sin_tt*sin(phi_tt), -pt*sin_tt*cos(phi_tt), -pt*cos_tt);

    // Boost to w rest frames
    pu.fromRest(pwpls);
    pd.fromRest(pwpls);
    pe.fromRest(pwmin);
    pv.fromRest(pwmin);

    // Boost to t1 rest frames
    pu.fromRest(pt1);
    pd.fromRest(pt1);
    pbpls.fromRest(pt1);

    // Boost to t2 rest frames
    pe.fromRest(pt2);
    pv.fromRest(pt2);
    pbmin.fromRest(pt2);

    result.push_back(pbpls);
    result.push_back(pbmin);
    result.push_back(pe);
    result.push_back(pv);
    result.push_back(pu);
    result.push_back(pd);

    return result;    
  }

  // Generated with e+ e- > t t~ > b b~ e- ve~ u d~
  double compute_cross_section_weight(vector<double> &sample) {
    vector<Vec4> sample_phase_space = convert_to_phase_space(sample);
    
    if (sample_phase_space.empty()) {return 0;}

    vector<double*> pMG;
    for (int i=0; i<8; i++){
      pMG.push_back(new double[4]);
    }
    pMG[0][0] = sqrt(s)/2., pMG[0][1] = 0., pMG[0][2] = 0., pMG[0][3] = sqrt(s)/2.;
    pMG[1][0] = sqrt(s)/2., pMG[1][1] = 0., pMG[1][2] = 0., pMG[1][3] = -sqrt(s)/2.;

    pMG[2][0] = sample_phase_space[0].E, pMG[2][1] = sample_phase_space[0].px, pMG[2][2] = sample_phase_space[0].py, pMG[2][3] = sample_phase_space[0].pz;
    pMG[3][0] = sample_phase_space[1].E, pMG[3][1] = sample_phase_space[1].px, pMG[3][2] = sample_phase_space[1].py, pMG[3][3] = sample_phase_space[1].pz;
    pMG[4][0] = sample_phase_space[2].E, pMG[4][1] = sample_phase_space[2].px, pMG[4][2] = sample_phase_space[2].py, pMG[4][3] = sample_phase_space[2].pz;
    pMG[5][0] = sample_phase_space[3].E, pMG[5][1] = sample_phase_space[3].px, pMG[5][2] = sample_phase_space[3].py, pMG[5][3] = sample_phase_space[3].pz;
    pMG[6][0] = sample_phase_space[4].E, pMG[6][1] = sample_phase_space[4].px, pMG[6][2] = sample_phase_space[4].py, pMG[6][3] = sample_phase_space[4].pz;
    pMG[7][0] = sample_phase_space[5].E, pMG[7][1] = sample_phase_space[5].px, pMG[7][2] = sample_phase_space[5].py, pMG[7][3] = sample_phase_space[5].pz;
    
    // Set momenta for this event
    process.setMomenta(pMG);

    // Evaluate matrix element
    process.sigmaKin();

    double ME = process.getMatrixElements()[0];

    // Prevent memory leaks due to stupid MG5 architecture
    for (int i=0; i < pMG.size(); i++) {
      delete pMG[i];
    }

    // Weight consists of matrix element..
    double weight = ME;
    
    // .. some numerical factors
    weight *= pow(s,3)*pow(2, -20)*pow(pi, -9);
    // .. and some kallen functions
    weight *= sqrt_kallen(1., sample[0], sample[1]);
    weight *= sqrt_kallen(1., sample[2]/sample[0], mub2);
    weight *= sqrt_kallen(1., sample[3]/sample[1], mub2);

    return weight;
  }

  CPPProcess process;

  double mt, mb, mw;
  double mt2, mb2, mw2;
  double gt, gw;

  double mub2;
  double s;

  double pi = 3.14159265359;
};
