#include <iostream>
#include <iomanip> 
#include <random>
#include <cmath>
#include <string>
#include <fstream>
#include <algorithm>

#include "common_elements.hh"

// Compute and print total cross section
void cross_section(vector<double> &w) {
  int n = w.size();
  double w_sum = 0;
  double w_sum2 = 0;
  for (int i=0; i<n; i++) {
    w_sum += w[i];
    w_sum2 += w[i]*w[i];
  }
  
  cout << "Cross section: " << w_sum/n << " += " << sqrt((w_sum2 - w_sum*w_sum/n))/n << " GeV^-2" << endl;
}

int main(int argc, char *argv[]) {
  if (argc != 3) {
    cout << "Error: Provide an event file and a likelihoods file." << endl;
    exit(0);
  }

  ttbar_base evaluator;

  fstream fin_events; 
  fstream fin_likelihoods;
  
  cout << argv[1] << endl;
  fin_events.open(argv[1], ios::in);
  fin_likelihoods.open(argv[2], ios::in);

  vector<double> sample;
  vector<double> weights;

  double likelihood;
  string event_line, event_num, likelihood_num;

  int count = 0;
  while(fin_events >> event_line) {
    fin_likelihoods >> likelihood_num;
    sample.clear();

    stringstream ss(event_line);
    while(getline(ss, event_num, ',')) {
      sample.push_back(stod(event_num));
    }
    double likelihood = stod(likelihood_num);

    double ME = evaluator.compute_cross_section_weight(sample);
    
    double weight = ME/likelihood;

    weights.push_back(weight);
    count++;
  }
  
  compute_metrics(weights);
}