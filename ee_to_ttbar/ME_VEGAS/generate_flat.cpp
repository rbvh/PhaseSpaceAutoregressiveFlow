#include <iostream>
#include <iomanip> 
#include <random>
#include <cmath>
#include <string>
#include <fstream>
#include <algorithm>
#include <array>

#include "common_elements.hh"

// Write events to file
void write_events(vector<vector<double> > &samples, string file_name) {
  ofstream outf_sample;
  outf_sample.open("../data/" + file_name + "_samples.csv");
  for (int i=0; i<samples.size(); i++) {
    for (int j=0; j<samples[i].size(); j++) {
      if (j == samples[i].size() - 1) {
        outf_sample << scientific << setprecision(10) << samples[i][j] << endl;
      }
      else {
        outf_sample << scientific << setprecision(10) << samples[i][j] << ",";
      }
    }
  }
}

// Write weights to file
void write_weights(vector<double> &weights, string file_name) {
  ofstream outf_weight;
  outf_weight.open("../data/" + file_name + "_weights.csv");
  for (int i=0; i<weights.size(); i++) {
    outf_weight << scientific << setprecision(10) << weights[i] << endl;
  }
}

struct flat_generator : public ttbar_base {
  flat_generator(): urd(0,1) {}

  double rng() {return urd(gen);}

  void sample_flat_events(int n) {
    // Fill samples
    while(samples.size() < n) {
      vector<double> sample;
      for (int i=0; i<14; i++) {
        sample.push_back(rng());
      }

      double weight = compute_cross_section_weight(sample);

      if (weight > 0) {
        samples.push_back(sample);
        weights.push_back(weight);
      }
    }

    // And save the events
    write_events(samples, "flat");
    write_weights(weights, "flat");
  }



  // Set of samples (n_samples, dims)
  vector<vector<double> > samples;

  // Weights (n_samples)
  vector<double> weights;

  mt19937 gen;
  uniform_real_distribution<double> urd;

};

int main() {
  flat_generator flat;

  flat.sample_flat_events(1e7);
}
