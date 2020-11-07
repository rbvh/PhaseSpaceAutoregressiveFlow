#include <iostream>
#include <iomanip> 
#include <random>
#include <cmath>
#include <string>
#include <fstream>
#include <algorithm>

#include "common_elements.hh"

int main(int argc, char *argv[]) {
  if (argc != 2) {
    cout << "Error: Provide a weights file." << endl;
    exit(0);
  }

  fstream fin_weights;
  fin_weights.open(argv[1], ios::in);

  vector<double> weights;

  string weight_num;

  fin_weights >> weight_num;

  while(fin_weights >> weight_num) {
    weights.push_back(stod(weight_num));
  }

  compute_metrics(weights);
}