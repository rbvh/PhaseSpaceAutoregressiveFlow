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
    cout << "Error: Provide an event file." << endl;
    exit(0);
  }

  ttbar_base ttbar;

  fstream fin;
  fin.open(argv[1], ios::in);

  string file_string = string(argv[1]);
  string file_string_out = file_string.substr(0, file_string.find(".csv", 0)) + "_phase_space.csv";

  ofstream fout(file_string_out);

  vector<double> sample;

  int count = 0;
  string event_line, event_num;
  while(fin >> event_line) {
    sample.clear();

    stringstream ss(event_line);
    while(getline(ss, event_num, ',')) {
      sample.push_back(stod(event_num));
    }

    vector<Vec4> sample_phase_space = ttbar.convert_to_phase_space(sample);

    if (sample_phase_space.empty()) continue;

    for (int i=0; i<sample_phase_space.size(); i++) {
      if (i == sample_phase_space.size()-1) {
        fout << scientific << sample_phase_space[i].E << "," << sample_phase_space[i].px << "," << sample_phase_space[i].py << "," << sample_phase_space[i].pz;  
      }
      else {
        fout << scientific << sample_phase_space[i].E << "," << sample_phase_space[i].px << "," << sample_phase_space[i].py << "," << sample_phase_space[i].pz << ",";
      }
    }
    fout << endl;

    count++;
  }

}