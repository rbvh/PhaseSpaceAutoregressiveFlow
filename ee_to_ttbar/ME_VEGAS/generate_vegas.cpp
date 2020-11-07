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

// ------------------------------------------------------------------------
// --------------------------- ttbar_generator ----------------------------
// ------------------------------------------------------------------------
struct vegas_generator : public ttbar_base {
  vegas_generator(bool BW_prior=true): urd(0,1) {
    if (BW_prior) {
      // Use a mass-dependent BW for the tops and a regular one for the Ws
      // Get n_grid * n_per_bin points
      int n_per_bin = 1e5;
      vector<double> bw_points;

      // Sample lots of points for the tops
      cout << "Setting up Breit-Wigner prior for VEGAS grid" << endl;
      for (int i=0; i<n_grid*n_per_bin; i++) {
        double w_temp;
        bw_points.push_back(sample_bw_mass_dep(6, w_temp));
      }
      // Sort points
      sort(bw_points.begin(), bw_points.end());
      // Get bin edges
      for (int i=0; i<n_grid-1; i++) {
        grid[0][i] = bw_points[n_per_bin*i];
        grid[1][i] = bw_points[n_per_bin*i];
      }
      grid[0][n_grid-1] = 1.;
      grid[1][n_grid-1] = 1.;
      
      bw_points.clear();
      // Sample lots of points for the Ws
      for (int i=0; i<n_grid*n_per_bin; i++) {
        double w_temp;
        bw_points.push_back(sample_bw_mass(24, w_temp));
      }
      // Sort points
      sort(bw_points.begin(), bw_points.end());
      // Get bin edges
      for (int i=0; i<n_grid-1; i++) {
        grid[2][i] = bw_points[n_per_bin*i];
        grid[3][i] = bw_points[n_per_bin*i];
      }
      grid[2][n_grid-1] = 1.;
      grid[3][n_grid-1] = 1.;

      // Fill the rest uniformly
      double increment = 1/(double)n_grid;
      for (int i=4; i<grid.size(); i++) {
        double x_now = 0;
        for (int j=0; j<grid[i].size(); j++) {
          x_now += increment;
          grid[i][j] = x_now;
        }
      }
    }

    else {
      // Initialize grid uniformly
      double increment = 1/(double)n_grid;
      for (int i=0; i<grid.size(); i++) {
        double x_now = 0;
        for (int j=0; j<grid[i].size(); j++) {
          x_now += increment;
          grid[i][j] = x_now;
        }
      }
    }
  }

  double rng() {return urd(gen);}

  double ess() {
    double weights_sum = 0;
    double weights_2_sum = 0;
    for (int i=0; i<weights.size(); i++) {
      weights_sum += weights[i];
      weights_2_sum += weights[i]*weights[i];
    }
    return weights_sum*weights_sum/weights_2_sum;
  }

  double unweighting_efficiency() {
    double weight_max = 0;
    double weight_sum = 0;
    for(int i=0; i<weights.size(); i++) {
      if (weights[i] > weight_max) {weight_max = weights[i];}
      weight_sum += weights[i];
    }

    return weight_sum/weight_max/weights.size();
  }

  // ------------------------------------------------------------------------
  // ------------------------ Breit-Wigner samplers -------------------------
  // ------------------------------------------------------------------------
  double sample_bw_mass(int id, double &w) {
    // Sample from a Breit-Wigner
    double m = (id == 6) ? mt : mw;
    double g = (id == 6) ? gt : gw;
    double m2 = m*m;
    double g2 = g*g;

    double r = rng();
    
    double m2out = m2 + m*g*tan(r*atan((s - m2)/m/g) - (1.-r)*atan(m/g));
    double xOut = m2out/s;

    w *= pi* ( (m2out - m2)*(m2out - m2) + m2*g2 ) / m / g / s;

    return xOut;
  }

  double sample_bw_mass_dep(int id, double &w) {
    // Sample from a mass-dependent width Breit-Wigner
    // Follows algorithm of AcerMC
    // Returns x = m^2/s with an appropriate weight for that transform

    // Note that this function doesn't care about the phase space boundaries 
    // but this is effectively achieved by setting the matrix element to zero in those regions
    double m = (id == 6) ? mt : mw;
    double g = (id == 6) ? gt : gw;
    double m2 = m*m;
    double g2 = g*g;

    double eta_min = -m/g;
    double eta_max = (s - m2)/m/g;
    double del1 = (m/g)*(atan(eta_max) - atan(eta_min));
    double del2 = 0.5*log(1. + eta_max*eta_max) - 0.5*log(1. + eta_min*eta_min);
    double del = del1 + del2;

    double eta;
    while(true) {
      if (rng() < del2/del) {
        double x = del2*rng() + 0.5*log(1. + eta_min*eta_min);
        eta = sqrt( exp(2*x) - 1 );
      }
      else {
        double x = del1*rng() + (m/g)*atan(eta_min);
        eta = tan(x*g/m);
        
        double p1 = (m/g) / (1 + eta*eta);
        double p  = (m/g + eta) / (1 + eta*eta);
        if (rng() > p/p1) {eta = -eta;}
      }

      // Check if inside range
      if (eta > eta_min && eta < eta_max) {break;}
    }

    double m2out = m*g*eta + m2;
    double xOut = m2out/s;
    w *= del*( (m2out - m2)*(m2out - m2) + m2*g2 ) / m2out/s;

    return xOut;
  }

  // ------------------------------------------------------------------------
  // ------------------------------- Samplers -------------------------------
  // ------------------------------------------------------------------------
  void sample_from_grid(int m, bool store_events = true) {
    if (store_events) {
      samples.clear();
      weights.clear();
    }

    // Set func val sums to zero
    array<double, n_grid> zeros;
    zeros.fill(0.);
    est_sums.fill(zeros);

    int fraction_empty = 0;

    int count = 0;
    while(count < m) {
      vector<double> new_sample;
      vector<int> new_sample_bins;

      // Weight of the sampled point
      double weight = 1; 

      for (int d=0; d<dim_grid; d++) {
        // Select a bin
        int bin_now = rng()*n_grid;

        // Find bin edges
        double x_low = (bin_now == 0) ? 0. : grid[d][bin_now-1];
        double x_high = grid[d][bin_now];
        double delta = x_high - x_low;

        // Sample a point
        double x_now = x_low + delta*rng();

        // Update likelihood
        weight *= delta*n_grid;

        // Store point and bins
        new_sample.push_back(x_now);
        new_sample_bins.push_back(bin_now);
      }
      
      // Compute function value
      weight *= compute_cross_section_weight(new_sample);

      // Add square of weights to every appropriate bin
      for (int d=0; d<dim_grid; d++) {
        est_sums[d][new_sample_bins[d]] += weight*weight;
      }

      // Store
      if (store_events) {
        if (weight != 0.) {
          samples.push_back(new_sample);
          weights.push_back(weight);
          count++;
        }
      }
      else {
        count++;
      }
    }
  }

  // Recompute bin edge assignments - returns loss function
  double refine_grid() {
    double performance = 0;

    for (int d=0; d<dim_grid; d++) {
      // Stores the estimators of this dimension
      array<double, n_grid> grid_weights;

      // Square root of all estimators (they are a sum of squared weights)
      for (int i=0; i<n_grid; i++) {
        grid_weights[i] = sqrt(est_sums[d][i]);
      }
      
      // Smooth all weights
      // 0th bin
      double weight_prev = grid_weights[0];
      double weight_cur  = grid_weights[1];
      grid_weights[0] = (weight_prev + weight_cur)/2.;

      // 1st till n_grid-1-th bin
      for (int i=1; i<n_grid-1; i++) {
        double weight_temp = weight_prev + weight_cur;
        weight_prev = weight_cur;
        weight_cur = grid_weights[i+1];
        grid_weights[i] = (weight_temp + weight_cur)/3.;
      }

      // n_grid-th bin
      grid_weights[n_grid-1] = (weight_prev + weight_cur)/2.;

      // Promote small weights for stability
      // Find maximum weight 

      double weight_max = 0;
      for (int i=0; i<n_grid; i++) {
        if (grid_weights[i] > weight_max) {weight_max = grid_weights[i];}
      }
      // Adjust weights
      for (int i=0; i<n_grid; i++) {
        double w_now = grid_weights[i]/weight_max;
        grid_weights[i] = pow(w_now*(1 + c*(1. - w_now)), alpha)*weight_max;
      }

      // Now normalize weights
      double weight_sum = 0;
      for (int i=0; i<n_grid; i++) {
        weight_sum += grid_weights[i];
      }
      for (int i=0; i<n_grid; i++) {
        grid_weights[i] /= weight_sum;
      }
    
      // Recompute size of each bin
      array<double, n_grid> new_bins;
      double weight_average = 1./n_grid;
      double weight_accumulated = 0;
      int old_bin = -1;
      for (int i=0; i<n_grid-1; i++) {
        // Merge bins with weights below the average
        while(weight_accumulated < weight_average) {
          old_bin++;
          weight_accumulated += grid_weights[old_bin];
          //cout << "Add old_bin = " << grid_weights[old_bin] << " total = " << weight_accumulated << endl;
        }
        // Shrink the rightmost accumulated bin
        double shrink_fac = (weight_accumulated - weight_average)/grid_weights[old_bin];
        //cout << i << " " << old_bin << " " << shrink_fac << endl;
        double old_width = (old_bin == 0) ? grid[d][old_bin] : grid[d][old_bin] - grid[d][old_bin-1];
        new_bins[i] = grid[d][old_bin] - old_width*shrink_fac;

        weight_accumulated -= weight_average;
      }
      new_bins[n_grid-1] = 1.;
      grid[d] = new_bins;

      // Compute performance
      for (int i=0; i<n_grid; i++) {
        performance += fabs(grid_weights[i]-weight_average);
      }
    }
    return performance/dim_grid/n_grid;
  }

  void sample_unweighted_events(int n, bool refine_during_sampling = false) {
    // Make n unweighted samples
    
    // Refine grid after m samples
    int m = 1e6;
    // Burn in VEGAS n_burn times
    int n_burn = 10;

    cout << "Burning in VEGAS for " << n_burn << " cycles of " << m << " events" << endl;
    // Burn in the vegas grid for 5 cycles
    for (int i=0; i<n_burn; i++) {
      sample_from_grid(m, false);
      double loss_function = refine_grid();
      cout << "Burn-in cycle " << i << " VEGAS loss function " << loss_function << endl;
    }

    cout << "Generating " << n << " unweighted samples." << endl;
    vector<vector<double> > unweighted_samples;
    double weight_max = 0;

    int vegas_counter = 0;
    while(unweighted_samples.size() < n) {
      vegas_counter++;
      // Sample from grid and store 
      sample_from_grid(m, true);

      // Run over stored events
      for (int i=0; i<samples.size(); i++) {
        // Check if we found a new highest weight
        if (weights[i] > weight_max) {
          cout << "New highest weight " << weights[i] << endl;

          // Rejection sampling on the previous set
          double p_accept_max = weight_max/weights[i];

          vector<vector<double> > replace_samples;
          for (int i=0; i<unweighted_samples.size(); i++) {
            if (rng() < p_accept_max) {
              replace_samples.push_back(unweighted_samples[i]);
            }
          }
          unweighted_samples = replace_samples;
          weight_max = weights[i];
        }

        // Accept probability
        double p_accept = weights[i]/weight_max;
        if (rng() < p_accept) {
          unweighted_samples.push_back(samples[i]);
        }
      }

      // Refine the grid
      if (refine_during_sampling) {
        double loss_function = refine_grid();
        cout << "Sampled " << m*vegas_counter << " events, saved " << unweighted_samples.size();
        cout << " VEGAS loss function = " << loss_function << endl;
      }
      else {
        samples.clear();
        weights.clear();
        cout << "Sampled " << m*vegas_counter << " events, saved " << unweighted_samples.size() << endl;
      }
    }

    // Make sure it's exactly n events
    while(unweighted_samples.size() > n) {
      unweighted_samples.pop_back();
    }

    write_events(unweighted_samples, "unweighted");
  }

  void sample_weighted_events(int n) {
    // Make n weighted samples
    
    // Refine grid after m samples
    int m = 1e6;
    // Burn in VEGAS n_burn times (was originally 15)
    int n_burn = 10;

    cout << "Burning in VEGAS for " << n_burn << " cycles of " << m << " events" << endl;
    // Burn in the vegas grid for n_burn cycles
    for (int i=0; i<n_burn; i++) {
      sample_from_grid(m, false);
      double loss_function = refine_grid();
      cout << "Burn-in cycle " << i << " VEGAS loss function " << loss_function << endl;
    }

    // Now sample from the grid
    cout << "Sampling " << n << " events" << endl;
    sample_from_grid(n, true);

    double w_max = 0;
    double w_sum = 0;
    for (int i=0; i<weights.size(); i++) {
      if (weights[i] > w_max) w_max = weights[i];
      w_sum += weights[i];
    }
    cout << "Unweighting efficiency " << w_sum/w_max/weights.size() << endl;

    // And save the events
    write_events(samples, "weighted");
    write_weights(weights, "weighted");
  }

  void print_grid() {
    for (int i=0; i<grid.size(); i++) {
      for (int j=0; j<grid[i].size(); j++) {
        cout << grid[i][j] << " ";
      }
      cout << endl;
    }
  }

  // Set of samples (n_samples, dims)
  vector<vector<double> > samples;

  // Weights (n_samples)
  vector<double> weights;

  // Bins per dimension
  static const int n_grid = 500;
  // dim of function
  static const int dim_grid = 14;
  // Constant alpha
  constexpr static const double alpha = 1.;
  constexpr static const double c = 0.5;

  // Grid stores the right side of the edges (dims, n_grid)
  // So it doesn't store the 0, but does store the 1
  array<array<double, n_grid>, dim_grid> grid;

  // Function sums per bin (dims, n_grid)
  array<array<double, n_grid>, dim_grid> est_sums;

  mt19937 gen;
  uniform_real_distribution<double> urd;

  bool print = false;
};

int main() {
  // Initialize VEGAS with BW prior distribution
  vegas_generator vegas;

  int n_events = 1e6;

  // Sample unweighted events
  vegas.sample_unweighted_events(n_events, false);

  // Sample weighted events
  // vegas.sample_weighted_events(n_events);
}