// delta_hedge_sim.cpp
// Compile: g++ -std=c++17 -O2 -o delta_hedge_sim delta_hedge_sim.cpp
// Or:     clang++ -std=c++17 -O2 -o delta_hedge_sim delta_hedge_sim.cpp
// Run:    ./delta_hedge_sim

#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <string>

using namespace std;

// Portable PI constant
const double PI = acos(-1.0);

/* ---------------------------
   Math helpers (normal pdf/cdf)
   --------------------------- */
double norm_pdf(double x) {
    return exp(-0.5 * x * x) / sqrt(2.0 * PI);
}
double norm_cdf(double x) {
    // Phi(x) = 0.5 * (1 + erf(x / sqrt(2)))
    return 0.5 * (1.0 + erf(x / sqrt(2.0)));
}

/* ---------------------------
   Black-Scholes call price & delta
   --------------------------- */
double bs_call_price(double S, double K, double r, double sigma, double tau) {
    if (tau <= 0.0) {
        return max(0.0, S - K);
    }
    double sqt = sqrt(tau);
    double d1 = (log(S / K) + (r + 0.5 * sigma * sigma) * tau) / (sigma * sqt);
    double d2 = d1 - sigma * sqt;
    return S * norm_cdf(d1) - K * exp(-r * tau) * norm_cdf(d2);
}
double bs_call_delta(double S, double K, double r, double sigma, double tau) {
    if (tau <= 0.0) {
        return (S > K) ? 1.0 : 0.0;
    }
    double sqt = sqrt(tau);
    double d1 = (log(S / K) + (r + 0.5 * sigma * sigma) * tau) / (sigma * sqt);
    return norm_cdf(d1);
}

/* ---------------------------
   Simulation + hedging function
   --------------------------- */
struct SimParams {
    double S0 = 100.0;
    double T = 0.4;       // years
    double mu = 0.05;
    double sigma = 0.24;
    double r = 0.025;
    int N = 100;          // steps
    int n_paths = 1000;
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    SimParams P;

    double dt = P.T / P.N;
    double sqrt_dt = sqrt(dt);

    // RNG (deterministic seed for reproducibility during development)
    // change the seed number (12345) to random_device rd() if you want non-deterministic runs
    mt19937 rng(12345);
    normal_distribution<double> nd(0.0, 1.0);

    int save_paths = min(100, P.n_paths);
    vector<vector<double>> saved_paths(save_paths, vector<double>(P.N + 1));
    vector<double> hedge_errors;
    hedge_errors.reserve(P.n_paths);

    for (int path = 0; path < P.n_paths; ++path) {
        // simulate S series
        vector<double> S(P.N + 1);
        S[0] = P.S0;
        for (int i = 0; i < P.N; ++i) {
            double Z = nd(rng);
            double dS = P.mu * S[i] * dt + P.sigma * S[i] * sqrt_dt * Z;
            S[i + 1] = S[i] + dS;
            if (S[i+1] <= 0) S[i+1] = 1e-8;
        }

        if (path < save_paths) {
            for (int i = 0; i <= P.N; ++i) saved_paths[path][i] = S[i];
        }

        // compute BS prices and deltas
        vector<double> V(P.N + 1), delta(P.N + 1);
        const double K = 105.0;
        for (int i = 0; i <= P.N; ++i) {
            double tau = P.T - i * dt;
            V[i] = bs_call_price(S[i], K, P.r, P.sigma, max(0.0, tau));
            delta[i] = bs_call_delta(S[i], K, P.r, P.sigma, max(0.0, tau));
        }
        V[P.N] = max(0.0, S[P.N] - K);
        delta[P.N] = (S[P.N] > K) ? 1.0 : 0.0;

        double B_prev = V[0] - delta[0] * S[0];
        double HE = 0.0;
        for (int i = 1; i <= P.N; ++i) {
            double grow = exp(P.r * dt);
            double B_i = delta[i-1] * S[i] + B_prev * grow - delta[i] * S[i];
            HE = delta[i-1] * S[i] + B_prev * grow - V[i];
            B_prev = B_i;
        }
        hedge_errors.push_back(HE);
    }

    // save sample paths
    {
        ofstream fout("paths_first100.csv");
        fout << "time";
        for (int p = 0; p < save_paths; ++p) fout << ",path_" << p;
        fout << '\n';
        for (int i = 0; i <= P.N; ++i) {
            double t = i * dt;
            fout << fixed << setprecision(6) << t;
            for (int p = 0; p < save_paths; ++p) fout << ',' << fixed << setprecision(6) << saved_paths[p][i];
            fout << '\n';
        }
    }

    // save hedging errors
    {
        ofstream fout("hedge_errors.csv");
        fout << "path,HE_T\n";
        for (size_t i = 0; i < hedge_errors.size(); ++i) {
            fout << i << ',' << fixed << setprecision(10) << hedge_errors[i] << '\n';
        }
    }

    // summary stats
    double sum = 0.0;
    for (double x : hedge_errors) sum += x;
    double mean = sum / hedge_errors.size();
    double accum = 0.0;
    for (double x : hedge_errors) accum += (x - mean) * (x - mean);
    double stddev = sqrt(accum / hedge_errors.size());

    cout << "Wrote paths_first100.csv and hedge_errors.csv\n";
    cout << "Hedging error at T across " << hedge_errors.size() << " paths:\n";
    cout << "  mean = " << mean << "\n";
    cout << "  stddev = " << stddev << "\n";

    return 0;
}
