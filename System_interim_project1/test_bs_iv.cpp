// test_bs_iv.cpp
// Small unit test for Black-Scholes price + implied-vol bisection.
// Works on macOS (clang) and Linux (g++). Does NOT use bits/stdc++.h.

#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;

double norm_cdf(double x) {
    return 0.5 * (1.0 + erf(x / sqrt(2.0)));
}

double bs_call_price(double S, double K, double r, double sigma, double tau) {
    if (tau <= 0.0) return max(0.0, S - K);
    if (S <= 0.0) return 0.0;
    double s = sqrt(tau);
    double d1 = (log(S / K) + (r + 0.5 * sigma * sigma) * tau) / (sigma * s);
    double d2 = d1 - sigma * s;
    return S * norm_cdf(d1) - K * exp(-r * tau) * norm_cdf(d2);
}

double implied_vol_bisect(double market_price, double S, double K, double r, double tau,
                          double tol = 1e-8, int maxIter = 200) {
    if (!(market_price == market_price)) return NAN;
    if (tau <= 0.0) return 0.0;
    double intrinsic = max(0.0, S - K * exp(-r * tau));
    if (market_price < intrinsic - 1e-12) return NAN;
    double lo = 1e-8, hi = 5.0;
    double flo = bs_call_price(S, K, r, lo, tau) - market_price;
    double fhi = bs_call_price(S, K, r, hi, tau) - market_price;
    // expand hi if needed
    for (int it = 0; it < 10 && flo * fhi > 0.0; ++it) { hi *= 2.0; fhi = bs_call_price(S, K, r, hi, tau) - market_price; }
    if (flo * fhi > 0.0) return NAN;
    double mid = 0.0;
    for (int it = 0; it < maxIter; ++it) {
        mid = 0.5 * (lo + hi);
        double fmid = bs_call_price(S, K, r, mid, tau) - market_price;
        if (fabs(fmid) < tol) return mid;
        if (fmid * flo > 0.0) { lo = mid; flo = fmid; } else { hi = mid; fhi = fmid; }
    }
    return mid;
}

int main() {
    cout << fixed << setprecision(8);

    // Known test case
    double S = 100.0;
    double K = 100.0;
    double r = 0.01;
    double sigma_true = 0.25;
    double tau = 0.5;

    double market = bs_call_price(S, K, r, sigma_true, tau);
    double implied = implied_vol_bisect(market, S, K, r, tau);

    cout << "market_price = " << market << "\n";
    cout << "true_sigma  = " << sigma_true << "\n";
    cout << "implied_iv  = " << implied << "\n";

    if (std::isnan(implied)) {
        cerr << "Implied vol routine returned NaN (failure)\n";
        return 2;
    }

    double diff = fabs(implied - sigma_true);
    cout << "abs(implied - true) = " << diff << "\n";

    if (diff < 1e-5) {
        cout << "Unit test PASSED\n";
        return 0;
    } else {
        cout << "Unit test FAILED (tolerance 1e-5)\n";
        return 1;
    }
}
