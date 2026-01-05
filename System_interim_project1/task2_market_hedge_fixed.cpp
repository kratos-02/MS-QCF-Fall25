// task2_market_hedge_fixed.cpp
// Compile:
//   g++ -std=c++17 -O2 -o task2_market_hedge_fixed task2_market_hedge_fixed.cpp
// or
//   clang++ -std=c++17 -O2 -o task2_market_hedge_fixed task2_market_hedge_fixed.cpp
//
// Usage examples:
//   ./task2_market_hedge_fixed
//   ./task2_market_hedge_fixed interest.csv sec_GOOG.csv op_GOOG.csv
//   ./task2_market_hedge_fixed interest.csv sec_GOOG.csv op_GOOG.csv 560 2011-07-29
//   ./task2_market_hedge_fixed interest.csv sec_GOOG.csv op_GOOG.csv 560 2011-07-29 2011-07-06 2011-07-29
//
// Notes:
//  - If you provide only the three filenames, the program auto-selects the first CALL contract.
//  - If you provide strike + exdate, it uses those. Optional start_date and end_date (YYYY-MM-DD) restrict the output range.
//  - The program will automatically start hedging at the first date that has a market option mid for the chosen contract.
//  - Forward-fill uses only the *previous* known mid (no future quotes used for earlier dates).

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

/* ----------------- tiny CSV helpers ----------------- */
static inline string trim(const string &s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}
static inline vector<string> split_csv_line(const string &line) {
    vector<string> out; string cur; bool in_quote=false;
    for (size_t i=0;i<line.size();++i){
        char c=line[i];
        if (c=='"'){ in_quote=!in_quote; continue; }
        if (!in_quote && c==','){ out.push_back(trim(cur)); cur.clear(); }
        else cur.push_back(c);
    }
    out.push_back(trim(cur));
    return out;
}
static inline double safe_stod(const string &s) {
    try { size_t idx; double v = stod(s, &idx); (void)idx; return v; }
    catch(...) { return numeric_limits<double>::quiet_NaN(); }
}

/* ----------------- Black-Scholes helpers ----------------- */
const double PI = acos(-1.0);
double norm_cdf(double x) { return 0.5*(1.0+erf(x/sqrt(2.0))); }

double bs_call_price(double S,double K,double r,double sigma,double tau){
    if (tau <= 0.0) return max(0.0, S - K);
    if (S <= 0.0) return 0.0;
    double s = sqrt(tau);
    double d1 = (log(S/K) + (r + 0.5*sigma*sigma)*tau) / (sigma * s);
    double d2 = d1 - sigma*s;
    return S * norm_cdf(d1) - K * exp(-r * tau) * norm_cdf(d2);
}
double bs_call_delta(double S,double K,double r,double sigma,double tau){
    if (tau <= 0.0) return (S > K) ? 1.0 : 0.0;
    if (S <= 0.0) return 0.0;
    double s = sqrt(tau);
    double d1 = (log(S/K) + (r + 0.5*sigma*sigma)*tau) / (sigma * s);
    return norm_cdf(d1);
}

/* ----------------- Implied vol (bisection) ----------------- */
double implied_vol_bisect(double market_price,double S,double K,double r,double tau,
                          double tol=1e-6,int maxIter=200){
    if (!(market_price == market_price)) return numeric_limits<double>::quiet_NaN();
    if (tau <= 0.0) return 0.0;
    double intrinsic = max(0.0, S - K * exp(-r * tau));
    if (market_price < intrinsic - 1e-12) return numeric_limits<double>::quiet_NaN();
    double lo = 1e-8, hi = 5.0;
    double flo = bs_call_price(S,K,r,lo,tau) - market_price;
    double fhi = bs_call_price(S,K,r,hi,tau) - market_price;
    for (int k=0; k<10 && flo * fhi > 0; ++k) { hi *= 2.0; fhi = bs_call_price(S,K,r,hi,tau) - market_price; }
    if (flo * fhi > 0) return numeric_limits<double>::quiet_NaN();
    for (int it=0; it<maxIter; ++it){
        double mid = 0.5 * (lo + hi);
        double fmid = bs_call_price(S,K,r,mid,tau) - market_price;
        if (fabs(fmid) < tol) return mid;
        if (fmid * flo > 0) { lo = mid; flo = fmid; } else { hi = mid; fhi = fmid; }
    }
    return 0.5 * (lo + hi);
}

/* ----------------- CSV readers (simple) ----------------- */
bool read_interest_csv(const string &fname, map<string,double> &out_rates, vector<string> &ordered_dates){
    ifstream fin(fname); if (!fin.is_open()) return false;
    string line; bool header=false;
    while (getline(fin,line)){
        if (!header){ header=true; continue; }
        if (trim(line).empty()) continue;
        auto cols = split_csv_line(line);
        if (cols.size() < 2) continue;
        string date = cols[0];
        double v = safe_stod(cols[1]);
        if (!(v == v)) continue;
        out_rates[date] = v / 100.0; // percent -> decimal
        ordered_dates.push_back(date);
    }
    fin.close(); return true;
}
bool read_stock_csv(const string &fname, map<string,double> &out_close, vector<string> &ordered_dates){
    ifstream fin(fname); if (!fin.is_open()) return false;
    string line; bool header=false;
    while (getline(fin,line)){
        if (!header){ header=true; continue; }
        if (trim(line).empty()) continue;
        auto cols = split_csv_line(line);
        if (cols.size() < 2) continue;
        string date = cols[0];
        double v = safe_stod(cols[1]);
        if (!(v == v)) continue;
        out_close[date] = v;
        ordered_dates.push_back(date);
    }
    fin.close(); return true;
}
struct OptionRow { string date, exdate; char cp_flag; double strike, bid, ask; };
bool read_options_csv(const string &fname, vector<OptionRow> &rows){
    ifstream fin(fname); if (!fin.is_open()) return false;
    string line; bool header=false;
    while (getline(fin,line)){
        if (!header){ header=true; continue; }
        if (trim(line).empty()) continue;
        auto cols = split_csv_line(line);
        if (cols.size() < 6) continue;
        OptionRow r;
        r.date = cols[0];
        r.exdate = cols[1];
        string cp = cols[2];
        r.cp_flag = cp.empty() ? ' ' : cp[0];
        r.strike = safe_stod(cols[3]);
        r.bid = safe_stod(cols[4]);
        r.ask = safe_stod(cols[5]);
        if (!(r.strike == r.strike)) continue;
        rows.push_back(r);
    }
    fin.close(); return true;
}

/* ----------------- helper to find index of date in ordered vector ----------------- */
int find_index(const vector<string> &vec, const string &date){
    for (size_t i=0;i<vec.size();++i) if (vec[i] == date) return (int)i;
    return -1;
}

/* ----------------- main ----------------- */
int main(int argc, char **argv){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string interest_file = "interest.csv";
    string stock_file    = "sec_GOOG.csv";
    string options_file  = "op_GOOG.csv";
    double desired_strike = numeric_limits<double>::quiet_NaN();
    string desired_exdate = "";
    string start_date = "";
    string end_date = "";
    bool start_date_provided = false;

    // Arg parsing:
    // Support:
    //  ./prog
    //  ./prog interest stock options
    //  ./prog interest stock options strike exdate
    //  ./prog interest stock options strike exdate start_date end_date
    if (argc >= 4) {
        interest_file = argv[1];
        stock_file = argv[2];
        options_file = argv[3];
    }
    if (argc >= 6) {
        desired_strike = safe_stod(argv[4]);
        desired_exdate = string(argv[5]);
    }
    if (argc >= 7) {
        start_date = string(argv[6]);
        start_date_provided = true;
    }
    if (argc >= 8) {
        end_date = string(argv[7]);
    }

    // Read files
    map<string,double> rate_by_date, stock_by_date;
    vector<string> rate_dates, stock_dates;
    vector<OptionRow> option_rows;
    if (!read_interest_csv(interest_file, rate_by_date, rate_dates)) { cerr<<"Error reading "<<interest_file<<"\n"; return 1; }
    if (!read_stock_csv(stock_file, stock_by_date, stock_dates)) { cerr<<"Error reading "<<stock_file<<"\n"; return 1; }
    if (!read_options_csv(options_file, option_rows)) { cerr<<"Error reading "<<options_file<<"\n"; return 1; }

    if (stock_dates.empty()) { cerr<<"No stock dates loaded.\n"; return 1; }

    // Choose option contract
    OptionRow chosen; bool chosen_ok = false;
    if (desired_strike == desired_strike && !desired_exdate.empty()) {
        for (const auto &r : option_rows) {
            if (r.cp_flag == 'C' && fabs(r.strike - desired_strike) < 1e-8 && r.exdate == desired_exdate) {
                chosen = r; chosen_ok = true; break;
            }
        }
        if (!chosen_ok) cerr<<"Warning: specified strike+exdate not found; auto-selecting first CALL.\n";
    }
    if (!chosen_ok) {
        for (const auto &r : option_rows) {
            if (r.cp_flag == 'C') { chosen = r; chosen_ok = true; break; }
        }
    }
    if (!chosen_ok) { cerr<<"No CALL option found in options CSV.\n"; return 1; }

    double K = chosen.strike;
    string exdate = chosen.exdate;
    cout << "Selected contract: strike="<<K<<", exdate="<<exdate<<"\n";

    // Build mid map for chosen contract (may be sparse)
    unordered_map<string,double> option_mid_by_date;
    for (const auto &r : option_rows) {
        if (r.cp_flag != 'C') continue;
        if (fabs(r.strike - K) > 1e-8) continue;
        if (r.exdate != exdate) continue;
        if (!(r.bid == r.bid) || !(r.ask == r.ask)) continue;
        double mid = 0.5 * (r.bid + r.ask);
        option_mid_by_date[r.date] = mid;
    }

    // Map exdate to stock index (ex_index)
    int ex_index = find_index(stock_dates, exdate);
    if (ex_index == -1) {
        cerr<<"exdate "<<exdate<<" not found in stock dates; attempting nearest mapping\n";
        int found = -1;
        for (size_t i=0;i<stock_dates.size();++i) if (stock_dates[i] > exdate) { found = (int)i; break; }
        if (found == -1) for (int i=(int)stock_dates.size()-1;i>=0; --i) if (stock_dates[i] < exdate) { found = i; break; }
        if (found == -1) { cerr<<"Cannot find suitable exdate in stock dates.\n"; return 1; }
        ex_index = found; exdate = stock_dates[ex_index];
        cerr<<"Using stock date "<<exdate<<" as effective exdate\n";
    }

    // If user provided end_date, map it to index else use ex_index
    int t_end_index = ex_index;
    if (!end_date.empty()) {
        int idx = find_index(stock_dates, end_date);
        if (idx == -1) cerr<<"Warning: end_date not in stock_dates; using exdate as end\n";
        else t_end_index = idx;
    }

    // Find first index (stock_dates) that has a real option mid for chosen contract
    int first_quote_index = -1;
    for (size_t i = 0; i < stock_dates.size(); ++i) {
        if (option_mid_by_date.find(stock_dates[i]) != option_mid_by_date.end()) { first_quote_index = (int)i; break; }
    }
    if (first_quote_index == -1) {
        cerr<<"No option quotes found for chosen contract in the entire file. Cannot hedge.\n";
        return 1;
    }

    // Determine start index t0_index: if start_date provided, map it; else default to first_quote_index
    int t0_index = 0;
    if (start_date_provided) {
        int idx = find_index(stock_dates, start_date);
        if (idx == -1) {
            cerr<<"start_date not in stock_dates; using first quote date "<<stock_dates[first_quote_index]<<"\n";
            t0_index = first_quote_index;
        } else {
            // Do not allow start before first_quote_index (we must not use future quotes)
            t0_index = max(idx, first_quote_index);
            if (idx < first_quote_index) {
                cerr<<"start_date earlier than first option quote; starting at first quote date "<<stock_dates[first_quote_index]<<"\n";
            }
        }
    } else {
        t0_index = first_quote_index;
    }

    if (t_end_index < t0_index) { cerr<<"end_date occurs before start_date after adjustments. Exiting.\n"; return 1; }

    // Setup hedging loop variables
    const double DAYS_PER_YEAR = 252.0;
    const double dt = 1.0 / DAYS_PER_YEAR;

    vector<string> out_dates;
    vector<double> out_S, out_V, out_sigma, out_delta, out_HE, out_PNL, out_PNL_with_HE, out_delta_prev, out_Bprev, out_Bi;

    double last_known_mid = numeric_limits<double>::quiet_NaN();
    double last_known_sigma = numeric_limits<double>::quiet_NaN();

    bool first_row = true;
    double V0 = numeric_limits<double>::quiet_NaN();
    double B_prev = 0.0;
    double delta_prev = 0.0;

    // Loop from t0_index to t_end_index inclusive
    for (int i = t0_index; i <= t_end_index; ++i) {
        string date = stock_dates[i];
        // must have stock price
        if (stock_by_date.find(date) == stock_by_date.end()) continue;
        double S = stock_by_date[date];

        // Determine marketV (option mid) for this date:
        double marketV = numeric_limits<double>::quiet_NaN();
        auto it = option_mid_by_date.find(date);
        if (it != option_mid_by_date.end()) {
            marketV = it->second;
            last_known_mid = marketV; // update last known
        } else {
            // If no market quote today: forward-fill previous known mid if available
            if (last_known_mid == last_known_mid) {
                marketV = last_known_mid; // carry forward last known mid
            } else {
                // No previous mid known (shouldn't happen because we started at first_quote_index),
                // but guard against it: skip date.
                continue;
            }
        }

        // get interest rate for date if available
        double r = 0.0;
        auto itr = rate_by_date.find(date);
        if (itr != rate_by_date.end()) r = itr->second;

        // compute tau (business days left to ex_index)
        int days_left = max(0, ex_index - i);
        double tau = (double)days_left / DAYS_PER_YEAR;
        if (tau < 0.0) tau = 0.0;

        // invert implied vol (bisection) with fallback to last known sigma
        double sigma = implied_vol_bisect(marketV, S, K, r, tau);
        if (!(sigma == sigma)) {
            if (last_known_sigma == last_known_sigma) sigma = last_known_sigma;
            else sigma = 0.24;
        }
        last_known_sigma = sigma;

        double delta = bs_call_delta(S, K, r, sigma, tau);

        // initialize hedge on first row (we start at first quote or user-selected date >= first quote)
        if (first_row) {
            V0 = marketV;
            B_prev = V0 - delta * S; // initial cash after selling option V0 and buying delta*S
            delta_prev = delta;
            double HE = delta_prev * S + B_prev * exp(r * dt) - marketV;
            out_dates.push_back(date);
            out_S.push_back(S);
            out_V.push_back(marketV);
            out_sigma.push_back(sigma);
            out_delta.push_back(delta);
            out_HE.push_back(HE);
            out_PNL.push_back(V0 - marketV);
            out_PNL_with_HE.push_back(HE);
            out_delta_prev.push_back(numeric_limits<double>::quiet_NaN());
            out_Bprev.push_back(B_prev);
            out_Bi.push_back(B_prev);
            first_row = false;
            continue;
        }

        // hedging update using previous delta and B_prev
        double grow = exp(r * dt);
        double B_i = delta_prev * S + B_prev * grow - delta * S; // cash after rebalancing
        double HE = delta_prev * S + B_prev * grow - marketV;    // hedging error (cum.)
        double PNL = V0 - marketV;
        double PNL_with_HE = HE;

        // save
        out_dates.push_back(date);
        out_S.push_back(S);
        out_V.push_back(marketV);
        out_sigma.push_back(sigma);
        out_delta.push_back(delta);
        out_HE.push_back(HE);
        out_PNL.push_back(PNL);
        out_PNL_with_HE.push_back(PNL_with_HE);
        out_delta_prev.push_back(delta_prev);
        out_Bprev.push_back(B_prev);
        out_Bi.push_back(B_i);

        // update prevs
        delta_prev = delta;
        B_prev = B_i;
    }

    // Write CSV
    ofstream fout("result.csv");
    fout << fixed << setprecision(8);
    fout << "date,S,V,implied_sigma,delta,HE,PNL,PNL_with_HE,delta_prev,B_prev,B_i\n";
    for (size_t i = 0; i < out_dates.size(); ++i) {
        fout << out_dates[i] << "," << out_S[i] << "," << out_V[i] << "," << out_sigma[i] << "," << out_delta[i] << ","
             << out_HE[i] << "," << out_PNL[i] << "," << out_PNL_with_HE[i] << ",";
        if (out_delta_prev[i] == out_delta_prev[i]) fout << out_delta_prev[i];
        fout << ",";
        fout << out_Bprev[i] << "," << out_Bi[i] << "\n";
    }
    fout.close();

    cout << "Wrote result.csv with " << out_dates.size() << " rows (columns: date,S,V,implied_sigma,delta,HE,PNL,PNL_with_HE,delta_prev,B_prev,B_i)\n";

    // Print summary of hedging error across rows
    if (!out_PNL_with_HE.empty()) {
        double sum = 0.0; for (double x: out_PNL_with_HE) sum += x;
        double mean = sum / out_PNL_with_HE.size();
        double accum = 0.0; for (double x: out_PNL_with_HE) accum += (x - mean) * (x - mean);
        double stddev = sqrt(accum / out_PNL_with_HE.size());
        cout << "Hedging error summary over output rows: mean=" << mean << ", std=" << stddev << "\n";
    } else {
        cout << "No rows output. Check that option has quotes in the chosen date range.\n";
    }

    return 0;
}
