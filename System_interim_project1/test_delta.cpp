// test_delta.cpp - verify Black-Scholes delta numeric
#include <cmath>
#include <iostream>
#include <iomanip>
using namespace std;
double norm_cdf(double x){ return 0.5*(1.0+erf(x/sqrt(2.0))); }
double bs_call_delta(double S,double K,double r,double sigma,double tau){
    if (tau<=0.0) return (S>K)?1.0:0.0;
    double s = sqrt(tau);
    double d1=(log(S/K)+(r+0.5*sigma*sigma)*tau)/(sigma*s);
    return norm_cdf(d1);
}
int main(){
    cout<<fixed<<setprecision(8);
    double S=120,K=100,r=0.01,sigma=0.2,tau=0.5;
    double delta=bs_call_delta(S,K,r,sigma,tau);
    cout<<"delta="<<delta<<"\n";
    if (delta>0.9 && delta<1.0) { cout<<"delta test PASSED\n"; return 0; }
    else { cout<<"delta test FAILED\n"; return 1; }
}
