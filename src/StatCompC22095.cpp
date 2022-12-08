#include <Rcpp.h>
#include <cmath> 
using namespace Rcpp;
using namespace std;

//' @title B-S Option Pricing Model
//' @description use B-S option pricing model to compute the price of European call option
//' @param S underlying asset price
//' @param K exercise price
//' @param r risk free interest rate
//' @param sigma volatility of assets(year)
//' @param time term of contract
//' @return the price of one European call option
//' @examples
//' \dontrun{
//' option_price_call_black_scholes(4.019, 3.3, 0.022685, 0.4199, 0.05479)
//' }
//' @export
//[[Rcpp::export]]
double option_price_call_black_scholes(const double& S,
                                       const double& K,
                                       const double& r,
                                       const double& sigma,
                                       const double& time) 
{
  double N(const double& z);
  double time_sqrt = sqrt(time);
  double d1 = (log(S/K)+r*time)/(sigma*time_sqrt)+0.5*sigma*time_sqrt;
  double d2 = d1-sigma*time_sqrt;
  return S*N(d1)-K*exp(-r*time)*N(d2);
}

double N(const double& z) {
  if (z >  6.0) { return 1.0; }; // this guards against overflow
  if (z < -6.0) { return 0.0; };
  double b1 =  0.31938153;
  double b2 = -0.356563782;
  double b3 =  1.781477937;
  double b4 = -1.821255978;
  double b5 =  1.330274429;
  double p  =  0.2316419;
  double c2 =  0.3989423;
  double a=fabs(z);
  double t = 1.0/(1.0+a*p);
  double b = c2*exp((-z)*(z/2.0));
  double n = ((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
  n = 1.0-b*n;
  if ( z < 0.0 ) n = 1.0 - n;
  return n;
}



//' @title B-S Option Pricing Model
//' @description use B-S option pricing model to compute the price of European put option
//' @param S underlying asset price
//' @param K exercise price
//' @param r risk free interest rate
//' @param sigma volatility of assets(year)
//' @param time term of contract
//' @return the price of one European put option
//' @examples
//' \dontrun{
//' option_price_put_black_scholes(4.019, 4.0, 0.022685, 0.1717, 0.05479)
//' }
//' @export
//[[Rcpp::export]]
double option_price_put_black_scholes(const double& S,
                                      const double& K,
                                      const double& r,
                                      const double& sigma,
                                      const double& time)
{
  double time_sqrt = sqrt(time);
  double d1 = (log(S/K)+r*time)/(sigma*time_sqrt)+0.5*sigma*time_sqrt;
  double d2 = d1-sigma*time_sqrt;
  return K*exp(-r*time)*N(-d2)-S*N(-d1);
}