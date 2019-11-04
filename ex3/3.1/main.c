/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <tuple>
#include <cmath>
#include "rng/random.h"

using namespace std;

tuple<vector<double>,vector<double>> blocking_method(vector<double> v, int n_step, int n_blocks);
Random rng_load();
//Black-Scholes analytical solution for comparison
tuple <double, double> black_scholes(double S0, double K, double T, double r, double sigma);

int n_price = 100000;
int n_bin = 1000;
int n_block = 100;
double l_block = n_price/n_block;

vector <double> call_direct, call_discrete;
vector <double> put_direct, put_discrete;

vector <double> m_call_direct, m_call_discrete;
vector <double> err_call_direct, err_call_discrete;
vector <double> m_put_direct, m_put_discrete;
vector <double> err_put_direct, err_put_discrete;

int main (int argc, char *argv[])
{
  Random rnd = rng_load();
  double S0 =100;
  double S_direct, S_discrete;
  double K=100;
  double r=0.1;
  double T=1;
  double sigma=0.25;

  //each bin has the same size, so the time intervals are equal
  double dt=T * 1./n_bin;

  for( int i=0; i < n_price; i++)
  {
    S_discrete=0;
    for( int j=0; j < n_bin; j++)
    {
      S_discrete *= exp((r - 0.5*pow(sigma,2))*dt + sigma * rnd.Gauss(0,1) *sqrt(dt));
    }
    S_direct = S0 * exp((r - 0.5*pow(sigma,2))*T + sigma * rnd.Gauss(0,T));

    call_direct.push_back(exp(-r * T) * max(0., S_direct - K));
    call_discrete.push_back(exp(-r * T) * max(0., S_discrete - K));
    put_direct.push_back(exp(-r * T) * max(0., K - S_direct));
    put_discrete.push_back(exp(-r * T) * max(0., K - S_discrete));
  }

  m_call_direct = get<0>(blocking_method(call_direct, n_price, n_block));
  err_call_direct = get<1>(blocking_method(call_direct, n_price, n_block));
  m_call_discrete = get<0>(blocking_method(call_discrete, n_price, n_block));
  err_call_discrete = get<1>(blocking_method(call_discrete, n_price, n_block));
  m_put_direct = get<0>(blocking_method(put_direct, n_price, n_block));
  err_put_direct = get<1>(blocking_method(put_direct, n_price, n_block));
  m_put_discrete = get<0>(blocking_method(put_discrete, n_price, n_block));
  err_put_discrete = get<1>(blocking_method(put_discrete, n_price, n_block));

/********************************
* Saving to file
********************************/

  ofstream c_dir_out("call_direct.dat");
  if (c_dir_out.is_open())
  {
    for (int i = 0; i < n_block; i++)
      c_dir_out << i + 1 << " " << m_call_direct[i] << " " << err_call_direct[i] << endl;
  }
  else
    cerr <<"Unable to open one or more output file: data saving failed" <<endl;
  c_dir_out.close();

  ofstream c_dis_out("call_discrete.dat");
  if (c_dis_out.is_open())
  {
    for (int i = 0; i < n_block; i++)
      c_dis_out << i + 1 << " " << m_call_discrete[i] << " " << err_call_discrete[i] << endl;
  }
  else
    cerr <<"Unable to open one or more output file: data saving failed" <<endl;
  c_dis_out.close();

  ofstream p_dir_out("put_price_direct.dat");
  if (p_dir_out.is_open())
  {
    for (int i = 0; i < n_block; i++)
      p_dir_out << i + 1 << " " << m_put_direct[i] << " " << err_put_direct[i] << endl;
  }
  else
    cerr <<"Unable to open one or more output file: data saving failed" <<endl;
  p_dir_out.close();

  ofstream p_dis_out("put_price_discrete.dat");
  if (p_dis_out.is_open())
  {
    for (int i = 0; i < n_block; i++)
      p_dis_out << i + 1 << " " << m_put_discrete[i] << " " << err_put_discrete[i] << endl;
  }
  else
    cerr <<"Unable to open one or more output file: data saving failed" <<endl;
  p_dis_out.close();

  double call_price, put_price;
  call_price = get<0>(black_scholes(100, 100, 1, 0.1, 0.25));
  put_price = get<1>(black_scholes(100, 100, 1, 0.1, 0.25));

  cout<<"===== Black-Scholes analytical solution ====="<<endl;
  cout<<"Call-option price = "<<call_price<<endl;
  cout<<"Put-option price = "<<put_price<<endl;

  return 0;
}

/***********************************
*   Function implementation
***********************************/
double N(double x)
{
    return 0.5 * (1. + erf(x / sqrt(2.)));
}

tuple <double, double> black_scholes(double S0, double K, double T, double r, double sigma)
{
    double d1 = 1./(sigma * sqrt(T)) * (log(S0 / K) + (r + pow(sigma,2) / 2.) * T);
    double d2 = d1 - sigma * sqrt(T);
    double C = S0 * N(d1) - K * exp(-r * T) * N(d2);
    double P = S0 *(N(d1) - 1.) - K * exp(-r * T) * (N(d2)-1.);
    return make_tuple(C,P);
}
/**********************************************/
double error(double sum, double sum2, int n)
{
  if (n == 0)
    return 0;
  else
    return sqrt((sum2-sum*sum)/n);
}

/**********************************************/
/* use tuple to return multiple values
***********************************************/

tuple<vector<double>,vector<double>> blocking_method(vector<double> v, int n_step, int n_blocks) {
  vector<double> av, av2, err;
  double s_block;
  double sum = 0, sum2 = 0;
  int l_block = n_step / n_block;

  // loop over blocks
  for (int i = 0; i < n_block; i++) {
    s_block = 0;
    // sum elements in one block
    for (int j = 0; j < l_block; ++j)
      s_block += v[i*l_block + j];

    // compute progressive sum
    sum += s_block/l_block;
    sum2 += pow(s_block/l_block, 2);

    // average over throws and evaluate error
    av.push_back(sum /(i + 1));
    av2.push_back(sum2 /(i + 1));
    err.push_back(error(av[i], av2[i], i));
  }
  return make_tuple(av, err);
}

/*********************************/

Random rng_load() {

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("rng/Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("rng/seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   rnd.SaveSeed();
   return rnd;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
