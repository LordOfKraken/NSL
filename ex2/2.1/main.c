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
double function(double x);
Random rng_load();

vector <double> r_uni, m_u, err_u;
vector <double> r_imp, m_i, err_i;

int n_point = 1000;
int n_block = 100;
int n_rep = 100000;
double l_block = n_rep/n_block;
int x_min=0, x_max=1;

int main (int argc, char *argv[])
{
  Random rnd = rng_load();

  double sample;
  double uniform, importance;
  
  for( int i=0; i<n_rep; i++)
  {
    uniform=0, importance=0;
    for( int j=0; j< n_point;j++)
    {
      sample=rnd.Importance();
      uniform += M_PI_2 * cos(M_PI_2*rnd.Rannyu());
      //importance += M_PI_2 * cos(M_PI_2*sample)*5./(6*(1.-0.5*pow(sample,2)));
      importance += M_PI_2 * cos(M_PI_2*sample) * 0.5/(1-sample);
    }
    r_uni.push_back(uniform*1./n_point);
    r_imp.push_back(importance*1./n_point);
  }
  
  // use blocking method to evaluate average and error
  m_u = get<0>(blocking_method(r_uni, n_rep, n_block));
  m_i = get<0>(blocking_method(r_imp, n_rep, n_block));
  err_u = get<1>(blocking_method(r_uni, n_rep, n_block));
  err_i = get<1>(blocking_method(r_imp, n_rep, n_block));

/********************************
* Saving to file
********************************/

  ofstream u_out("uniform.dat");
  if (u_out.is_open()) 
  {
    for (int i = 0; i < n_block; i++)
      u_out << i * l_block << " " << m_u[i] << " " << err_u[i] << endl;
  } 
  else
    cerr <<"Unable to open integral output file: data saving failed" <<endl;
  u_out.close();
  
  ofstream i_out("importance.dat");
  if (i_out.is_open()) 
  {
    for (int i = 0; i < n_block; i++)
      i_out << i * l_block << " " << m_i[i] << " " << err_i[i] << endl;
  } 
  else
    cerr <<"Unable to open importance sampling output file: data saving failed" <<endl;
  i_out.close();
  return 0; 
}

/***********************************
*   Function implementation
***********************************/

double error(double sum, double sum2, int n){

  if (n == 0)
    return 0;
  
  else
    return sqrt((sum2-sum*sum)/n);
}


class Integral
{
  double evaluate(double x)
  {
    return M_PI_2*cos(M_PI_2*x);
  }

  private:

  double max=M_PI_2;
};

double function(double x)
{
  return M_PI_2*cos(M_PI_2*x);
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
