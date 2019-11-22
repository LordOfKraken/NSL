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

double error(double sum, double sum2, int n);
Random rng_load();
tuple<vector<double>,vector<double>> blocking_method(vector<double> v, int n_throw, int n_block);

int n_throw = 100000;
int n_block = 100;
int l_block = n_throw/n_block;

vector <double> r, m_r, err_r, s, m_s, err_s;

int main (int argc, char *argv[])
{
  Random rnd = rng_load();

//  fill r with random numbers

  for(int i = 0; i < n_throw; i++) {
      r.push_back(rnd.Rannyu());
      s.push_back(pow(r[i]-0.5,2));
  }

  // use blocking method to evaluate average and error
  m_r = get<0>(blocking_method(r, n_throw, n_block));
  err_r = get<1>(blocking_method(r, n_throw, n_block));

  m_s = get<0>(blocking_method(s, n_throw, n_block));
  err_s = get<1>(blocking_method(s, n_throw, n_block));

/*****************************
* Chi squared evaluation
*****************************/

  int n_interval=100;
  int n_rep=100;
  vector <int> count;
  vector <double> chi, chi_sum;
  count.resize(n_interval);
  chi.resize(n_rep);
  chi_sum.resize(n_rep);
  double chi_avg=0;

// enlarge r<> to 1E6
  n_throw=1E6;
  r.resize(0);

  for( int i = 0; i < n_throw; i++)
  {
    r.push_back(rnd.Rannyu());
  }

  for(int i = 0; i < n_rep; i++)
  {
    for( int j = 0; j < n_throw/n_rep; j++)
    {
      count[int(r[i*n_throw/n_rep + j]*100)] ++;
    }

    for( int k = 0; k< n_interval; k++)
    {
      chi[i]+=pow(count[k]-n_throw/(n_rep*n_interval),2)/(n_throw/(n_rep*n_interval));
      count[k]=0;
    }
    chi_avg+=chi[i];
    chi_sum[i]=chi_avg/i;
  }
  cout<< "Chi squared average = " << chi_avg/n_rep << endl;

/********************************
* Saving to file
********************************/

  ofstream r_out("r.dat");
  if (r_out.is_open())
  {
    for (int i = 0; i < n_block; i++)
      r_out << i * l_block << " " << m_r[i] - 0.5 << " " << err_r[i] << endl;
  }
  else
    cerr <<"Unable to open r output file: data saving failed" <<endl;
  r_out.close();

  ofstream s_out("s.dat");
  if (s_out.is_open())
  {
    for (int i = 0; i < n_block; i++)
      s_out << i * l_block << " " << m_s[i] - 1./12 << " " << err_s[i] << endl;
  }
  else
    cerr <<"Unable to open sigma output file: data saving failed" <<endl;
  s_out.close();

  ofstream chi_out("chi.dat");
  if (chi_out.is_open())
  {
    for (int i = 0; i < n_interval; i++)
      chi_out << i  << " " << chi[i] << endl;
  }
  else
    cerr <<"Unable to open chi output file: data saving failed" <<endl;
  chi_out.close();

  ofstream chis_out("chi_avg.dat");
  if (chis_out.is_open())
  {
    for (int i = 0; i < n_interval; i++)
      chis_out << i  << " " << chi_sum[i] << endl;
  }
  else
    cerr <<"Unable to open chi average output file: data saving failed" <<endl;
  chis_out.close();
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

/**********************************************/
/* tuple are the common way to return multiple values, of different types
*/

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
