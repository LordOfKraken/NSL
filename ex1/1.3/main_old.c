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

int drop(double L, double d, Random *r);
double error(double sum, double sum2, int n);
tuple<vector<double>,vector<double>> blocking_method(vector<double> v, int n_step, int n_block);
Random rng_load();

int main(int argc, char *argv[])
{
  Random rnd=rng_load();
  double L=0.5;
  double d=1;
  int n_hit=0;
  int n_throw = 10000;
  int n_block = 100;
  int n_rep = 10000;
  double l_block = n_throw*1./n_block;
  vector <double> pi, ave, err;

  for( int j=0; j< n_rep; j++)
  {
    n_hit=0;
    for( int i=0; i< n_throw; i++)
    {
      n_hit+=drop(L,d,&rnd);
    }
    pi.push_back(2*L*n_throw/(d*n_hit));
  }

  ave=get<0>(blocking_method(pi,n_rep,n_block));
  err=get<1>(blocking_method(pi,n_rep,n_block));

  ofstream out("pi.dat");
  if (out.is_open())
    for (int i = 0; i < n_block; ++i)
      out << i << " " << pi[i] << " " << err[i] << endl;
  else
    cerr << "PROBLEM: Unable to open output file" << endl;

  out.close();
  return 0;
}

int drop(double L, double d, Random *r)
{
  // X coord of first needle end
  double x1=r->Rannyu(0,d);
  // cos_t have a non-uniform distribution in [-1,+1]
  double x2, y2, cos_t;
  do{
    x2=r->Rannyu(-L,L);
    y2=r->Rannyu(-L,L);
  }while(pow(x2,2)+pow(y2,2)>pow(L,2));

  cos_t = x2*1./sqrt(x2*x2+y2*y2);

  if(x1 + L*cos_t > d || x1 + L*cos_t < 0)
    return 1;
  else
    return 0;
}

double error(double sum, double sum2, int n)
{
  if (n == 0)
    return 0;
  else
    return sqrt((sum2-sum*sum)/n);
}

tuple<vector<double>,vector<double>> blocking_method(vector<double> v, int n_step, int n_block)
{
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
