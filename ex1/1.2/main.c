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
#include <cmath>
#include "rng/random.h"

using namespace std;

double error(double sum, double sum2, int n);
Random rng_load();

int n_rep = 10000;
double av_r, av_e, av_l;

int main (int argc, char *argv[])
{
  Random rnd = rng_load();

  //fill r with random numbers
  for(int n : {1, 2, 10, 100})
  {
    ofstream r_out("data/uniform_average_" + to_string(n) + ".dat");
    ofstream e_out("data/exponential_average_" + to_string(n) + ".dat");
    ofstream l_out("data/lorentz_average_" + to_string(n) + ".dat");

    for(int i = 0; i< n_rep; i++)
    {
      av_r=0;
      av_e=0;
      av_l=0;
      for( int j = 0; j < n; j++)
      {
        av_r += rnd.Rannyu()/n;
        av_e += rnd.Exp(1)/n;
        av_l += rnd.Lorentz(0,1)/n;
      }
      /*if (r_out.is_open() + l_out.is_open() +  e_out.is_open() <3) 
        std::cerr << "PROBLEM: Unable to open output file" << std::endl;
*/
      r_out << av_r << endl;
      e_out << av_e << endl;
      l_out << av_l << endl;
    }
    r_out.close();
    l_out.close();
    e_out.close();
  }
  return 0;
}

double error(double sum, double sum2, int n){

  if (n == 0)
    return 0;
  
  else
    return sqrt((sum2-sum*sum)/n);
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
