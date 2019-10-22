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

int n_throw = 100000;
int n_block = 100;
int l_block = n_throw/n_block;

// _r variables are related to the average of r requested 
// _s variables are related to the sigma requested

vector <double> av_r, av2_r, av_s, av2_s;
vector <double> r, m_r, err_r, m_s, err_s;

int main (int argc, char *argv[])
{
  Random rnd = rng_load();

  //fill r with random numbers
  for(int i = 0; i < n_throw; i++) {
      r.push_back(rnd.Rannyu());
  }
  
  //average on single blocks, and save in vector
  for( int i = 0; i < n_block; i++)
  {
    double sum_r = 0, sum_s = 0;
    for(int j = 0; j < l_block; j++)
    {
      sum_r += r[i*l_block + j];
      sum_s += pow(r[i*l_block + j] - 0.5,2);
    }
    av_r.push_back(sum_r/l_block);
    av2_r.push_back(pow(sum_r/l_block,2));
    
    av_s.push_back(sum_s/l_block);
    av2_s.push_back(pow(sum_s/l_block,2));
  }
  

  //computing cumulative average and error
  for( int i = 0; i < n_block; i++)
  {
    double sum_r = 0, sum2_r = 0;
    double sum_s = 0, sum2_s = 0;
    for(int j = 0; j < i+1; j++)
    {
      sum_r += av_r[j];
      sum2_r += av2_r[j];
      
      sum_s += av_s[j];
      sum2_s += av2_s[j];
    }
    m_r.push_back(sum_r/(i+1));
    err_r.push_back(error(sum_r/(i+1), sum2_r/(i+1),i+1));
    
    m_s.push_back(sum_s/(i+1));
    err_s.push_back(error(sum_s/(i+1), sum2_s/(i+1),i+1));
  }

  //saving to file
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
