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
#include <cmath>
#include "rng/random.h"

using namespace std;

int drop(double L, double d, Random *r);
Random rng_load();

// Blocking method
void Measure(int);
void Reset(int);
void Accumulate(void);
void Averages(int);
double Error(double, double, int);

Random rnd=rng_load();

double L=0.7;
double d=1;
int n_block = 100;
int n_step = 1;

double walker;
double glob_av,glob_av2;
double blk_norm, blk_av;

int main()
{
  double n_throw=10000;

  for( int iblk = 1; iblk <= n_block; iblk++ )
  {
    Reset(iblk);      // Reset block averages
    for( int istep = 1; istep <= n_step; ++istep )
    {
      Measure(n_throw); //Pi measurement -> save in walker
      Accumulate();     //Sum walker for each block
    }
    Averages(iblk); // Evaluate block averages and print to file
  }
  return 0;
}

int drop(double L, double d)
{
  // X coord of first needle end
  double x1=rnd.Rannyu(0,d);
  // cos_t have a non-uniform distribution in [-1,+1]
  double x2, y2, cos_t;
  do{
    x2=rnd.Rannyu(-1,1);
    y2=rnd.Rannyu(-1,1);
  }while(pow(x2,2) + pow(y2,2) > 1);

  cos_t = x2 * pow(x2*x2 + y2*y2, -0.5);

  if(x1 + L*cos_t >= d || x1 + L*cos_t <= 0)
    return 1;
  else
    return 0;
}


void Measure(int n_throw){ //Properties measurement
  double n_hit;

  n_hit=0;
  for( int i=0; i< n_throw; i++)
  {
    n_hit+=drop(L,d);
  }
  walker = L*n_throw*2./(d*n_hit);
  return;
}

void Accumulate(void) //Update block averages
{
  blk_av += walker;
  blk_norm += 1;
}

void Averages(int iblk) //Print results for current block
{
  ofstream Pi;
  double stima_pi, err_pi;
  Pi.open("output/pi.out",ios::app);

  cout << "Block number " << iblk << endl;

  stima_pi = blk_av/blk_norm;
  glob_av += stima_pi;
  glob_av2 += stima_pi*stima_pi;
  err_pi = Error(glob_av,glob_av2,iblk);

  //Potential energy per particle
  Pi << iblk << " " << glob_av/(double)iblk << " " << err_pi << endl;

  cout << " Pi = " << glob_av/(double)iblk << " +- " << err_pi << endl;
  cout << "----------------------------" << endl << endl;

  Pi.close();
}

void Reset(int iblk) //Reset block averages
{
  if(iblk == 1)
  {
    glob_av = 0;
    glob_av2 = 0;
  }

  blk_av = 0;
  blk_norm = 0;
  walker = 0;
}

double Error(double sum, double sum2, int iblk)
{
  return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
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
