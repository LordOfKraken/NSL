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

Random rng_load();

// Blocking method
void Measure(int);
void Reset(int);
void Accumulate(void);
void Averages(int);
double Integral(double);
double Error(double, double, int);

double walker_i, walker_u;
double glob_av_uni, glob_av_imp, glob_av2_uni, glob_av2_imp;
double blk_norm, blk_av_uni, blk_av_imp;

Random rnd = rng_load();
int n_block = 100;
int n_step = 1;

int main (int argc, char *argv[])
{
  int n_point = 1000;

  for( int iblk = 1; iblk <= n_block; iblk++ )
  {
    Reset(iblk);      // Reset block averages
    for(int istep=1; istep <= n_step; ++istep)
    {
      Measure(n_point);
      //Accumulate();     //Sum for each
    }
    Averages(iblk); // Evaluate block averages and print to file
  }
  return 0;
}

/***********************************
*   Function implementation
***********************************/

double Integral(double x)
{
  return M_PI_2*cos(M_PI_2*x);
}

void Measure(int n_point)
{
  double sample, s_uni;

  for( int i = 0; i < n_point; i++ )
  {
    sample=rnd.Importance();
    s_uni=rnd.Rannyu();
    walker_u += Integral(s_uni);
    walker_i += Integral(sample) /(2.0 * (1.0 - sample));
  }

  walker_u /= n_point;
  walker_i /= n_point;
  //cout << "walker_u " << walker_u << endl;
  //cout << "f( "<< s_uni <<" ) = " << Integral(s_uni)<< endl;

  blk_av_uni += walker_u;
  blk_av_imp += walker_i;
  blk_norm += 1;

  //cout << "blk av " << blk_av_uni << endl;
  //cout << "blk norm " << blk_norm << endl;
}

void Accumulate(void) //Update block averages
{
  blk_av_uni += walker_u/n_step;
  blk_av_imp += walker_i/n_step;
  blk_norm += 1;
  cout << "blk av " << blk_av_uni << endl;
}

void Averages(int iblk) //Print results for current block
{
  ofstream Imp,Uni;
  double stima_uni, stima_imp;
  double err_uni, err_imp;
  Uni.open("output/uniform.out",ios::app);
  Imp.open("output/importance.out", ios::app);

  cout << "Block number " << iblk << endl;

  stima_uni = blk_av_uni/blk_norm;
  cout << "stima "<< stima_uni << endl;
  glob_av_uni += stima_uni;
  cout << "globale " << glob_av_uni << endl;
  glob_av2_uni += stima_uni*stima_uni;
  err_uni = Error(glob_av_uni,glob_av2_uni,iblk);

  stima_imp = blk_av_imp/blk_norm;
  glob_av_imp += stima_imp;
  glob_av2_imp += stima_imp*stima_imp;
  err_imp = Error(glob_av_imp,glob_av2_imp,iblk);

  Uni << iblk << " " << glob_av_uni/iblk << " "<< err_uni << endl;
  Imp << iblk << " " << glob_av_imp/iblk << " "<< err_imp << endl;

  cout << iblk << " " << glob_av_uni/iblk << " "<< err_uni << endl;
  cout << iblk << " " << glob_av_imp/iblk << " "<< err_imp << endl;

  cout << "----------------------------" << endl << endl;

  Uni.close();
  Imp.close();
}

void Reset(int iblk) //Reset block averages
{
  if(iblk == 1)
  {
    glob_av_uni = 0;
    glob_av_imp = 0;
    glob_av2_uni = 0;
    glob_av2_imp = 0;
  }

  blk_av_uni = 0;
  blk_av_imp = 0;
  blk_norm = 0;

  walker_i=0;
  walker_u=0;
}

double Error(double sum, double sum2, int iblk)
{
  return sqrt((sum2/iblk - pow(sum/iblk,2))/iblk);
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
