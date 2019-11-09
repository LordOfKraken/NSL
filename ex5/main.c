/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include "rng/random.h"

using namespace std;

double error(double, double, int);
Random rng_load();

double psi_1s(point);
double psi_2p(point);
double metropolis(double, double);

const double a0 = 5.29E-11;

int main(int argc, char *argv[])
{
  // Step calculation
  Random rnd = rng_load();
  int n_step=1E6;
  int n_block = 100;
  int l_block = n_step / n_block;
  point p_1s, new_1s;
  point p_2p, new_2p;
  double step_1s, step_2p;
  double r_1s;
  double r_2p;

  // Check acceptance ratio
  int count = 0;
  int good_1s = 0, good_2p = 0;

  // Blocking method
  double sum_1s = 0, sum2_1s = 0, av_1s = 0;
  double sum_2p = 0, sum2_2p = 0, av_2p = 0;

  // output files
  ofstream out_1s("r_ave_1s.dat");
  ofstream out_2p("r_ave_2p.dat");
  ofstream out_point_1s("point_1s.dat");
  ofstream out_point_2p("point_2p.dat");

  // Uniform
  step_1s = 1.35;
  step_2p = 3.53;

  for ( int j = 0; j < n_block; ++j)
  {
    p_1s = {0, 0, 0};
    p_2p = {0, 0, 0};
    av_1s = 0;
    av_2p = 0;

    for ( int i = 0; i < l_block; i++)
    {
      count++;
      // Uniform probability
      //new_1s = p_1s + rnd.Walk(step_1s);
      //new_2p = p_2p + rnd.Walk(step_2p);

      // Multivariate
      // the acceptance ratio changes to 0.40 for 1s and 0.70 for 2p
      new_1s = p_1s + rnd.WalkGauss(0.5,1.);
      new_2p = p_2p + rnd.WalkGauss(0.5,1.);

      if (rnd.Rannyu() <= metropolis(psi_1s(new_1s*a0), psi_1s(p_1s*a0)))
      {
        p_1s = new_1s;
        good_1s++;

        if (out_point_1s.is_open())
          out_point_1s << p_1s.x << " " << p_1s.y << " " << p_1s.z << endl;
        else
          cerr << "ERROR: Unable to open output file" << endl;
      }

      if (rnd.Rannyu() <= metropolis(psi_2p(new_2p*a0), psi_2p(p_2p*a0)))
      {
        p_2p = new_2p;
        good_2p++;

        if (out_point_2p.is_open())
          out_point_2p << p_2p.x << " " << p_2p.y << " " << p_2p.z << endl;
        else
          cerr << "ERROR: Unable to open output file" << endl;
      }

      r_1s = sqrt(pow(p_1s.x,2) + pow(p_1s.y,2) + pow(p_1s.z,2));
      r_2p = sqrt(pow(p_2p.x,2) + pow(p_2p.y,2) + pow(p_2p.z,2));

      av_1s += r_1s;
      av_2p += r_2p;
    }
    //cout << " Acceptance ratio ground state = " << good_1s * 1./ count << endl;
    //cout << " Acceptance ratio 2p state = " << good_2p * 1./ count << endl;

    av_1s /= l_block;
    sum_1s += av_1s;
    sum2_1s += pow(av_1s, 2);

    av_2p /= l_block;
    sum_2p += av_2p;
    sum2_2p += pow(av_2p, 2);

    if (out_1s.is_open())
      out_1s << j + 1 << " " << sum_1s * 1./(j + 1) << " " << error(sum_1s, sum2_1s, j+1) << endl;
    else
      cerr << "ERROR: Unable to open output file" << endl;

    if (out_2p.is_open())
      out_2p << j + 1 << " " << sum_2p * 1./(j + 1) << " " << error(sum_2p, sum2_2p, j+1) << endl;
    else
      cerr << "ERROR: Unable to open output file" << endl;
  }
  out_1s.close();
  out_2p.close();
  out_point_1s.close();
  out_point_2p.close();

  return 0;
}

/***********************************
*   Function implementation
***********************************/

double psi_1s(point p)
{
  double r = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
  return pow(a0, -3.) * M_1_PI * exp(-r *2. / a0);
}

/*****************************************/

double psi_2p(point p)
{
  double r = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
  double c_theta = p.x / sqrt(p.x*p.x + p.y*p.y);
  return pow(a0, -5.) * M_1_PI * pow(r,2) * exp(-r / a0) * pow(c_theta,2)/ 32.;
  //return pow(2. * a0, -5.*0.5) * r * exp(-r / (2. * a0)) * c_theta / sqrt(M_PI);
}

/***********************************************/

double metropolis(double p_new, double p_old) {
  return min(1., p_new / p_old);
}

/**********************************************/

double error(double ave, double ave2, int n)
{
  if (n == 1)
    return 0;
  else
    return sqrt((ave2/n-pow(ave/n,2))/(n-1));
}

/************************************************/

Random rng_load()
{
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
   if (input.is_open())
   {
      while ( !input.eof() )
      {
         input >> property;
         if( property == "RANDOMSEED" )
         {
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   }
   else cerr << "PROBLEM: Unable to open seed.in" << endl;

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
