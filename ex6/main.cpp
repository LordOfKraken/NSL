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
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{
  Input(); //Inizialization
  while(t <= tstep)
  {
    //fix temperature with 2 decimal digits to solve approximation problem
    temp = round(100*(tmin + t*(tmax-tmin)/tstep))*0.01;
    cout << "Temperature = " << temp <<endl;

    beta=1/temp;

    for(int term=1; term <= nstep; ++term) // Thermalization
          Move(metro);

    for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
    {
      Reset(iblk);   //Reset block averages
      for(int istep=1; istep <= nstep; ++istep)
      {
        Move(metro);
        Measure();
        Accumulate(); //Update block averages
      }
      Averages(iblk);   //Calculate results for current block
    }
    t++;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    cout << "----------------------------" << endl << endl;
  }

  ConfFinal(); //Write final configuration
  return 0;
}


void Input(void)
{
  ifstream ReadInput;
  string oldspin = "conf.old";

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

  rnd=rng_load();

  //Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;

  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  ReadInput >> old;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  if(old == 1) cout << "Loading an old spin configuration" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;

  ReadInput.close();

  //Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility

  n_props = 4; //Number of observables

  if(old == 0)
  {
    //new configuration
    for (int i=0; i<nspin; ++i)
    {
      //s[i]=trunc(rnd.Rannyu()*2 + 1)
      if(rnd.Rannyu() >= 0.5) s[i] = 1;
      else s[i] = -1;
    }
  }
  else Reload(oldspin);

  //Evaluate energy etc. of the initial configuration
  Measure();

  //Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}

void Reload(string load_file)
{
  ifstream SpinConf(load_file);
  for (int i=0; i<nspin; ++i)
  {
    SpinConf >> s[i];
  }
}

void Move(int metro)
{
  int o;
  //double p, energy_old, energy_new, sm;
  //double energy_up, energy_down;
  double delta_energy;

  for(int i=0; i<nspin; ++i)
  {
    //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);
    delta_energy = 2 * s[o] * (J * (s[Pbc(o - 1)] + s[Pbc(o + 1)]) + h);

    if(metro==1) //Metropolis
    {
      attempted++;
      if( delta_energy < 0 || rnd.Rannyu() < exp(-beta * delta_energy) )
      {
        s[o] *= -1;
        accepted++;
      }
    }
    else //Gibbs sampling
    {
      attempted = 1;
      accepted = 1;
      //s[o] = 1;
      delta_energy = 2 * s[o] * (J * (s[Pbc(o - 1)] + s[Pbc(o + 1)]) + h);
      if (rnd.Rannyu() > (1. / (1. + exp(-beta * delta_energy))))
        s[o] *= -1;
    }
  }
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure()
{
  int bin;
  double u = 0.0, m = 0.0;

  //cycle over spins
  for (int i=0; i<nspin; ++i)
  {
    u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
    m += s[i];
  }
  walker[iu] = u;
  walker[ic] = pow(u,2);
  walker[im] = m;
  walker[ix] = pow(m,2);
}


void Reset(int iblk) //Reset block averages
{

  if(iblk == 1)
  {
    for(int i=0; i<n_props; ++i)
    {
      glob_av[i] = 0;
      glob_av2[i] = 0;
    }
  }

  for(int i=0; i<n_props; ++i)
  {
    blk_av[i] = 0;
  }
  blk_norm = 0;
  attempted = 0;
  accepted = 0;
}


void Accumulate(void) //Update block averages
{

  for(int i=0; i<n_props; ++i)
  {
    blk_av[i] = blk_av[i] + walker[i];
  }
  blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
  ofstream Ene, Heat, Mag, Chi;
  const int wd=12;

  //cout << "Block number " << iblk << endl;
  //cout << "Acceptance rate " << accepted/attempted << endl << endl;

  //Energy
  stima_u = blk_av[iu]/blk_norm/(double)nspin;
  glob_av[iu]  += stima_u;
  glob_av2[iu] += stima_u*stima_u;
  err_u=Error(glob_av[iu],glob_av2[iu],iblk);

  // Heat Capacity
  stima_c = pow(beta,2)*(blk_av[ic]/(blk_norm*nspin) - pow(stima_u,2)*nspin);
  glob_av[ic] += stima_c;
  glob_av2[ic] += stima_c*stima_c;
  err_c=Error(glob_av[ic],glob_av2[ic],iblk);

  // Magnetization
  stima_m = blk_av[im]/(blk_norm*nspin);
  glob_av[im] += stima_m;
  glob_av2[im] += stima_m*stima_m;
  err_m=Error(glob_av[im],glob_av2[im],iblk);

  // Magnetic susceptibility
  stima_x = beta*blk_av[ix]/(blk_norm*nspin);
  glob_av[ix] += stima_x;
  glob_av2[ix] += stima_x*stima_x;
  err_x=Error(glob_av[ix],glob_av2[ix],iblk);

  string rtype="gibbs";
  if(metro==1)
    rtype="metropolis";

  if(iblk == nblk)
  {
    Ene.open("energy_" + rtype + ".out",ios::app);
    Ene << setw(wd) << t <<   setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    Ene.close();

    Heat.open("heat_capacity_" + rtype + ".out",ios::app);
    Heat << setw(wd) << t << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
    Heat.close();

    Mag.open("magnetization_" + rtype + ".out",ios::app);
    Mag << setw(wd) << t <<  setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
    Mag.close();

    if(h!=0)
    {
      Chi.open("susceptibility_" + rtype + ".out",ios::app);
      Chi << setw(wd) << t << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
      Chi.close();
    }
  }
  //cout << "----------------------------" << endl << endl;
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
  if(i >= nspin) i = i - nspin;
  else if(i < 0) i = i + nspin;
  return i;
}

double Error(double sum, double sum2, int iblk)
{
  return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

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
