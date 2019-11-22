/****************************************************************
*****************************************************************
_/    _/  _/_/_/  _/       Numerical Simulation Laboratory
_/_/  _/ _/       _/       Physics Department
_/  _/_/    _/    _/       Universita' degli Studi di Milano
_/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <iomanip>      // output formatting
#include "MolDyn_NVE.h"

using namespace std;

int main()
{
  Input();             //Inizialization
  int nconf = 1;
  for( int iblk = 0; iblk < nblk; iblk++ )
  {
    Reset(iblk);      // Reset block averages

    for(int istep=1; istep <= nstep; ++istep)
    {
      Move();           //Move particles with Verlet algorithm

      if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
      //if(istep%mstep == 0)
      //{ //measure every mstep steps
        Measure();      //Properties measurement -> save in each walker[]
        Accumulate();   //Sum for each block
        //        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
        nconf += 1;
      //}
    }
    Averages(iblk); // Evaluate block averages and print to file
  }
  ConfFinal();         //Write final configuration to restart

  return 0;
}


void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf,ReadOldConf;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator

  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nblk;
  ReadInput >> nstep;
  ReadInput >> mstep;
  ReadInput >> iprint;
  ReadInput >> old;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Measurements are performed every " << mstep << " steps " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps per block = " << nstep << endl << endl;

  ReadInput.close();

  //Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  ip = 4; //Pressure
  n_props = 5; //Number of observables

  //measurement of g(r)
  igofr = n_props;
  nbins = 100;
  n_props = n_props + nbins;
  bin_size = (box/2.0)/(double)nbins;

  //Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

  switch(old)
  {
    case 0:
    //Prepare random initial velocities
    cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
    for (int i=0; i<npart; ++i){
      vx[i] = rand()/double(RAND_MAX) - 0.5;
      vy[i] = rand()/double(RAND_MAX) - 0.5;
      vz[i] = rand()/double(RAND_MAX) - 0.5;
    }
    break;

    case 1:
    // Read old configuration
    cout << "Read old configuration from file old.0 " << endl << endl;
    ReadOldConf.open("old.0");
    for (int i=0; i<npart; ++i){
      ReadOldConf >> xold[i] >> yold[i] >> zold[i];
      xold[i] = xold[i] * box;
      yold[i] = yold[i] * box;
      zold[i] = zold[i] * box;
    }
    ReadOldConf.close();
    cout << "Calculate velocities with Verlet" << endl << endl;
    Move();
    break;
  }

  //scale velocities to match temperature
  double sumv[3] = {0.0, 0.0, 0.0};
  for (int i=0; i<npart; ++i){
    sumv[0] += vx[i];
    sumv[1] += vy[i];
    sumv[2] += vz[i];
  }
  for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;

  double sumv2 = 0.0, fs;
  for (int i=0; i<npart; ++i){
    vx[i] = vx[i] - sumv[0];
    vy[i] = vy[i] - sumv[1];
    vz[i] = vz[i] - sumv[2];

    sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
  }
  sumv2 /= (double)npart;

  fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
  for (int i=0; i<npart; ++i){
    vx[i] *= fs;
    vy[i] *= fs;
    vz[i] *= fs;
    //estimate old conf
    xold[i] = Pbc(x[i] - vx[i] * delta);
    yold[i] = Pbc(y[i] - vy[i] * delta);
    zold[i] = Pbc(z[i] - vz[i] * delta);
  }
  return;
}

void Move(void){ //Move particles with Verlet algorithm, calculate vx,vy,vz
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  return f;
}

void Measure(){ //Properties measurement
  double v, t, p, vij, pij;
  double dx, dy, dz, dr;

  v = 0.0; //reset observables
  t = 0.0;
  p = 0.0;

  //cycle over pairs of particles
  for (int i=0; i<npart-1; ++i)
  {
    for (int j=i+1; j<npart; ++j)
    {

      dx = Pbc( x[i] - x[j] );
      dy = Pbc( y[i] - y[j] );
      dz = Pbc( z[i] - z[j] );

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      //filling g(r) histogram;
      walker[int(igofr + trunc(dr/bin_size))] += 2;

      if(dr < rcut)
      {
        //Potential energy
        vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
        v += vij;
        // Pressure
        pij = 48 * (pow(1 / dr, 12) - 0.5 * pow(1 / dr, 6));
        p += pij;
      }
    }
  }

  for (int k = igofr; k < igofr + nbins; ++k)
  {
    walker[k] /= rho * npart * (4. * M_PI / 3.) * (pow(bin_size * (k - igofr + 1), 3) - pow(bin_size * (k - igofr), 3));
  }


  //Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

  walker[iv] = v/(double)npart;
  walker[ik] = t/(double)npart;
  walker[ie] = (v+t)/(double)npart;
  walker[it] = (2.0 / 3.0) * t/(double)npart;
  walker[ip] = (rho * t + 1./3. * p)/(double)npart;
  /*
  stima_pot = v/(double)npart; //Potential energy per particle
  stima_kin = t/(double)npart; //Kinetic energy per particle
  stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
  stima_etot = (t+v)/(double)npart; //Total energy per particle
  stima_pres = (rho * t + 1. / 3. * p) / (double)npart; //Pressure per particle

  Epot << stima_pot  << endl;
  Ekin << stima_kin  << endl;
  Temp << stima_temp << endl;
  Etot << stima_etot << endl;
  Pres << stima_pres << endl;
  */
  return;
}

void Accumulate(void) //Update block averages
{

  for(int i=0; i<n_props; ++i)
  {
    blk_av[i] = blk_av[i] + walker[i];
  }
  blk_norm += mstep; //measure and accumulate every 10 step (mstep)
}

void Averages(int iblk) //Print results for current block
{
  double r, gdir;
  ofstream Epot, Ekin, Etot, Temp, Pres, Gofr, Gave;
  int blk=iblk+1;

  Epot.open("output/epot.out",ios::app);
  Ekin.open("output/ekin.out",ios::app);
  Temp.open("output/temp.out",ios::app);
  Etot.open("output/etot.out",ios::app);
  Pres.open("output/pres.out",ios::app);
  Gofr.open("output/gofr.out",ios::app);
  Gave.open("output/gave.out",ios::app);
  const int wd=12;

  cout << "Block number " << iblk << endl;
  //cout << "Acceptance rate " << accepted/attempted << endl << endl;

  stima_pot = blk_av[iv]/blk_norm * (double)npart; //Potential energy
  glob_av[iv] += stima_pot;
  glob_av2[iv] += stima_pot*stima_pot;
  err_pot = Error(glob_av[iv],glob_av2[iv],blk);

  stima_kin = blk_av[ik]/blk_norm * (double)npart; //Kinetic energy
  glob_av[ik] += stima_kin;
  glob_av2[ik] += stima_kin*stima_kin;
  err_kin = Error(glob_av[ik],glob_av2[ik],blk);

  stima_temp = blk_av[it]/blk_norm * (double)npart; //Kinetic energy
  glob_av[it] += stima_temp;
  glob_av2[it] += stima_temp*stima_temp;
  err_temp = Error(glob_av[it],glob_av2[it],blk);

  stima_etot = blk_av[ie]/blk_norm * (double)npart; //Total energy
  glob_av[ie] += stima_etot;
  glob_av2[ie] += stima_etot*stima_etot;
  err_etot = Error(glob_av[ie],glob_av2[ie],blk);

  stima_pres = blk_av[ip]/blk_norm * (double)npart; //Pressure
  glob_av[ip] += stima_pres;
  glob_av2[ip] += stima_pres*stima_pres;
  err_pres = Error(glob_av[ip],glob_av2[ip],blk);

  //Potential energy per particle
  Epot << setw(wd) << iblk << setw(wd) << glob_av[iv]/(double)blk << setw(wd) << err_pot << endl;

  //Kinetic energy per particle
  Ekin << setw(wd) << iblk << setw(wd) << glob_av[ik]/(double)blk << setw(wd) << err_kin << endl;

  //Temperature
  Temp << setw(wd) << iblk << setw(wd) << glob_av[it]/(double)blk << setw(wd) << err_temp << endl;

  //Total energy per particle
  Etot << setw(wd) << iblk << setw(wd) << glob_av[ie]/(double)blk << setw(wd) << err_etot << endl;

  //Pressure
  Pres << setw(wd) << iblk << setw(wd) << glob_av[ip]/(double)blk << setw(wd) << err_pres << endl;

  //g(r)

  for (int i = igofr; i < igofr + nbins; i++)
  {
    stima_g = blk_av[i] / blk_norm;
    glob_av[i] += stima_g;
    glob_av2[i] += pow(stima_g,2);
    err_g = Error(glob_av[i], glob_av2[i], iblk+1);
    Gofr << setw(wd) << iblk << setw(wd) << i - igofr << setw(wd) << stima_g << endl;
    if( iblk == nblk - 1)
    {
      Gave << setw(wd) << i - igofr << setw(wd) << glob_av[i]/(double)blk << setw(wd) << err_g << endl;
    }
  }
  cout << "----------------------------" << endl << endl;

  Epot.close();
  Ekin.close();
  Temp.close();
  Etot.close();
  Pres.close();
  Gofr.close();
  Gave.close();
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
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

void ConfFinal(void){ //Write final configuration
  ofstream WriteConf, WriteOldConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();

  if(old == 1)
  {
    cout << "Print old final configuration to file old.final " << endl << endl;
    WriteOldConf.open("old.final");

    for (int i=0; i<npart; ++i){
      WriteOldConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
    }
    WriteOldConf.close();
  }
  return;
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
  return r - box * rint(r/box);
}


/**********************************************/
double Error(double sum, double sum2, int iblk)
{
  return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
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
