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
#include <tuple>        // blocking_method
#include <vector>       // blocking_method
#include "md.h"

using namespace std;

int main(){
  Input();             //Inizialization
  int nconf = 1;
  for(int istep=1; istep <= nstep; ++istep){
    Move();           //Move particles with Verlet algorithm
    if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
    if(istep%10 == 0){
      Measure();     //Properties measurement
      //        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
      nconf += 1;
    }
  }
  if(nblock>0 && nblock<nstep)
    Blocking();

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
  ReadInput >> nstep;
  ReadInput >> nblock;
  ReadInput >> iprint;
  ReadInput >> old;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  cout << "Number of blocks = " << nblock << endl << endl;
  ReadInput.close();

  //Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

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
  int bin;
  double v, t, p, vij, pij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp, Pres;

  Epot.open("output/output_epot.dat",ios::app);
  Ekin.open("output/output_ekin.dat",ios::app);
  Temp.open("output/output_temp.dat",ios::app);
  Etot.open("output/output_etot.dat",ios::app);
  Pres.open("output/output_pres.dat",ios::app);

  v = 0.0; //reset observables
  t = 0.0;
  p = 0.0;

  //cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

      dx = Pbc( x[i] - x[j] );
      dy = Pbc( y[i] - y[j] );
      dz = Pbc( z[i] - z[j] );

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      if(dr < rcut){
        //Potential energy
        vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
        v += vij;
        // Pressure
        pij = 48 * (pow(1 / dr, 12) - 0.5 * pow(1 / dr, 6));
        p += pij;
      }
    }
  }

  //Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

  stima_epot = v/(double)npart; //Potential energy per particle
  stima_ekin = t/(double)npart; //Kinetic energy per particle
  stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
  stima_etot = (t+v)/(double)npart; //Total energy per particle
  stima_pres = (rho * t + 1. / 3. * p) / (double)npart; //Pressure per particle

  v_epot.push_back(stima_epot);
  v_ekin.push_back(stima_ekin);
  v_etot.push_back(stima_ekin + stima_epot);
  v_temp.push_back(stima_temp);
  v_pres.push_back(stima_pres);

  Epot << stima_epot  << endl;
  Ekin << stima_ekin  << endl;
  Temp << stima_temp << endl;
  Etot << stima_etot << endl;
  Pres << stima_pres << endl;

  Epot.close();
  Ekin.close();
  Temp.close();
  Etot.close();
  Pres.close();

  return;
}

void ConfFinal(void){ //Write final configuration
  ofstream WriteConf, WriteOldConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();

  //if(old == 1)
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

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
  return r - box * rint(r/box);
}

void Blocking(void){
  vector<double> m_epot, m_ekin, m_etot, m_temp, m_pres;
  vector<double> err_epot, err_ekin, err_etot, err_temp, err_pres;
  cout<< "roba" <<endl;
  m_epot=get<0>(blocking_method(v_epot,nstep*0.1,nblock));
    cout<< "roba" <<endl;
  m_ekin=get<0>(blocking_method(v_ekin,nstep*0.1,nblock));
  m_etot=get<0>(blocking_method(v_etot,nstep*0.1,nblock));
  m_temp=get<0>(blocking_method(v_temp,nstep*0.1,nblock));
  m_pres=get<0>(blocking_method(v_pres,nstep*0.1,nblock));

  err_epot=get<1>(blocking_method(v_epot,nstep*0.1,nblock));
  err_ekin=get<1>(blocking_method(v_ekin,nstep*0.1,nblock));
  err_etot=get<1>(blocking_method(v_etot,nstep*0.1,nblock));
  err_temp=get<1>(blocking_method(v_temp,nstep*0.1,nblock));
  err_pres=get<1>(blocking_method(v_pres,nstep*0.1,nblock));

  ofstream ave_epot("output/ave_epot.out");
  if (ave_epot.is_open())
  {
    for (int i = 0; i < nblock; i++)
      ave_epot << i + 1 << " " << m_epot[i] << " " << err_epot[i] << endl;
  }
  ave_epot.close();

  ofstream ave_ekin("output/ave_ekin.out");
  if (ave_ekin.is_open())
  {
    for (int i = 0; i < nblock; i++)
      ave_ekin << i + 1 << " " << m_ekin[i] << " " << err_ekin[i] << endl;
  }
  ave_ekin.close();

  ofstream ave_etot("output/ave_etot.out");
  if (ave_etot.is_open())
  {
    for (int i = 0; i < nblock; i++)
      ave_etot << i + 1 << " " << m_etot[i] << " " << err_etot[i] << endl;
  }
  ave_etot.close();

  ofstream ave_temp("output/ave_temp.out");
  if (ave_temp.is_open())
  {
    for (int i = 0; i < nblock; i++)
      ave_temp << i + 1 << " " << m_temp[i] << " " << err_temp[i] << endl;
  }
  ave_temp.close();

  ofstream ave_pres("output/ave_pres.out");
  if (ave_pres.is_open())
  {
    for (int i = 0; i < nblock; i++)
      ave_pres << i + 1 << " " << m_pres[i] << " " << err_pres[i] << endl;
  }
  ave_pres.close();
  return;
}

/**********************************************/
double error(double ave, double ave2, int n)
{
  if (n == 0)
    return 0;
  else
    return sqrt((ave2-ave*ave)/n);
}

tuple<vector<double>,vector<double>> blocking_method(vector<double> v, int n_step, int n_block) {
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

/****************************************************************
*****************************************************************
_/    _/  _/_/_/  _/       Numerical Simulation Laboratory
_/_/  _/ _/       _/       Physics Department
_/  _/_/    _/    _/       Universita' degli Studi di Milano
_/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
