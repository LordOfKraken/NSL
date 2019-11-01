#include "integral.h"
#include <algorithm>
#include <cmath>

using namespace std;

Integral::Integral(double a, double b, FunzioneBase *function)
{
	_integrand= function;
	_a = min(a,b);
	_b = max(a,b);
	if(a>b)		_sign= -1;
	else 		_sign=1;
}	

double Integral::midpoint(int nstep)
{
	_sum = 0.;
	_h = (_b - _a )/nstep;
	
	for(int i=0; i<nstep; i++)
	{
		double x= _a+(i+0.5)*_h;
		_sum+= _integrand->Eval(x);
	}
	
	_integral=_sign*_sum*_h;
	return _integral;
}

double Integral::simpson(int nstep)
{
	if(nstep%2!=0)
	{
		cout<<"il numero di intervalli e' stato aumentato di 1"<<endl;
		nstep=nstep+1;
	}
	
	_sum=0.;
	_h = (_b - _a )/nstep;
	
	double x=_a;
	
	_sum+=_integrand->Eval(x)/3.;
	
	for( int i=1;i<nstep;i++)
	{
		x=_a+i*_h;
		if(i%2==0)
		{
			_sum+=_integrand->Eval(x)*2/3;
		}
		else
		{
			_sum+=_integrand->Eval(x)*4/3;
		}
	}
	
	_sum+=_integrand->Eval(_b)/3;
	
	_integral=_sign*_sum*_h;
	return _integral;

}

double Integral::trapezi(int prec)
{
	int n=1;
	double sumt=0.;
	
	_sum=(_integrand->Eval(_a)-_integrand->Eval(_b))/2;
	
	
	do
	{
		n=2*n;
		_h=(_b-_a)/n;
		_sum=_sum+sumt;
		sumt=0.;
				
		for(int i=0;i<n;i++)
		{
			double x=_a+(i+0.5)*_h;
			sumt+=_integrand->Eval(x);		
		}
		
		cout<<"===== "<<sumt<<endl;
		
	}while((sumt-_sum)*_h>pow(10,-prec));
	

	_integral=_sign*_sum*_h;

	return _integral;
}


