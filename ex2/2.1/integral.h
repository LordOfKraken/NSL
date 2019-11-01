#ifndef integral_h
#define integral_h

#include "funzionebase.h"

class Integral {

	public:
		Integral(double a, double b, FunzioneBase *function);
		double midpoint(int nstep);
		double simpson(int nstep);
		double trapezi(int prec);
		
	private:
		double _a, _b;
		double _sum;
		double _h;
		int _sign;
		double _integral;
		FunzioneBase *_integrand;
		
};

#endif 
