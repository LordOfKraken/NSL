#include "function.h"
#include <cmath>

using namespace std;

Function::Function()
{
  _max=M_PI_2;
}

Function::~Function(){}

Function::Function(double a)
{
  _max=a;
}

double Function::Eval(double x)
{
  return M_PI_2*cos(M_PI_2*x);
}

double Function::Max()
{
  return _max;
}
