#ifndef myfunc_h
#define myfunc_h

class Function
{
  public:

  Function();
  Function(double a);
  ~Function();

  double Eval(double x);
  double Max();
  
  private:

  double _max;
};

#endif
