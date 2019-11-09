struct point
{
  double x;
  double y;
  double z;

  point(double a=0, double b=0, double c=0)
  : x(a), y(b), z(c)
  {
  }

  point& operator=(const point& a)
  {
    x=a.x;
    y=a.y;
    z=a.z;
    return *this;
  }

  point operator+(const point& a) const
  {
      return point(a.x + x, a.y + y, a.z + z);
  }
  point operator*(const double a) const
  {
      return point(x*a, y*a, z*a);
  }
  point operator/(const double a) const
  {
      return point(x*1./a, y*1./a, z*1./a);
  }
};
