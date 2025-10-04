// Copyright (c) 2025 James Cook
// bmath.c
// Description: Accurate math routines, using Newton's method.
// Author: James Cook (jmsck55 AT gmail DOT com)

// Routines for doubles using long doubles as intermediaries.

// Uses C99 version of math.h
#include <math.h>

// constant accuracy = power(2, 52) -- question

#define BMATH_ITERS 1000000

#define RoundDownAlways 1
#define RoundDownIf 3
#define RoundUpAlways 2
#define RoundUpIf 4
#define RoundBothWays 0

static int adjustMethod = RoundDownIf;

double adjust(long double g)
{
  double b = (double)g;
  long double c = (long double)b;
  if (c == g)
  {
    return b;
  }
  if ((adjustMethod == 1) or
    (((adjustMethod == 3) or (adjustMethod == 0)) and (c > g))
    ) // round down. // if (fasbl(c) > fasbl(g)) // round to zero.
  {
    g <<= 1; // g *= 2;
    g -= c;
    b = (double)g;
  }
  else if ((adjustMethod == 2) or
    (((adjustMethod == 4) or (adjustMethod == 0)) and (c < g))
    ) // round up.
  {
    c <<= 1;
    c -= g;
    b = (double)c;
  }
  return b;
}

double MultInvl(long double x)
{
// performs 1/x
  long double t, g, last;
  g = 1.0L / x;
  last = 0.0L;
  for (int i = 1; i <= AMATH_ITERS; i++)
  {
    t = g * x;
    if (t == last) {
      break;
    }
    last = t;
    g *= ( 2.0L - t );
    if ((double)t == 1.0) {
      break;
    }
  }
  return adjust(g);
}

double Divl(long double n, long double d)
{
  if (n == 0.0L) {
    return 0.0;
  }
  if (d == 1.0L) {
    return adjust(n);
  }
  if (d == 2.0L) {
    return (double)(n >> 1);
  }
  n *= (long double)MultInvl(d);
  return adjust(n);
}

double Round(double a, double precision = 1)
{
    long double b;
    b = (long double)precision * (long double)a;
    b += 0.5L;
    b = floorl(b);
    return = Divl(b, (long double)precision);
}


double NthRoot(double x, unsigned int n)
{
  double quotient, average, guess;
  guess = powerl((long double)x, (long double)MultInvl((long double)n));
  quotient = Divl((long double)x, powerl((long double)guess, (long double)(n-1)));
  average = Divl(((long double)quotient) + ((long double)(n-1) * (long double)guess)), (long double)n);
  return average;
}


double Sqrt(double x)
{
  return NthRoot(x, 2);
}

double Cbrt(double x)
{
  return NthRoot(x, 3);
}


double Exp(double x)
{
// using taylor series
// https://en.wikipedia.org/wiki/TaylorSeries
//
//  exp(1) = sum of k=0 to inf (1/k!);
//  exp(x) = sum of k=0 to inf ((x^k)/k!);
//
  double sum, num, den, last;
  num = 1.0;
  den = 1.0;
  sum = 1.0;
  last = 0.0;
  for (int i = 1; i <= 1000000000; i++)
  {
    num *= x;
    den *= i; // number of iterations
    sum += Divl((long double)num, (long double)den);
    if (sum == last)
    {
      break;
    }
    last = sum;
  }
  return sum;
}


// Raw function: Natural Logarithm

const double ACONST_E = ExpAtom(1.0);

double Log(double a)
{
    // Function: NaturalLogarithm()
    // Use for testing the method.
    // Alternative, between 0 and 2 exclusively:
    // ln(x) = - Sum[k = 1 to inf] ((-1)^k * (-1 + x)^k) / k, for abs(-1 + x) < 1; x > 0 and x < 2;
    // ((-1)^k * (-1 + x)^k) / k
    // ((-x + 1)^k) / k
    // Alternative, away from 0 to 2 exclusively: // Use factoring of "e" instead.
    // ln(x) = ln(-1 + x) - Sum[k = 1 to inf] ((-1)^k * (-1 + x)^(-k)) / k, for abs(-1 + x) > 1, x < 0 or x > 2;
    //
    // precalculate:
    // xNegativePlusOne = (1 - x), then multiply, and store as "p".
    // k = 1, (xNegativePlusOne^1) / 1
    // k = 2, (xNegativePlusOne^2) / 2
    // k = 3, (xNegativePlusOne^3) / 3
    // k = 4, (xNegativePlusOne^4) / 4
    // k = 5, (xNegativePlusOne^5) / 5
    // Then, summate and negate, then return sum:
    // while 1 do
    //  p *= xNegativePlusOne
    //  sum += p / k
    //  if k == inf then -- as k approaches infinity.
    //    exit -- break;
    //  end if
    //  k += 1
    // end while
    // return - (sum)
    //
    double x, p, sum, last, f;
    if (a <= 0.0)
    {
      exit(1);
    }
    else if (a >= 2.0)
    {
      // factor, first find n, then calculate e^n
      // log(m/e^n) + n = log(m)
      if (a > ACONST_E)
      {
        f = floor(log(a));
      }
      else
      {
        f = 1.0;
      }
      a = Divl((long double)a, (long double)power(ACONST_E, f));
    }
    else
    {
      f = 0.0;
    }
    x = 1 - a;
    last = 0.0;
    p = 1.0;
    sum = 0.0;
    for (int k = 1; k <= 1000000000; k++)
    {
      p *= x;
      sum += Divl((long double)p, (long double)k);
      if (sum == last)
      {
        break;
      }
      last = sum;
    }
    sum = f - (sum);
    return sum;
}

double Power(double base, double raisedTo)
{
    // b^x = e^(x * ln(b));
    double r;
    r = Exp(Log(base) * raisedTo);
    return r;
}

double GeneralRoot(double rooted, double nyNumber)
{
    double r;
    r = Power(rooted, MultInvl((long double)anyNumber));
    return r;
}

// Trig functions

// cos
// sin
// tan
// arctan
// arccos
// arcsin
// ACONST_PI

double Cos(double a)
{
// cos(x) = 1 - ((x^2)/(2!)) + ((x^4)/(4!)) - ((x^6)/(6!)) + ((x^8)/(8!)) - ...
  double f, r, d, c, x;
  f = 1.0;
  r = 1.0;
  d = -1.0;
  c = 0.0;
  x = a;
  for (int i = 2; i <= 1000000000; i++)
  {
    f *= i;
    x *= a;
    r += (d * Divl((long double)x, (long double)f));
    if (r == c)
    {
      break;
    }
    c = r;
    i++; // i, increment by 2, each loop
    f *= i; // (i + 1)
    x *= a;
    d *= (-1.0);
  }
  return r;
}

double Sin(double a)
{
// sine(x) = x - ((x^3)/(3!)) + ((x^5)/(5!)) - ((x^7)/(7!)) + ((x^9)/(9!)) - ...
  double f, r, d, c, x;
  f = 2.0;
  r = a;
  d = -1.0;
  c = 0.0;
  x = a * a;
  for (int i = 3; i <= 1000000001; i++)
  {
    f *= i;
    x *= a;
    r += (d * Divl((long double)x, (long double)f));
    if (r == c)
    {
      break;
    }
    c = r;
    i++; // by 2
    f *= i; // (i + 1)
    x *= a;
    d *= (-1.0);
  }
  return r;
}


double Tan(double a)
{
  double r;
  r = Divl((long double)Sin(a), (long double)Cos(a));
  return r;
}


// arc functions

double ArcTan(double a)
{
// Begin ArcTanExpA()
//                z        +inf         n             2kz^2
// arctan(z) = ------- * Sumation of Product of --------------------
//             1 + z^2     n=0         k=1      (2k + 1) * (1 + z^2)
//
// (The term in the sum for n = 0 is the empty product, so is 1.)
//
// ans = (z / (1 + z*z)) * Sum [ n=0 to +inf ] Prod [ k=1 to n ] (2kz*z) / ((2k + 1) * (1 + z*z))
//
  double b, r, s, p, c;
  b = a * a + 1.0;
  s = 1.0;
  for (int n = 1; n <= 1000000000; n++)
  {
    p = 1.0;
    for (int k = 1; k <= n; k++)
    {
      p *= Divl((long double)(k * 2.0 * (b - 1.0)), (long double)((k * 2.0 + 1.0) * b));
    }
    c = s;
    s += p;
    if (s = c)
    {
      break;
    }
  }
  r = Divl((long double)a, (long double)b) * s;
  return r;
}

const double ACONST_PI = 4.0 * ArcTan(1.0); // 4 * arctan(1)

double ArcTan2(double y, double x)
{
// ArcTan2(y, x) = ArcTan(y/x)
//
// The comments below use UTF-8
//
// atan2(y,x) = arctan(y/x) if x > 0,
// atan2(y,x) = arctan(y/x) + π if x < 0 and y≥0,
// atan2(y,x) = arctan(y/x) - π if x < 0 and y < 0,
// atan2(y,x) = +π / 2 if x=0 and y > 0,
// atan2(y,x) = - π / 2 if x=0 and y < 0,
// atan2(y,x) = undefined if x=0 and y=0
//
// Table:
//      x     arctan(x) (°)  arctan(x) (rad.)
//      -∞    -90°   -π/2
//      -√3   -60°   -π/3
//      -1    -45°   -π/4
//      -1/√3 -30°   -π/6
//      0       0°     0
//      1/√3  30°    π/6
//      1     45°    π/4
//      √3    60°    π/3
//      +∞    90°    π/2
//
    double tmp;
    if (x == 0.0) // x == 0
    {
        if (y == 0.0) // y == 0
        {
            exit(1);
        }
        tmp = Divl((long double)ACONST_PI, (long double)2.0); // half of PI.
        if (y > 0.0) // y > 0
        {
            return tmp;
        }
        else if (y < 0.0) // y < 0
        {
            return -(tmp); // negated tmp
        }
    }
    tmp = ArcTan(Divl((long double)y, (long double)x));
    if (x > 0.0) // x > 0
    {
        return tmp;
    }
    if (y < 0.0) // y < 0
    {
        return adjust((long double)(tmp - ACONST_PI));
    }
    else
    {
        return adjust((long double)(tmp + ACONST_PI));
    }
}

double ArcSin(double a)
{
// arcsin(x) = arctan( x / sqrt(1 - x^2) )
  double r;
  r = Sqrt(1.0 - a * a);
  r = Divl((long double)a, (long double)r);
  r = ArcTan(r);
  return r;
}
  
double ArcCos(double a)
{
// arccos(x) = arctan( sqrt(1 - x^2) / x )
// Limited domain: -1 to 1
// Also:
//   arccos(x) = arcsin(1) - arcsin(x)
//   arccos(x) = (EunPi / 2) - arcsin(x)
  double r;
  r = Sqrt(1.0 - a * a);
  r = Divl((long double)r, (long double)a);
  r = ArcTan(r);
  return r;
}

// other trig functions

// Cosh
// Sinh
// Tanh
// ArcCosh
// ArcCot
// ArcCoth
// ArcCsc
// ArcCsch
// ArcSec
// ArcSech
// ArcSinh
// ArcTanh
// Cot
// Coth
// Csc
// Csch
// Sec
// Sech


double Cosh(double a)
{
// cosh(x) = (e^(x) + e^(-x)) / 2
    double r;
    r = Divl((long double)(Exp(a) + Exp(-a)), (long double)2.0);
    return r;
}

double Sinh(double a)
{
// sinh(x) = (e^(x) - e^(-x)) / 2
    double r;
    r = Divl((long double)(Exp(a) - Exp(-a)), (long double)2.0);
    return r;
}

double Tanh(double a)
{
// tanh(x) = e^(2*x) => a; (a - 1) / (a + 1)
    double r;
    r = Exp(a * 2.0);
    r = Divl((long double)(r - 1.0), (long double)(r + 1.0));
    return r;
}

double ArcCosh(double a)
{
// arccosh(x) = x >= 1; ln(x + sqrt(x^2 - 1))
  double r;
  if (a < 1.0)
  {
    exit(1);
  }
  r = Sqrt(a * a - 1.0);
  r = Log(a + r);
  return r;
}

double ArcCot(double a)
{
    double r;
    if (a == 0.0)
    {
        return Divl((long double)ACONST_PI, (long double)2.0);
    }
    r = ArcTan(MultInvl((long double)a));
    if (a < 0.0)
    {
        r += ACONST_PI;
        r = adjust(r);
    }
    return r;
}

double ArcCoth(double a)
{
// arccoth(x) = abs(x) > 1; ln((x + 1)/(x - 1)) / 2
    double r;
    if (abs(a) <= 1.0)
    {
        exit(1);
    }
    r = Divl((long double)(a + 1.0), (long double)(a - 1.0));
    r = Divl((long double)LnAtom(r), (long double)2.0);
    return r;
}

/*

global function ArcCscAtom(atom a)
  atom r
  r = ArcSinAtom(MultInvAtom(a))
  return r
end function

global function ArcCsc(object x)
  sequence s
  if atom(x) then
    return ArcCscAtom(x)
  end if
  s = repeat(0, length(x))
  for i = 1 to length(x) do
    s[i] = ArcCsc(x[i])
  end for
  return s
end function

global function ArcCschAtom(atom a)
-- arccsch(x) = x != 0; 1 / x => a; ln(a + sqrt(a^2 + 1))
    atom r
    if a = 0 then
        abort(1)
    end if
    r = MultInvAtom(a)
    r = LnAtom(r + SqrtAtom(r * r + 1))
    return r
end function

global function ArcCsch(object x)
  sequence s
  if atom(x) then
    return ArcCschAtom(x)
  end if
  s = repeat(0, length(x))
  for i = 1 to length(x) do
    s[i] = ArcCsch(x[i])
  end for
  return s
end function

global function ArcSecAtom(atom a)
    atom r
    r = ArcCosAtom(MultInvAtom(a))
    return r
end function

global function ArcSec(object x)
  sequence s
  if atom(x) then
    return ArcSecAtom(x)
  end if
  s = repeat(0, length(x))
  for i = 1 to length(x) do
    s[i] = ArcSec(x[i])
  end for
  return s
end function

global function ArcSechAtom(atom a)
-- arcsech(x) = 0 < x <= 1; 1 / x => a; ln(a + sqrt(a^2 - 1)) :: ln((1 + sqrt(1 - x^2)) / x)
    atom r
    if a <= 0 or a > 1 then
        abort(1)
    end if
    r = MultInvAtom(a)
    r = LnAtom(r + SqrtAtom(r * r - 1))
    return r
end function

global function ArcSech(object x)
  sequence s
  if atom(x) then
    return ArcSechAtom(x)
  end if
  s = repeat(0, length(x))
  for i = 1 to length(x) do
    s[i] = ArcSech(x[i])
  end for
  return s
end function

global function ArcSinhAtom(atom a)
-- arcsinh(x) = ln(x + sqrt(x^2 + 1))
    atom r
    r = SqrtAtom(a * a + 1)
    r = LnAtom(a + r)
    return r
end function

global function ArcSinh(object x)
  sequence s
  if atom(x) then
    return ArcSinhAtom(x)
  end if
  s = repeat(0, length(x))
  for i = 1 to length(x) do
    s[i] = ArcSinh(x[i])
  end for
  return s
end function

global function ArcTanhAtom(atom a)
-- arctanh(x) = abs(x) < 1; ln((1 + x)/(1 - x)) / 2
    atom r
    if abs(a) >= 1 then
        abort(1)
    end if
    r = DivAtom(a + 1, 1 - a)
    r = LnAtom(r) / 2
    return adjust_atom(r)
end function

global function ArcTanh(object x)
  sequence s
  if atom(x) then
    return ArTanhAtom(x)
  end if
  s = repeat(0, length(x))
  for i = 1 to length(x) do
    s[i] = ArcTanh(x[i])
  end for
  return s
end function

global function CotAtom(atom a)
  atom r
  r = MultInvAtom(TanAtom(a))
  return r
end function

global function Cot(object x)
  sequence s
  if atom(x) then
    return CotAtom(x)
  end if
  s = repeat(0, length(x))
  for i = 1 to length(x) do
    s[i] = Cot(x[i])
  end for
  return s
end function

global function CothAtom(atom a)
-- coth(x) = x != 0; 1 / tanh(x)
    atom r
    if a = 0 then
        abort(1)
    end if
    r = MultInvAtom(TanhAtom(a))
    return r
end function

global function Coth(object x)
  sequence s
  if atom(x) then
    return CothAtom(x)
  end if
  s = repeat(0, length(x))
  for i = 1 to length(x) do
    s[i] = Coth(x[i])
  end for
  return s
end function

global function CscAtom(atom a)
    atom r
    r = MultInvAtom(SinAtom(a))
    return r
end function

global function Csc(object x)
  sequence s
  if atom(x) then
    return CscAtom(x)
  end if
  s = repeat(0, length(x))
  for i = 1 to length(x) do
    s[i] = Csc(x[i])
  end for
  return s
end function

global function CschAtom(atom a)
-- csch(x) = x != 0; 1 / sinh(x)
    atom r
    if a = 0 then
        abort(1)
    end if
    r = MultInvAtom(SinhAtom(a))
    return r
end function

global function Csch(object x)
  sequence s
  if atom(x) then
    return CschAtom(x)
  end if
  s = repeat(0, length(x))
  for i = 1 to length(x) do
    s[i] = Csch(x[i])
  end for
  return s
end function

global function SecAtom(atom a)
    atom r
    r = MultInvAtom(CosAtom(a))
    return r
end function

global function Sec(object x)
  sequence s
  if atom(x) then
    return SecAtom(x)
  end if
  s = repeat(0, length(x))
  for i = 1 to length(x) do
    s[i] = Sec(x[i])
  end for
  return s
end function

global function SechAtom(atom a)
-- sech(x) = 1 / cosh(x)
    atom r
    r = MultInvAtom(CoshAtom(a))
    return r
end function

global function Sech(object x)
  sequence s
  if atom(x) then
    return SechAtom(x)
  end if
  s = repeat(0, length(x))
  for i = 1 to length(x) do
    s[i] = Sech(x[i])
  end for
  return s
end function

global function RadiansToDegreesAtom(atom r)
    atom d
    d = DivAtom(r, (ACONST_PI / 2)) * 90
    return adjust_atom(d)
end function

global function RadiansToDegrees(object x)
  sequence s
  if atom(x) then
    return RadiansToDegreesAtom(x)
  end if
  s = repeat(0, length(x))
  for i = 1 to length(x) do
    s[i] = RadiansToDegrees(x[i])
  end for
  return s
end function

global function DegreesToRadiansAtom(atom d)
    atom r
    r = DivAtom(d, 90) * (ACONST_PI / 2)
    return adjust_atom(r)
end function

global function DegreesToRadians(object x)
  sequence s
  if atom(x) then
    return DegreesToRadiansAtom(x)
  end if
  s = repeat(0, length(x))
  for i = 1 to length(x) do
    s[i] = DegreesToRadians(x[i])
  end for
  return s
end function

-- end of file.
*/
