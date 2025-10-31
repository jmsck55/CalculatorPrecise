// Copyright (c) 2025 James Cook
// bcomplex.c, accuracate complex calculations
// uses Newton's method.

#include "bmath.h"
#include <complex.h>


double complex CMultInv(double complex n2)
{
// Reciprocal
    // Eun a, b, c
    // (a+bi)(a-bi) <=> a*a + b*b
    // n2 = (a+bi)
    // a = n2[1]
    // b = n2[2]
    // 1 / n2 <=> (a-bi) / (a*a + b*b)
    // <=> (a / (a*a + b*b)) - (b / (a*a + b*b))i
    // c = (a*a + b*b)
    // <=> (a / c) - (b / c)i
    double a, b, c, real, imag;
    double complex n1;
    a = creal(n2);
    b = cimag(n2);
    c = MultInvl((long double)((a * a) + (b * b)));
    real = (a * c);
    imag = -(b * c);
    n1 = real + imag * I;
    return n1;
}

double complex CDiv(double complex n1, double complex n2)
{
// correct!
    return n1 * CMultInv(n2);
}

double CAbsoluteValue(double complex z)
{
// same as: ComplexModulus or ComplexMagnitude

// Abs[z] == Sqrt[z Conjugate[z]]

// abs(z) = sqrt(z * conj(z))
// abs(x + iy) = sqrt((x + iy)(x - iy))

// (x + iy)(x - iy) = x*x -x*iy +x*iy -iy*iy
// x^2 - (iy)^2
// x^2 + y^2

// abs(x + iy) = sqrt(x^2 + y^2)

    double a, b, r;
    a = creal(z);
    b = cimag(z);
    r = Sqrt((a * a) + (b * b));
    return r;
}

double complex CSqrt(double complex z)
{
// sqrt(x + iy) <=> (1/2) * sqrt(2) * [ sqrt( sqrt(x*x + y*y) + x )  +  i*sign(y) * sqrt( sqrt(x*x + y*y) - x ) ]
  double x, y, a, b, c;
  x = creal(z);
  y = cimag(z);
  c = Sqrt((x * x) + (y * y));
  a = Sqrt(c + x);
  b = Sqrt(c - x);
  if (y < 0)
  {
    b = -(b);
  }
  z = (0.5 * Sqrt(2.0) * a) + (b * I);
  return z;
}

double complex CExp(double complex z)
{
// ComplexExponent()
    double complex r;
    double a, b;
    a = creal(z);
    b = cimag(z);
    a = Exp(a);
    r = (a * Cos(b)) + (a * Sin(b) * I);
    return r;
}

double complex CLog(double complex z)
{
// Natural Logarithm
// ln(z) = (ln(x^2 + y^2)/2) + arctan(y/x)i
    double complex r;
    double a, b;
    a = creal(z);
    b = cimag(z);
    r = (0.5 * Log((a * a) + (b * b))) + (ArcTan2(b, a) * I);
    return r;
}

double complex CPower(double complex z, double complex raisedTo)
{
    double complex r;
    r = CExp(raisedTo * CLog(z));
    return r;
}

// For z = x + iy,

// sin(z) = sin(x) * cosh(y) + i * cos(x) * sinh(y)
// {\displaystyle \sin {z}=\sin {x}\cosh {y}+i\cos {x}\sinh {y}}

double complex CSin(double complex z)
{
// sin(z) = (sin(x) * cosh(y)) + (cos(x) * sinh(y))i
    double complex r;
    double a, b;
    a = creal(z);
    b = cimag(z);
    r = (Sin(a) * Cosh(b)) + (Cos(a) * Sinh(b) * I);
    return r;
}

// cos(z) = cos(x) * cosh(y) − i * sin(x) * sinh(y)
// {\displaystyle \cos {z}=\cos {x}\cosh {y}-i\sin {x}\sinh {y}}

double complex CCos(double complex z)
{
// cos(z) = (cos(x) * cosh(y)) - (sin(x) * sinh(y))i
    double complex r;
    double a, b;
    a = creal(z);
    b = cimag(z);
    r = (Cos(a) * Cosh(b)) - (Sin(a) * Sinh(b) * I);
    return r;
}

// tan(z) = (tan(x) + i * tanh(y)) / (1 − i * tan(x) * tanh(y))
// {\displaystyle \tan {z}={\frac {\tan {x}+i\tanh {y}}{1-i\tan {x}\tanh {y}}}}

double complex CTan(double complex z)
{
// z = Real(x) + Imaginary(y)
// f(x,y) = cos(2x) + cosh(2y)
// tan(z) = (sin(2x)/f(x,y)) + (sinh(2y)/f(x,y))i
  double complex r;
  double x, y, f;
  x = creal(z) * 2.0;
  y = cimag(z) * 2.0;
  f = Cos(x) + Cosh(y);
  r = Div(Sin(x), f) + (Div(Sinh(y), f) * I);
  return r;
}

// cot(z) = −((1 + i * cot(x) * coth(y)) / (cot(x) − i * coth(y)))
// {\displaystyle \cot {z}=-{\frac {1+i\cot {x}\coth {y}}{\cot {x}-i\coth {y}}}}

double complex CCot(double complex z)
{
    // cot(z) = −((1 + i * cot(x) * coth(y)) / (cot(x) − i * coth(y)))
    double complex r;
    double a, b;
    a = creal(z);
    b = cimag(z);
    a = Cot(a);
    b = Coth(b);
    r = - CDiv(1.0 + I * a * b, a - I * b);
    return r;
}

// sinh(z) = sinh(x) * cos(y) + i * cosh(x) * sin(y)
// {\displaystyle \sinh {z}=\sinh {x}\cos {y}+i\cosh {x}\sin {y}}

double complex CSinh(double complex z)
{
// Sinus hyperbolic (sinh)
// sinh(z) = (sinh(x) * cos(y)) - (cosh(x) * sin(y))i
    double complex r;
    double a, b;
    a = creal(z);
    b = cimag(z);
    r = (Sinh(a) * Cos(b)) - (Cosh(a) * Sin(b) * I); // ? plus or minus?
    return r;
}

// cosh(z) = cosh(x) * cos(y) + i * sinh(x) * sin(y)
// {\displaystyle \cosh {z}=\cosh {x}\cos {y}+i\sinh {x}\sin {y}}

double complex CCosh(double complex z)
{
// Cosine hyperbolic (cosh)
// cosh(z) = (cosh(x) * cos(y)) - (sinh(x) * sin(y))i
    double complex r;
    double a, b;
    a = creal(z);
    b = cimag(z);
    r = (Cosh(a) * Cos(b)) - (Sinh(a) * Sin(b) * I); // plus or minus?
    return r;
}

// tanh(z) = (tanh(x) + i * tan(y)) / (1 + i * tanh(x) * tan(y))
// {\displaystyle \tanh {z}={\frac {\tanh {x}+i\tan {y}}{1+i\tanh {x}\tan {y}}}}

double complex CTanh(double complex z)
{
    // tanh(z) = (tanh(x) + i * tan(y)) / (1 + i * tanh(x) * tan(y))
    double complex r;
    double a, b;
    a = creal(z);
    b = cimag(z);
    a = Tanh(a);
    b = Tan(b);
    r = CDiv(a + I * b, 1.0 + I * a * b);
    return r;
}

// coth(z) = (1 − i * coth(x) * cot(y)) / (coth(x) − i * cot(y))
// {\displaystyle \coth {z}={\frac {1-i\coth {x}\cot {y}}{\coth {x}-i\cot {y}}}}

double complex CCoth(double complex z)
{
    // coth(z) = (1 − i * coth(x) * cot(y)) / (coth(x) − i * cot(y))
    double complex r;
    double a, b;
    a = creal(z);
    b = cimag(z);
    a = Coth(a);
    b = Cot(b);
    r = CDiv(1.0 - I * a * b, a - I * b);
    return r;
}

double complex CArcTan(double complex z)
{
// Given: arctan(x + iy), z = x + iy
// (1/2) * i * log(1 - i(x + iy)) - (1/2) * i * log(1 + i(x + iy))
// (1/2) * i * (log(-ix + y + 1) - log(ix - y + 1))
// (1/2) * i * (log(1 - i * z) - log(1 + i * z))
  double complex r;
  r = (0.5 * I) * (CLog(1.0 - (I * z)) - CLog(1.0 + (I * z)));
  return r;
}

double complex * CQuadraticEquation(double complex a, double complex b, double complex c)
{
    // The quadratic equation produces two answers (the answers may be the same)
    // ax^2 + bx + c
    // f(a,b,c) = (-b +-sqrt(b*b - 4*a*c)) / (2*a)
    // answer[0] = ((-b + sqrt(b*b - 4*a*c)) / (2*a))
    // answer[1] = ((-b - sqrt(b*b - 4*a*c)) / (2*a))
    double complex r[2];
    double complex s, t;
    t = (b * b) - (4.0 * a * c);
    s = CSqrt(t);
    t = 2.0 * a;
    r[0] = CDiv((-b) + s, t);
    r[1] = CDiv((-b) - s, t);
    return r;
}

// end of file.
