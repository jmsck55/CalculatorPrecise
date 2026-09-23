// Copyright (c) 2026 James Cook
// 
// DllComplex128.c, accuracate complex calculations
// uses Newton's method.

#include "DllDiv128.h"
#include "DllComplex128.h"


__complex128 CMultInvq(__complex128 n2)
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
    __float128 a, b, c, real, imag;
    __complex128 n1;
    a = crealq(n2);
    b = cimagq(n2);
    c = MultInvq((__float128)((a * a) + (b * b)));
    real = (a * c);
    imag = -(b * c);
    CAssign(n1, real, imag);
    return n1;
}

__complex128 CDivq(__complex128 n1, __complex128 n2)
{
// correct!
    __complex128 reciprocal = CMultInvq(n2);
    __float128 a = crealq(n1);
    __float128 b = cimagq(n1);
    __float128 c = crealq(reciprocal);
    __float128 d = cimagq(reciprocal);
    __complex128 result;

    CAssign(result, (a * c) - (b * d), (a * d) + (b * c));
    return result;
}

__float128 CAbsoluteValueq(__complex128 z)
{
// same as: ComplexModulus or ComplexMagnitude

// Abs[z] == Sqrt[z Conjugate[z]]

// abs(z) = sqrt(z * conj(z))
// abs(x + iy) = sqrt((x + iy)(x - iy))

// (x + iy)(x - iy) = x*x -x*iy +x*iy -iy*iy
// x^2 - (iy)^2
// x^2 + y^2

// abs(x + iy) = sqrt(x^2 + y^2)

    __float128 a, b, r;
    a = crealq(z);
    b = cimagq(z);
    r = Sqrtq((a * a) + (b * b));
    return r;
}

__complex128 CSqrtq(__complex128 z)
{
// sqrt(x + iy) <=> (1/2) * sqrt(2) * [ sqrt( sqrt(x*x + y*y) + x )  +  i*sign(y) * sqrt( sqrt(x*x + y*y) - x ) ]
  __float128 x, y, a, b, c;
  x = crealq(z);
  y = cimagq(z);
  c = Sqrtq((x * x) + (y * y));
  a = Sqrtq(c + x);
  b = Sqrtq(c - x);
  if (y < 0)
  {
    b = -(b);
  }
  CAssign(z, 0.5 * Sqrtq((__float128)2.0) * a, b);
  return z;
}

__complex128 CExpq(__complex128 z)
{
// ComplexExponent()
    __complex128 r;
    __float128 a, b;
    a = crealq(z);
    b = cimagq(z);
    a = Expq(a);
    CAssign(r, a * Cosq(b), a * Sinq(b));
    return r;
}

__complex128 CLogq(__complex128 z)
{
// Natural Logarithm
// ln(z) = (ln(x^2 + y^2)/2) + arctan(y/x)i
    __complex128 r;
    __float128 a, b, real, imag;
    a = crealq(z);
    b = cimagq(z);
    real = (__float128)0.5 * Logq((a * a) + (b * b));
    imag = ArcTan2q(b, a);
    CAssign(r, real, imag);
    return r;
}

__complex128 CPowerq(__complex128 z, __complex128 raisedTo)
{
    __complex128 r;
    __complex128 logarithm;
    __complex128 product;
    __float128 real, imag;

    logarithm = CLogq(z);
    real = (crealq(raisedTo) * crealq(logarithm)) -
           (cimagq(raisedTo) * cimagq(logarithm));
    imag = (crealq(raisedTo) * cimagq(logarithm)) +
           (cimagq(raisedTo) * crealq(logarithm));
    CAssign(product, real, imag);
    r = CExpq(product);
    return r;
}

// For z = x + iy,

// sin(z) = sin(x) * cosh(y) + i * cos(x) * sinh(y)
// {\displaystyle \sin {z}=\sin {x}\cosh {y}+i\cos {x}\sinh {y}}

__complex128 CSinq(__complex128 z)
{
// sin(z) = (sin(x) * cosh(y)) + (cos(x) * sinh(y))i
    __complex128 r;
    __float128 a, b;
    a = crealq(z);
    b = cimagq(z);
    CAssign(r, Sinq(a) * Coshq(b), Cosq(a) * Sinhq(b));
    return r;
}

// cos(z) = cos(x) * cosh(y) − i * sin(x) * sinh(y)
// {\displaystyle \cos {z}=\cos {x}\cosh {y}-i\sin {x}\sinh {y}}

__complex128 CCosq(__complex128 z)
{
// cos(z) = (cos(x) * cosh(y)) - (sin(x) * sinh(y))i
    __complex128 r;
    __float128 a, b;
    a = crealq(z);
    b = cimagq(z);
    CAssign(r, Cosq(a) * Coshq(b), -(Sinq(a) * Sinhq(b)));
    return r;
}

// tan(z) = (tan(x) + i * tanh(y)) / (1 − i * tan(x) * tanh(y))
// {\displaystyle \tan {z}={\frac {\tan {x}+i\tanh {y}}{1-i\tan {x}\tanh {y}}}}

__complex128 CTanq(__complex128 z)
{
// z = Real(x) + Imaginary(y)
// f(x,y) = cos(2x) + cosh(2y)
// tan(z) = (sin(2x)/f(x,y)) + (sinh(2y)/f(x,y))i
    __complex128 r;
    __float128 x, y, f;
    x = crealq(z) * 2.0;
    y = cimagq(z) * 2.0;
    f = Cosq(x) + Coshq(y);
    CAssign(r, Divq(Sinq(x), f), Divq(Sinhq(y), f));
    return r;
}

// cot(z) = −((1 + i * cot(x) * coth(y)) / (cot(x) − i * coth(y)))
// {\displaystyle \cot {z}=-{\frac {1+i\cot {x}\coth {y}}{\cot {x}-i\coth {y}}}}

__complex128 CCotq(__complex128 z)
{
    // cot(z) = −((1 + i * cot(x) * coth(y)) / (cot(x) − i * coth(y)))
    __complex128 numerator, denominator, quotient, r;
    __float128 a, b;
    a = crealq(z);
    b = cimagq(z);
    a = Cotq(a);
    b = Cothq(b);
    CAssign(numerator, 1.0, a * b);
    CAssign(denominator, a, -b);
    quotient = CDivq(numerator, denominator);
    CAssign(r, -crealq(quotient), -cimagq(quotient));
    return r;
}

// sinh(z) = sinh(x) * cos(y) + i * cosh(x) * sin(y)
// {\displaystyle \sinh {z}=\sinh {x}\cos {y}+i\cosh {x}\sin {y}}

__complex128 CSinhq(__complex128 z)
{
// Sinus hyperbolic (sinh)
// sinh(z) = (sinh(x) * cos(y)) - (cosh(x) * sin(y))i
    __complex128 r;
    __float128 a, b;
    a = crealq(z);
    b = cimagq(z);
    CAssign(r, Sinhq(a) * Cosq(b), Coshq(a) * Sinq(b)); // plus or minus?
    return r;
}

// cosh(z) = cosh(x) * cos(y) + i * sinh(x) * sin(y)
// {\displaystyle \cosh {z}=\cosh {x}\cos {y}+i\sinh {x}\sin {y}}

__complex128 CCoshq(__complex128 z)
{
// Cosine hyperbolic (cosh)
// cosh(z) = (cosh(x) * cos(y)) - (sinh(x) * sin(y))i
    __complex128 r;
    __float128 a, b;
    a = crealq(z);
    b = cimagq(z);
    CAssign(r, Coshq(a) * Cosq(b), Sinhq(a) * Sinq(b)); // plus or minus?
    return r;
}

// tanh(z) = (tanh(x) + i * tan(y)) / (1 + i * tanh(x) * tan(y))
// {\displaystyle \tanh {z}={\frac {\tanh {x}+i\tan {y}}{1+i\tanh {x}\tan {y}}}}

__complex128 CTanhq(__complex128 z)
{
    // tanh(z) = (tanh(x) + i * tan(y)) / (1 + i * tanh(x) * tan(y))
    __complex128 numerator, denominator, r;
    __float128 a, b;
    a = crealq(z);
    b = cimagq(z);
    a = Tanhq(a);
    b = Tanq(b);
    CAssign(numerator, a, b);
    CAssign(denominator, (__float128)1.0, a * b);
    r = CDivq(numerator, denominator);
    return r;
}

// coth(z) = (1 − i * coth(x) * cot(y)) / (coth(x) − i * cot(y))
// {\displaystyle \coth {z}={\frac {1-i\coth {x}\cot {y}}{\coth {x}-i\cot {y}}}}

__complex128 CCothq(__complex128 z)
{
    // coth(z) = (1 − i * coth(x) * cot(y)) / (coth(x) − i * cot(y))
    __complex128 r, numerator, denominator;
    __float128 a, b;
    a = crealq(z);
    b = cimagq(z);
    a = Cothq(a);
    b = Cotq(b);
    CAssign(numerator, (__float128)1.0, -(a * b));
    CAssign(denominator, a, -b);
    r = CDivq(numerator, denominator);
    return r;
}

__complex128 CArcTanq(__complex128 z)
{
// Given: arctan(x + iy), z = x + iy
// (1/2) * i * log(1 - i(x + iy)) - (1/2) * i * log(1 + i(x + iy))
// (1/2) * i * (log(-ix + y + 1) - log(ix - y + 1))
// (1/2) * i * (log(1 - i * z) - log(1 + i * z))
        __complex128 r, minusIz, plusIz, logMinus, logPlus;
        __float128 x, y, realDifference, imagDifference;

        // MSVC represents __complex128 as a struct, so do not use complex
        // arithmetic operators here.
        x = crealq(z);
        y = cimagq(z);

        // 1 - i(x + iy) = (1 + y) - ix
        CAssign(minusIz, (__float128)1.0 + y, -x);
        // 1 + i(x + iy) = (1 - y) + ix
        CAssign(plusIz, (__float128)1.0 - y, x);

        logMinus = CLogq(minusIz);
        logPlus = CLogq(plusIz);
        realDifference = crealq(logMinus) - crealq(logPlus);
        imagDifference = cimagq(logMinus) - cimagq(logPlus);

        // (i / 2) * (realDifference + i * imagDifference)
        CAssign(r, -((__float128)0.5 * imagDifference),
            (__float128)0.5 * realDifference);
        return r;
}

__complex128 * CQuadraticEquationq(__complex128 a, __complex128 b, __complex128 c)
{
    // The quadratic equation produces two answers (the answers may be the same)
    // ax^2 + bx + c
    // f(a,b,c) = (-b +-sqrt(b*b - 4*a*c)) / (2*a)
    // answer[0] = ((-b + sqrt(b*b - 4*a*c)) / (2*a))
    // answer[1] = ((-b - sqrt(b*b - 4*a*c)) / (2*a))
    __complex128 r[2];
    __complex128 s, t, bSquared, ac;
    __float128 ar, ai, br, bi, cr, ci;
    ar = crealq(a);
    ai = cimagq(a);
    br = crealq(b);
    bi = cimagq(b);
    cr = crealq(c);
    ci = cimagq(c);
    CAssign(bSquared, (br * br) - (bi * bi), (__float128)2.0 * br * bi);
    CAssign(ac, (__float128)4.0 * ((ar * cr) - (ai * ci)),
        (__float128)4.0 * ((ar * ci) + (ai * cr)));
    CAssign(t, crealq(bSquared) - crealq(ac), cimagq(bSquared) - cimagq(ac));
    s = CSqrtq(t);
    CAssign(t, (__float128)2.0 * ar, (__float128)2.0 * ai);
    CAssign(r[0], -br + crealq(s), -bi + cimagq(s));
    r[0] = CDivq(r[0], t);
    CAssign(r[1], -br - crealq(s), -bi - cimagq(s));
    r[1] = CDivq(r[1], t);
    return r;
}

// end of file.
