y// Copyright (c) 2025 James Cook
// bmath.h
// For the header.
// Description: Accurate math routines, using Newton's method.
// Author: James Cook (jmsck55 AT gmail DOT com)

// Routines for doubles using long doubles as intermediaries.
#pragma once

// Uses C99 version of math.h
//#include <math.h>

// constant accuracy = power(2, 52) -- question

#define BMATH_ITERS 1000000

#define RoundDownIf 3
#define RoundUpIf 4
#define RoundToZero 0
#define RoundNot 15

static int adjustMethod;

double adjust(long double g);

double MultInvl(long double x);

double Divl(long double n, long double d);

double Div(double n, double d);

double Round(double a, double precision = 1);


double NthRoot(double x, unsigned int n);


double Sqrt(double x);

double Cbrt(double x);


double Exp(double x);


// Raw function: Natural Logarithm

const double ACONST_E;

double Log(double a);

double Power(double base, double raisedTo);

double GeneralRoot(double rooted, double nyNumber);

// Trig functions

// cos
// sin
// tan
// arctan
// arccos
// arcsin
// ACONST_PI

double Cos(double a);

double Sin(double a);


double Tan(double a);


// arc functions

double ArcTan(double a);

const double ACONST_PI;

double ArcTan2(double y, double x);

double ArcSin(double a);
  
double ArcCos(double a);

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


double Cosh(double a);

double Sinh(double a);

double Tanh(double a);

double ArcCosh(double a);

double ArcCot(double a);

double ArcCoth(double a);

double ArcCsc(double a);

double ArcCsch(double a);

double ArcSec(double a);

double ArcSech(double a);

double ArcSinh(double a);

double ArcTanh(double a);

double Cot(double a);

double Coth(double a);

double Csc(double a);

double Csch(double a);

double Sec(double a);

double Sech(double a);

double RadiansToDegrees(double r);

double DegreesToRadians(double d);

// end of file.
