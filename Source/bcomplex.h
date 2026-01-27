// Copyright (c) 2025 James Cook
// bcomplex.h


#pragma once

#include <complex.h>

double complex CMultInv(double complex n2);

double complex CDiv(double complex n1, double complex n2);

double CAbsoluteValue(double complex z);

double complex CSqrt(double complex z);

double complex CExp(double complex z);

double complex CLog(double complex z);

double complex CPower(double complex z, double complex raisedTo);

double complex CSin(double complex z);

double complex CCos(double complex z);

double complex CTan(double complex z);

double complex CCot(double complex z);

double complex CSinh(double complex z);

double complex CCosh(double complex z);

double complex CTanh(double complex z);

double complex CCoth(double complex z);

double complex CArcTan(double complex z);

double complex * CQuadraticEquation(double complex a, double complex b, double complex c);

// end of file.
