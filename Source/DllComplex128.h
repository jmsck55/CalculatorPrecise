// Copyright (c) 2026 James Cook
// 
// DllComplex128.h
//
// Best to use gcc compiler for this project.  MSVC does not support __complex128 type.
//

#pragma once

#ifndef __GNUC__
#include <complex.h>
#define __complex128 _Lcomplex
#define CAssign(complex, real, imag) complex = _LCOMPLEX_(real, imag)
#else
#include "quadmath.h"
#define CAssign(complex, real, imag) complex = real + imag * I
#endif

__complex128 CMultInvq(__complex128 n2);

__complex128 CDivq(__complex128 n1, __complex128 n2);

__float128 CAbsoluteValueq(__complex128 z);

__complex128 CSqrtq(__complex128 z);

__complex128 CExpq(__complex128 z);

__complex128 CLogq(__complex128 z);

__complex128 CPowerq(__complex128 z, __complex128 raisedTo);

__complex128 CSinq(__complex128 z);

__complex128 CCosq(__complex128 z);

__complex128 CTanq(__complex128 z);

__complex128 CCotq(__complex128 z);

__complex128 CSinhq(__complex128 z);

__complex128 CCoshq(__complex128 z);

__complex128 CTanhq(__complex128 z);

__complex128 CCothq(__complex128 z);

__complex128 CArcTanq(__complex128 z);

__complex128 * CQuadraticEquationq(__complex128 a, __complex128 b, __complex128 c);

// end of file.
