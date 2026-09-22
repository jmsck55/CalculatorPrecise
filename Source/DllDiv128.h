// Copyright (c) 2026 James Cook
// DllDiv128.h

// Notes: Do a "complex variable" version of this library later.

// Compile with:
// $>icx -I/usr/lib/gcc/x86_64-redhat-linux/14/include -lquadmath test_quadmath_snprintf.c

#pragma once

#ifndef __GNUC__
#define __float128 long double
#endif

__float128 Absq(__float128 a);

#define WAIT_FOR 10
#define BMATH_ITERS 1000000000

#define RoundDownIf 3
#define RoundUpIf 4
#define RoundToZero 0
#define RoundNot 15

static int adjustMethod = 0;

__float128 adjustq(__float128 g);

__float128 MultInvq(__float128 x);

__float128 Divq(__float128 n, __float128 d);

__float128 Roundq(__float128 a, __float128 precision);


__float128 NthRootq(__float128 x, unsigned int n);


__float128 Sqrtq(__float128 x);

__float128 Cbrtq(__float128 x);


__float128 Expq(__float128 x);


// Raw function: Natural Logarithm

void init_log();

__float128 GetE();

__float128 Logq(__float128 a);

__float128 Powerq(__float128 base, __float128 raisedTo);

__float128 GeneralRootq(__float128 rooted, __float128 anyNumber);


// Trig functions

// cos
// sin
// tan
// arctan
// arccos
// arcsin
// ACONST_PI

void init_trig();

__float128 GetPI();

__float128 Cosq(__float128 a);

__float128 Sinq(__float128 a);


__float128 Tanq(__float128 a);


// arc functions

__float128 ArcTanq(__float128 a);

__float128 ArcTan2q(__float128 y, __float128 x);

__float128 ArcSinq(__float128 a);

__float128 ArcCosq(__float128 a);

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


__float128 Coshq(__float128 a);

__float128 Sinhq(__float128 a);

__float128 Tanhq(__float128 a);

__float128 ArcCoshq(__float128 a);

__float128 ArcCotq(__float128 a);

__float128 ArcCothq(__float128 a);

__float128 ArcCscq(__float128 a);

__float128 ArcCschq(__float128 a);

__float128 ArcSecq(__float128 a);

__float128 ArcSechq(__float128 a);

__float128 ArcSinhq(__float128 a);

__float128 ArcTanhq(__float128 a);

__float128 Cotq(__float128 a);

__float128 Cothq(__float128 a);

__float128 Cscq(__float128 a);

__float128 Cschq(__float128 a);

__float128 Secq(__float128 a);

__float128 Sechq(__float128 a);

__float128 RadiansToDegreesq(__float128 r);

__float128 DegreesToRadiansq(__float128 d);

// end of file.
