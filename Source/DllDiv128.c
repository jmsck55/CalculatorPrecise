// Copyright (c) 2026 James Cook
// DllDiv128.c


// Include directories:
// C:\strawberry-perl-5.32.1.1-64bit-portable\c\lib\gcc\x86_64-w64-mingw32\8.3.0\include
// $>icx -I/usr/lib/gcc/x86_64-redhat-linux/14/include -lquadmath test_quadmath_snprintf.c

#include "quadmath.h"
#include "DllDiv128.h"


__float128 Absq(__float128 a)
{
        if (a < (__float128)0.0)
        {
                a = -a;
        }
        return a;
}

//#define BMATH_ITERS 1000000

//#define RoundDownIf 3
//#define RoundUpIf 4
//#define RoundToZero 0

static int adjustMethod = RoundToZero;

__float128 adjustq(__float128 g)
{
        __float128 b = (__float128)g;
        __float128 c = (__float128)b;
        if (c == g)
        {
                return b;
        }
        if (
                (((adjustMethod == 3) || (adjustMethod == 0)) && (c > g))
                ) // round down. // if (fasbl(c) > fasbl(g)) // round to zero.
        {
                g *= 2;
                g -= c;
                b = (__float128)g;
        }
        else
        {
                if ((((adjustMethod == 4) || (adjustMethod == 0)) && (c < g))) // round up.
                {
                        c *= 2;
                        c -= g;
                        b = (__float128)c;
                }
        }
        return b;
}

__float128 MultInvq(__float128 x)
{
        // performs 1/x
        __float128 t, g;
        __float128 n, r, a, b, c;
        g = (__float128)1.0 / x;
        if ((__float128)x == (__float128)0.0)
        {
                return (__float128)g; // infinity
        }
        t = g * x;
        if (t == (__float128)1.0)
        {
                return adjustq(g);
                // return (__float128)g;
        }
        r = (__float128)0.0;
        a = (__float128)0.0;
        b = (__float128)0.0;
        c = (__float128)0.0;
        for (int i = 1; i <= BMATH_ITERS; i++)
        {
                g *= ((__float128)2.0 - t);
                if (r == (__float128)1.0) {
                        break;
                }
                t = g * x;
                n = adjustq(t);
                if (n == r) {
                        break;
                }
                if (i > WAIT_FOR) {
                        if (n == a) {
                                break;
                        }
                        if (n == b) {
                                break;
                        }
                        if (n == c) {
                                break;
                        }
                        c = b;
                        b = a;
                        a = r;
                }
        }
        return adjustq(g);
}

__float128 Divq(__float128 n, __float128 d)
{
        if (n == (__float128)0.0) {
                return 0.0;
        }
        if (d == (__float128)1.0) {
                return adjustq(n);
        }
        if (d == (__float128)2.0) {
                return (__float128)(n / (__float128)2.0);
        }
        n *= (__float128)MultInvq(d);
        return adjustq(n);
}

__float128 Roundq(__float128 a, __float128 precision)
{
        __float128 b;
        b = (__float128)precision * (__float128)a;
        b += (__float128)0.5;
        b = floorq(b);
        return Divq(b, (__float128)precision);
}


__float128 NthRootq(__float128 x, unsigned int n)
{
        __float128 quotient, average, guess;
        guess = powq((__float128)x, (__float128)MultInvq((__float128)n));
        quotient = Divq((__float128)x, powq((__float128)guess, (__float128)(n - 1)));
        average = Divq((__float128)(n - 1) * (__float128)guess + (__float128)quotient, (__float128)n);
        return average;
}


__float128 Sqrtq(__float128 x)
{
        return NthRootq(x, 2);
}

__float128 Cbrtq(__float128 x)
{
        return NthRootq(x, 3);
}


__float128 Expq(__float128 x)
{
        // using taylor series
        // https://en.wikipedia.org/wiki/TaylorSeries
        //
        //  exp(1) = sum of k=0 to inf (1/k!);
        //  exp(x) = sum of k=0 to inf ((x^k)/k!);
        //
        __float128 sum, num, den, last;
        num = (__float128)1.0;
        den = (__float128)1.0;
        sum = (__float128)1.0;
        last = (__float128)0.0;
        for (int i = 1; i <= 1000000000; i++)
        {
                num *= x;
                den *= i; // number of iterations
                sum += Divq((__float128)num, (__float128)den);
                if (sum == last)
                {
                        break;
                }
                last = sum;
        }
        return sum;
}


// Raw function: Natural Logarithm

static __float128 ACONST_E = (__float128)0.0;

void init_log() {
        ACONST_E = Expq((__float128)1.0);
}

__float128 GetE()
{
        if (ACONST_E == (__float128)0.0)
        {
                init_log();
        }
        return ACONST_E;
}


__float128 Logq(__float128 a)
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
        __float128 x, p, sum, last, f;
        if (ACONST_E == (__float128)0.0)
        {
                init_log();
        }
        if (a <= (__float128)0.0)
        {
                exit(1);
        }
        else if (a >= (__float128)2.0)
        {
                // factor, first find n, then calculate e^n
                // log(m/e^n) + n = log(m)
                if (a > ACONST_E)
                {
                        f = floorq(logq(a));
                }
                else
                {
                        f = (__float128)1.0;
                }
                a = Divq((__float128)a, (__float128)powq(ACONST_E, f));
        }
        else
        {
                f = (__float128)0.0;
        }
        x = 1 - a;
        last = (__float128)0.0;
        p = (__float128)1.0;
        sum = (__float128)0.0;
        for (int k = 1; k <= 1000000000; k++)
        {
                p *= x;
                sum += Divq((__float128)p, (__float128)k);
                if (sum == last)
                {
                        break;
                }
                last = sum;
        }
        sum = f - (sum);
        return sum;
}

__float128 Powerq(__float128 base, __float128 raisedTo)
{
        // b^x = e^(x * ln(b));
        __float128 r;
        r = Expq(Logq(base) * raisedTo);
        return r;
}

__float128 GeneralRootq(__float128 rooted, __float128 anyNumber)
{
        __float128 r;
        r = Powerq(rooted, MultInvq((__float128)anyNumber));
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

static __float128 ACONST_PI = (__float128)0.0;

void init_trig() {
        ACONST_PI = (__float128)4.0 * ArcTanq((__float128)1.0); // 4 * arctan(1)
}

__float128 GetPI()
{
        if (ACONST_PI == (__float128)0.0)
        {
                init_trig();
        }
        return ACONST_PI;
}


__float128 Cosq(__float128 a)
{
        // Range: -PI/2 to PI/2, exclusive
        // cos(x) = 1 - ((x^2)/(2!)) + ((x^4)/(4!)) - ((x^6)/(6!)) + ((x^8)/(8!)) - ...
        __float128 f, r, d, c, x;
        if (a == (__float128)0.0)
        {
                return (__float128)1.0;
        }
        if (ACONST_PI == (__float128)0.0)
        {
                init_trig();
        }
        x = Divq((__float128)ACONST_PI, (__float128)2.0);
        if (a < -(x))
        {
                exit(1);
        }
        if (a > x)
        {
                exit(1);
        }
        f = (__float128)1.0;
        r = (__float128)1.0;
        d = (__float128)-1.0;
        c = (__float128)0.0;
        x = a;
        for (int i = 2; i <= 1000000000; i++)
        {
                f *= i;
                x *= a;
                r += (d * Divq((__float128)x, (__float128)f));
                if (r == c)
                {
                        break;
                }
                c = r;
                i++; // i, increment by 2, each loop
                f *= i; // (i + 1)
                x *= a;
                d *= ((__float128)-1.0);
        }
        return r;
}

__float128 Sinq(__float128 a)
{
        // Cases: 0 equals zero (0)
        // Range: -PI/2 to PI/2, inclusive
        // To find sin(x):
        // sin(x) = cos(x - PI/2)
        // y = sin(x)
        // y = cos(x - PI/2)
        // sine(x) = x - ((x^3)/(3!)) + ((x^5)/(5!)) - ((x^7)/(7!)) + ((x^9)/(9!)) - ...
        __float128 f, r, d, c, x;
        if (a == (__float128)0.0)
        {
                return (__float128)0.0;
        }
        if (ACONST_PI == (__float128)0.0)
        {
                init_trig();
        }
        x = Divq((__float128)ACONST_PI, (__float128)2.0);
        if (a < -(x))
        {
                exit(1);
        }
        if (a > x)
        {
                exit(1);
        }
        f = (__float128)2.0;
        r = a;
        d = (__float128)-1.0;
        c = (__float128)0.0;
        x = a * a;
        for (int i = 3; i <= 1000000001; i++)
        {
                f *= i;
                x *= a;
                r += (d * Divq((__float128)x, (__float128)f));
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


__float128 Tanq(__float128 a)
{
        __float128 r;
        r = Divq((__float128)Sinq(a), (__float128)Cosq(a));
        return r;
}


// arc functions

__float128 ArcTanq(__float128 a)
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
        __float128 b, r, s, p, c;
        b = a * a + (__float128)1.0;
        s = (__float128)1.0;
        for (int n = 1; n <= 1000000000; n++)
        {
                p = (__float128)1.0;
                for (int k = 1; k <= n; k++)
                {
                        p *= Divq((__float128)(k * (__float128)2.0 * (b - (__float128)1.0)), (__float128)((k * (__float128)2.0 + (__float128)1.0) * b));
                }
                c = s;
                s += p;
                if (s = c)
                {
                        break;
                }
        }
        r = Divq((__float128)a, (__float128)b) * s;
        return r;
}

__float128 ArcTan2q(__float128 y, __float128 x)
{
        // ArcTan2q(y, x) = ArcTanq(y/x)
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
        __float128 tmp;
        if (ACONST_PI == (__float128)0.0)
        {
                init_trig();
        }
        if (x == (__float128)0.0) // x == 0
        {
                if (y == (__float128)0.0) // y == 0
                {
                        exit(1);
                }
                tmp = Divq((__float128)ACONST_PI, (__float128)2.0); // half of PI.
                if (y > (__float128)0.0) // y > 0
                {
                        return tmp;
                }
                else if (y < (__float128)0.0) // y < 0
                {
                        return -(tmp); // negated tmp
                }
        }
        tmp = ArcTanq(Divq((__float128)y, (__float128)x));
        if (x > (__float128)0.0) // x > 0
        {
                return tmp;
        }
        if (y < (__float128)0.0) // y < 0
        {
                return adjustq((__float128)(tmp - ACONST_PI));
        }
        else
        {
                return adjustq((__float128)(tmp + ACONST_PI));
        }
}

__float128 ArcSinq(__float128 a)
{
        // arcsin(x) = arctan( x / sqrt(1 - x^2) )
        __float128 r;
        if (a == (__float128)1.0)
        {
                exit(1);
        }
        if (a == -(__float128)-1.0) // -1.0
        {
                exit(1);
        }
        r = Sqrtq((__float128)1.0 - a * a);
        r = Divq((__float128)a, (__float128)r);
        r = ArcTanq(r);
        return r;
}

__float128 ArcCosq(__float128 a)
{
        // arccos(x) = arctan( sqrt(1 - x^2) / x )
        // Limited domain: -1 to 1
        // Also:
        //   arccos(x) = arcsin(1) - arcsin(x)
        //   arccos(x) = (EunPi / 2) - arcsin(x)
        __float128 r;
        if (a < (__float128)-1.0)
        {
                exit(1);
        }
        if (a > (__float128)1.0)
        {
                exit(1);
        }
        r = Sqrtq((__float128)1.0 - a * a);
        r = Divq((__float128)r, (__float128)a);
        r = ArcTanq(r);
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


__float128 Coshq(__float128 a)
{
        // cosh(x) = (e^(x) + e^(-x)) / 2
        __float128 r;
        r = Divq((__float128)(Expq(a) + Expq(-a)), (__float128)2.0);
        return r;
}

__float128 Sinhq(__float128 a)
{
        // sinh(x) = (e^(x) - e^(-x)) / 2
        __float128 r;
        r = Divq((__float128)(Expq(a) - Expq(-a)), (__float128)2.0);
        return r;
}

__float128 Tanhq(__float128 a)
{
        // tanh(x) = e^(2*x) => a; (a - 1) / (a + 1)
        __float128 r;
        r = Expq(a * (__float128)2.0);
        r = Divq((__float128)(r - (__float128)1.0), (__float128)(r + (__float128)1.0));
        return r;
}

__float128 ArcCoshq(__float128 a)
{
        // arccosh(x) = x >= 1; ln(x + sqrt(x^2 - 1))
        __float128 r;
        if (a < (__float128)1.0)
        {
                exit(1);
        }
        r = Sqrtq(a * a - (__float128)1.0);
        r = Logq(a + r);
        return r;
}

__float128 ArcCotq(__float128 a)
{
        __float128 r;
        if (ACONST_PI == (__float128)0.0)
        {
                init_trig();
        }
        if (a == (__float128)0.0)
        {
                return Divq((__float128)ACONST_PI, (__float128)2.0);
        }
        r = ArcTanq(MultInvq((__float128)a));
        if (a < (__float128)0.0)
        {
                r += ACONST_PI;
                r = adjustq(r);
        }
        return r;
}

__float128 ArcCothq(__float128 a)
{
        // arccoth(x) = abs(x) > 1; ln((x + 1)/(x - 1)) / 2
        __float128 r;
        if (Absq(a) <= (__float128)1.0)
        {
                exit(1);
        }
        r = Divq((__float128)(a + (__float128)1.0), (__float128)(a - (__float128)1.0));
        r = Divq((__float128)Logq(r), (__float128)2.0);
        return r;
}

__float128 ArcCscq(__float128 a)
{
        __float128 r;
        r = ArcSinq(MultInvq((__float128)a));
        return r;
}

__float128 ArcCschq(__float128 a)
{
        // arccsch(x) = x != 0; 1 / x => a; ln(a + sqrt(a^2 + 1))
        __float128 r;
        if (a == (__float128)0.0)
        {
                exit(1);
        }
        r = MultInvq((__float128)a);
        r = Logq(r + Sqrtq(r * r + (__float128)1.0));
        return r;
}

__float128 ArcSecq(__float128 a)
{
        __float128 r;
        r = ArcCosq(MultInvq((__float128)a));
        return r;
}

__float128 ArcSechq(__float128 a)
{
        // arcsech(x) = 0 < x <= 1; 1 / x => a; ln(a + sqrt(a^2 - 1)) :: ln((1 + sqrt(1 - x^2)) / x)
        __float128 r;
        if ((a <= (__float128)0.0) || (a > (__float128)1.0))
        {
                exit(1);
        }
        r = MultInvq((__float128)a);
        r = Logq(r + Sqrtq(r * r - (__float128)1.0));
        return r;
}

__float128 ArcSinhq(__float128 a)
{
        // arcsinh(x) = ln(x + sqrt(x^2 + 1))
        __float128 r;
        r = Sqrtq(a * a + (__float128)1.0);
        r = Logq(a + r);
        return r;
}

__float128 ArcTanhq(__float128 a)
{
        // arctanh(x) = abs(x) < 1; ln((1 + x)/(1 - x)) / 2
        __float128 r;
        if (Absq(a) >= (__float128)1.0)
        {
                exit(1);
        }
        r = Divq((__float128)(a + (__float128)1.0), (__float128)((__float128)1.0 - a));
        r = Divq((__float128)Logq(r), (__float128)2.0);
        return r;
}

__float128 Cotq(__float128 a)
{
        __float128 r;
        r = MultInvq((__float128)Tanq(a));
        return r;
}

__float128 Cothq(__float128 a)
{
        // coth(x) = x != 0; 1 / tanh(x)
        __float128 r;
        if (a == (__float128)0.0)
        {
                exit(1);
        }
        r = MultInvq((__float128)Tanh(a));
        return r;
}

__float128 Cscq(__float128 a)
{
        __float128 r;
        r = MultInvq((__float128)Sinq(a));
        return r;
}

__float128 Cschq(__float128 a)
{
        // csch(x) = x != 0; 1 / sinh(x)
        __float128 r;
        if (a == (__float128)0.0)
        {
                exit(1);
        }
        r = MultInvq((__float128)Sinh(a));
        return r;
}

__float128 Secq(__float128 a)
{
        __float128 r;
        r = MultInvq((__float128)Cosq(a));
        return r;
}

__float128 Sechq(__float128 a)
{
        // sech(x) = 1 / cosh(x)
        __float128 r;
        r = MultInvq((__float128)Cosh(a));
        return r;
}

__float128 RadiansToDegreesq(__float128 r)
{
        __float128 d;
        if (ACONST_PI == (__float128)0.0)
        {
                init_trig();
        }
        d = Divq((__float128)ACONST_PI, (__float128)2.0);
        d = Divq((__float128)r, (__float128)d) * (__float128)90;
        return d;
}

__float128 DegreesToRadiansq(__float128 d)
{
        __float128 r;
        if (ACONST_PI == (__float128)0.0)
        {
                init_trig();
        }
        r = Divq((__float128)ACONST_PI, (__float128)2.0);
        r *= DivAtom(d, (__float128)90);
        return r;
}

// end of file.
