-- Copyright (c) 2025 James Cook
-- amath.e
-- Description: Accurate math routines, using Newton's method.
-- Author: James Cook (jmsck55 AT gmail DOT com)

-- Routines for atoms and objects.

include std/math.e

global integer AMATH_DEBUG
AMATH_DEBUG = 0

ifdef BITS64 then
constant accuracy ≈ power(2, 64) -- question
elsedef
constant accuracy = power(2, 53) -- question
end ifdef

global constant AMATH_ITERS = 1000000

global function MultInvAtom(atom x)
-- performs 1/x
  atom t, g, last
  g = 1 / x
  last = 0
  for i = 1 to AMATH_ITERS do
    t = g * x
    if AMATH_DEBUG then
      if round(t, accuracy) = 1 then
        exit
      end if
      printf(2, "calc [%d] corrected %g\n", {i, t - 1})
    elsif t = 1 then
      exit
    end if
    if t = last then
      exit
    end if
    g *= ( 2 - t )
    last = t
  end for
  return g
end function

global function MultInv(object x)
    sequence s
    if atom(x) then
        return MultInvAtom(x)
    end if
    s = repeat(0, length(x))
    for i = 1 to length(x) do
        s[i] = MultInv(x[i])
    end for
    return s
end function

global function DivAtom(object n, atom d)
  if d = 2 then
    return n / 2
  end if
  n *= MultInvAtom(d)
  return n
end function

global function Div(object n, object d)
  sequence s
  if atom(d) then
    return DivAtom(n, d)
  end if
  s = repeat(0, length(d))
  for i = 1 to length(d) do
    s[i] = Div(n, d[i])
  end for
  return s
end function

global function NthRoot(object x, integer n)
  object quotient, average, guess
  guess = power(x, MultInvAtom(n))
  quotient = Div(x, power(guess, n-1))
  average = Div(quotient + ((n-1) * guess)), n)
  return average
end function

global function NthRootAtom(atom x, integer n)
  return NthRoot(x, n)
end function

global function Sqrt(object x)
  return NthRoot(x, 2)
end function

global function SqrtAtom(atom x)
  return Sqrt(x)
end function

global function Cbrt(object x)
  return NthRoot(x, 3)
end function

global function CbrtAtom(atom x)
  return Cbrt(x)
end function

global function ExpAtom(atom x)
-- using taylor series
-- https://en.wikipedia.org/wiki/TaylorSeries
--
-- -- exp(1) = sum of k=0 to inf (1/k!)
-- -- exp(x) = sum of k=0 to inf ((x^k)/k!)
--
  atom sum, num, den, last
  num = 1
  den = 1
  sum = 1
  last = 0
  for i = 1 to 1000000000 do
    num *= x
    den *= i -- number of iterations
    sum += DivAtom(num, den)
    if sum = last then
      exit
    end if
    last = sum
  end for
  return sum
end function

global function Exp(object x)
  sequence s
  if atom(x) then
    return ExpAtom(x)
  end if
  s = repeat(0, length(x))
  for i = 1 to length(x) do
    s[i] = Exp(x[i])
  end for
  return s
end function

-- Raw function: Natural Logarithm

global constant ACONST_E = ExpAtom(1)

global function LnAtom(atom a)
    -- Function: NaturalLogarithm()
    -- Use for testing the method.
    -- Alternative, between 0 and 2 exclusively:
    -- ln(x) = - Sum[k = 1 to inf] ((-1)^k * (-1 + x)^k) / k, for abs(-1 + x) < 1; x > 0 and x < 2;
    -- ((-1)^k * (-1 + x)^k) / k
    -- ((-x + 1)^k) / k
    -- Alternative, away from 0 to 2 exclusively: // Use factoring of "e" instead.
    -- ln(x) = ln(-1 + x) - Sum[k = 1 to inf] ((-1)^k * (-1 + x)^(-k)) / k, for abs(-1 + x) > 1, x < 0 or x > 2;
    --
    -- precalculate:
    -- xNegativePlusOne = (1 - x), then multiply, and store as "p".
    -- k = 1, (xNegativePlusOne^1) / 1
    -- k = 2, (xNegativePlusOne^2) / 2
    -- k = 3, (xNegativePlusOne^3) / 3
    -- k = 4, (xNegativePlusOne^4) / 4
    -- k = 5, (xNegativePlusOne^5) / 5
    -- Then, summate and negate, then return sum:
    -- while 1 do
    --  p *= xNegativePlusOne
    --  sum += p / k
    --  if k == inf then -- as k approaches infinity.
    --    exit -- break;
    --  end if
    --  k += 1
    -- end while
    -- return - (sum)
    atom x, p, sum, last
    integer f
    if a <= 0 then
      abort(1)
    elsif a >= 2 then
      -- factor, first find n, then calculate e^n
      -- log(m/e^n) + n = log(m)
      if a > ACONST_E then
        f = floor(log(a))
      else
        f = 1
      end if
      a = DivAtom(a, power(ACONST_E, f))
    else
      f = 0
    end if
    x = 1 - a
    last = 0
    p = 1
    sum = 0
    for k = 1 to 1000000000 do
      p *= x
      sum += DivAtom(p, k)
      if sum = last then
        exit
      end if
      last = sum
    end for
    return f - (sum)
end function

global function Ln(object x)
  sequence s
  if atom(x) then
    return LnAtom(x)
  end if
  s = repeat(0, length(x))
  for i = 1 to length(x) do
    s[i] = Ln(x[i])
  end for
  return s
end function

-- Trig functions

-- cos
-- sin
-- tan
-- arctan
-- arccos
-- arcsin
-- ACONST_PI

global function CosAtom(atom a)
-- cos(x) = 1 - ((x^2)/(2!)) + ((x^4)/(4!)) - ((x^6)/(6!)) + ((x^8)/(8!)) - ...
  atom f, r, d, c, x
  f = 1
  r = 1
  d = -1
  c = 0
  x = a
  for i = 2 to 1000000000 by 2 do
    f *= i
    x *= a
    r += (d * DivAtom(x, f))
    if r = c then
      exit
    end if
    c = r
    f *= (i + 1)
    x *= a
    d *= (-1)
  end for
  return r
end function

global function Cos(object x)
  sequence s
  if atom(x) then
    return CosAtom(x)
  end if
  s = repeat(0, length(x))
  for i = 1 to length(x) do
    s[i] = Cos(x[i])
  end for
  return s
end function

global function SinAtom(atom a)
-- sine(x) = x - ((x^3)/(3!)) + ((x^5)/(5!)) - ((x^7)/(7!)) + ((x^9)/(9!)) - ...
  atom f, r, d, c, x
  integer i
  f = 2
  r = a
  d = -1
  c = 0
  x = a * a
  for i = 3 to 1000000001 by 2 do
    f *= i
    x *= a
    r += (d * DivAtom(x, f))
    if r = c then
      exit
    end if
    c = r
    f *= (i + 1)
    x *= a
    d *= (-1)
  end for
  return r
end function

global function Sin(object x)
  sequence s
  if atom(x) then
    return SinAtom(x)
  end if
  s = repeat(0, length(x))
  for i = 1 to length(x) do
    s[i] = Sin(x[i])
  end for
  return s
end function

global function TanAtom(atom a)
  atom r
  r = DivAtom(SinAtom(a), CosAtom(a))
  return r
end function

global function Tan(object x)
  sequence s
  if atom(x) then
    return TanAtom(x)
  end if
  s = repeat(0, length(x))
  for i = 1 to length(x) do
    s[i] = Tan(x[i])
  end for
  return s
end function

-- arc functions

global function ArcTanAtom(atom a)
-- Begin ArcTanExpA()
--                z        +inf         n             2kz^2
-- arctan(z) = ------- * Sumation of Product of --------------------
--             1 + z^2     n=0         k=1      (2k + 1) * (1 + z^2)
--
-- (The term in the sum for n = 0 is the empty product, so is 1.)

-- ans = (z / (1 + z*z)) * Sum [ n=0 to +inf ] Prod [ k=1 to n ] (2kz*z) / ((2k + 1) * (1 + z*z))

  atom b, r, s, p, c
  b = a * a + 1
  s = 1
  for n = 1 to 1000000000 do
    p = 1
    for k = 1 to n do
      p *= DivAtom(k * 2 * (b - 1), (k * 2 + 1) * b)
    end for
    c = s
    s += p
    if s = c then
      exit
    end if
  end for
  r = DivAtom(a, b) * s
  return r
end function

global constant ACONST_PI = ArcTanAtom(1) * 4

global function ArcTan(object x)
  sequence s
  if atom(x) then
    return ArcTanAtom(x)
  end if
  s = repeat(0, length(x))
  for i = 1 to length(x) do
    s[i] = ArcTan(x[i])
  end for
  return s
end function

global function ArcSinAtom(atom a)
-- arcsin(x) = arctan( x / sqrt(1 - x^2) )
  atom r
  r = SqrtAtom(1 - a * a)
  r = DivAtom(a, r)
  r = ArcTanAtom(r)
  return r
end function

global function ArcSin(object x)
  sequence s
  if atom(x) then
    return ArcSinAtom(x)
  end if
  s = repeat(0, length(x))
  for i = 1 to length(x) do
    s[i] = ArcSin(x[i])
  end for
  return s
end function

global function ArcCosAtom(atom a)
-- arccos(x) = arctan( sqrt(1 - x^2) / x )
-- Limited domain: -1 to 1
-- Also:
--   arccos(x) = arcsin(1) - arcsin(x)
--   arccos(x) = (EunPi / 2) - arcsin(x)
  atom r
  r = SqrtAtom(1 - a * a)
  r = DivAtom(r, a)
  r = ArcTanAtom(r)
  return r
end function

global function ArcCos(object x)
  sequence s
  if atom(x) then
    return ArcCosAtom(x)
  end if
  s = repeat(0, length(x))
  for i = 1 to length(x) do
    s[i] = ArcCos(x[i])
  end for
  return s
end function

-- other trig functions

-- Cosh
-- Sinh
-- Tanh
-- ArcCosh
-- ArcCot
-- ArcCoth
-- ArcCsc
-- ArcCsch
-- ArcSec
-- ArcSech
-- ArcSinh
-- ArcTanh
-- Cot
-- Coth
-- Csc
-- Csch
-- Sec
-- Sech

global function CoshAtom(atom a)
-- cosh(x) = (e^(x) + e^(-x)) / 2
    atom r
    r = DivAtom(ExpAtom(a) + ExpAtom(-a), 2)
    return r
end function

global function Cosh(object x)
  sequence s
  if atom(x) then
    return CoshAtom(x)
  end if
  s = repeat(0, length(x))
  for i = 1 to length(x) do
    s[i] = Cosh(x[i])
  end for
  return s
end function

global function SinhAtom(atom a)
-- sinh(x) = (e^(x) - e^(-x)) / 2
    atom r
    r = DivAtom(ExpAtom(a) - EunExp(-a), 2)
    return r
end function

global function Sinh(object x)
  sequence s
  if atom(x) then
    return SinhAtom(x)
  end if
  s = repeat(0, length(x))
  for i = 1 to length(x) do
    s[i] = Sinh(x[i])
  end for
  return s
end function

global function TanhAtom(atom a)
-- tanh(x) = e^(2*x) => a; (a - 1) / (a + 1)
    atom r
    r = ExpAtom(a * 2)
    r = DivAtom(r - 1, r + 1)
    return r
end function

global function Tanh(object x)
  sequence s
  if atom(x) then
    return TanhAtom(x)
  end if
  s = repeat(0, length(x))
  for i = 1 to length(x) do
    s[i] = Tanh(x[i])
  end for
  return s
end function

global function ArcCoshAtom(atom a)
-- arccosh(x) = x >= 1; ln(x + sqrt(x^2 - 1))
  atom r
  if a < 1 then
    abort(1)
  end if
  r = SqrtAtom(a * a - 1)
  r = LnAtom(a + r)
  return r
end function

global function ArcCosh(object x)
  sequence s
  if atom(x) then
    return ArcCoshAtom(x)
  end if
  s = repeat(0, length(x))
  for i = 1 to length(x) do
    s[i] = ArcCosh(x[i])
  end for
  return s
end function

global function ArcCotAtom(atom a)
    atom r
    if a = 0 then
        return DivAtom(ACONST_PI, 2)
    end if
    r = ArcTanAtom(MultInvAtom(a))
    if a < 0 then
        r += ACONST_PI
    end if
    return r
end function

global function ArcCot(object x)
  sequence s
  if atom(x) then
    return ArcCotAtom(x)
  end if
  s = repeat(0, length(x))
  for i = 1 to length(x) do
    s[i] = ArcCot(x[i])
  end for
  return s
end function

global function ArcCothAtom(atom a)
-- arccoth(x) = abs(x) > 1; ln((x + 1)/(x - 1)) / 2
    atom r
    if abs(a) <= 1 then
        abort(1)
    end if
    r = DivAtom(a + 1, a - 1)
    r = DivAtom(LnAtom(r), 2)
    return r
end function

global function ArcCoth(object x)
  sequence s
  if atom(x) then
    return ArcCothAtom(x)
  end if
  s = repeat(0, length(x))
  for i = 1 to length(x) do
    s[i] = ArcCoth(x[i])
  end for
  return s
end function

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
    r = DivAtom(LnAtom(r), 2)
    return r
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

-- end of file.
