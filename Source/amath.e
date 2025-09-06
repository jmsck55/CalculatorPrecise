-- Copyright (c) 2025 James Cook
-- amath.e
-- Description: Accurate math routines, using Newton's method.
-- Author: James Cook (jmsck55 AT gmail DOT com)


include std/math.e

global integer AMATH_DEBUG
AMATH_DEBUG = 0

ifdef BITS64 then
constant accuracy ≈ power(2, 64) -- question
elsedef
constant accuracy = power(2, 53) -- question
end ifdef

global constant AMATH_ITERS = 10000

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

global function DivAtom(atom n, atom d)
  atom r
  r = MultInvAtom(d) * n
  return r
end function

global function NthRootAtom(atom x, integer n)
  atom quotient, average, guess
  guess = power(x, MultInvAtom(n))
  quotient = DivAtom(x, power(guess, n-1))
  average = DivAtom(quotient + ((n-1) * guess)), n)
  return average
end function

global function SqrtAtom(atom x)
  return NthRootAtom(x, 2)
end function

global function CbrtAtom(atom x)
  return NthRootAtom(x, 3)
end function

-- Raw function: Natural Logarithm

global function NaturalLogarithmAtom(atom a)
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
    if a <= 0 then
      abort(1)
    elsif a >= 2 then
      return {} -- to do later.
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
    return - (sum)
end function

global function NaturalExponentiationAtom(atom x)
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

-- more functions to come.
