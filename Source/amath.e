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
  atom r, t, g
  g = 1 / x
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
    r = g * ( 2 - t )
    -- if r = g then -- question, is this necessary?
    --   exit
    -- end if
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

