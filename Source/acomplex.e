-- Copyright (c) 2025 James Cook
-- acomplex.e, accuracate complex calculations using sequences of two atoms for real and imaginary parts.
-- uses Newton's method.

-- Check your work twice before trusting it.
-- When in doubt about formulas, use WolframAlpha:
-- https://www.wolframalpha.com/

include amath.e
include std/math.e

global constant AREAL = 1, AIMAG = 2

global type AComplex(sequence s)
    return length(s) = 2 and atom(s[1]) and atom(s[2])
end type

global function CAdd(AComplex a, AComplex b)
    return a + b
end function

global function CSum(sequence s)
  AComplex r
  r = {0, 0}
  for i = 1 to length(s) do
    r = CAdd(r, s[i])
  end for
  return r
end function

global function CNegateImag(AComplex a)
    a[AIMAG] = -(a[AIMAG])
    return a
end function

global function ComplexConjugate(AComplex z)
    return CNegateImag(z)
end function

global function CMult(AComplex n1, AComplex n2)
    -- n1 = (a+bi)
    -- n2 = (c+di)
    -- (a+bi)(c+di) <=> ac + adi + bci + bdii
    -- <=> (ac - bd) + (ad + bc)i
    sequence r
    r = {0, 0}
    r[AREAL] = n1[AREAL] * n2[AREAL] - n1[AIMAG] * n2[AIMAG]
    r[AIMAG] = n1[AREAL] * n2[AIMAG] + n1[AIMAG] * n2[AREAL]
    return r
end function

global function CMultInv(AComplex n2)
-- Reciprocal
    -- Eun a, b, c
    -- (a+bi)(a-bi) <=> a*a + b*b
    -- n2 = (a+bi)
    -- a = n2[1]
    -- b = n2[2]
    -- 1 / n2 <=> (a-bi) / (a*a + b*b)
    -- <=> (a / (a*a + b*b)) - (b / (a*a + b*b))i
    -- c = (a*a + b*b)
    -- <=> (a / c) - (b / c)i
    atom a, b, c, real, imag
    a = n2[AREAL]
    b = n2[AIMAG]
    c = MultInvAtom((a * a) + (b * b))
    real = (a * c)
    imag = -(b * c)
    return {real, imag}
end function

global function CDiv(AComplex n1, AComplex n2)
-- correct!
    return CMult(n1, CMultInv(n2))
end function

global function CAbsoluteValue(AComplex z)
-- same as: ComplexModulus or ComplexMagnitude

-- Abs[z] == Sqrt[z Conjugate[z]]

-- abs(z) = sqrt(z * conj(z))
-- abs(x + iy) = sqrt((x + iy)(x - iy))

-- (x + iy)(x - iy) = x*x -x*iy +x*iy -iy*iy
-- x^2 - (iy)^2
-- x^2 + y^2

-- abs(x + iy) = sqrt(x^2 + y^2)

    atom r
    r = z[AREAL] * z[AREAL] + z[AIMAG] * z[AIMAG]
    r = SqrtAtom(r)
    return r
end function

global function CSqrt(AComplex z)
-- sqrt(x + iy) <=> (1/2) * sqrt(2) * [ sqrt( sqrt(x*x + y*y) + x )  +  i*sign(y) * sqrt( sqrt(x*x + y*y) - x ) ]
  atom x, y, a, b, c
  x = z[AREAL]
  y = z[AIMAG]
  c = SqrtAtom(x * x + y * y)
  a = SqrtAtom(c + x)
  b = SqrtAtom(c - x)
  if y < 0 then
    b = -(b)
  end if
  z = (SqrtAtom(2) / 2) * {a, b}
  return z
end function

global function CExp(AComplex z)
-- ComplexExponent()
    sequence r
    r = repeat(ExpAtom(z[AREAL]), 2)
    r[AREAL] *= CosAtom(z[AIMAG])
    r[AIMAG] *= SinAtom(z[AIMAG])
    return r
end function

global function CLn(AComplex z)
-- Natural Logarithm
-- ln(z) = (ln(x^2 + y^2)/2) + arctan(y/x)i
    sequence r
    r = {0, 0}
    r[AREAL] = LnAtom(z[AREAL] * z[AREAL] + z[AIMAG] * z[AIMAG]) / 2
    r[AIMAG] = ArcTan2Atom(z[AIMAG], z[AREAL])
    return r
end function

global function CPower(AComplex z, AComplex raisedTo)
    return CExp(CMult(raisedTo, CLn(z)))
end function

-- Complex Trig functions

-- For z = x + iy,

-- sin(z) = sin(x) * cosh(y) + i * cos(x) * sinh(y)
-- {\displaystyle \sin {z}=\sin {x}\cosh {y}+i\cos {x}\sinh {y}}

global function CSin(AComplex z)
-- sin(z) = (sin(x) * cosh(y)) + (cos(x) * sinh(y))i
    sequence r
    r = {0, 0}
    r[AREAL] = SinAtom(z[AREAL]) * CoshAtom(z[AIMAG])
    r[AIMAG] = CosAtom(z[AREAL]) * SinhAtom(z[AIMAG])
    return r
end function

-- cos(z) = cos(x) * cosh(y) − i * sin(x) * sinh(y)
-- {\displaystyle \cos {z}=\cos {x}\cosh {y}-i\sin {x}\sinh {y}}

global function CCos(AComplex z)
-- cos(z) = (cos(x) * cosh(y)) - (sin(x) * sinh(y))i
    sequence r
    r = {0, 0}
    r[AREAL] = CosAtom(z[AREAL]) * CoshAtom(z[AIMAG])
    r[AIMAG] = -(SinAtom(z[AREAL]) * SinhAtom(z[AIMAG]))
    return r
end function

-- tan(z) = (tan(x) + i * tanh(y)) / (1 − i * tan(x) * tanh(y))
-- {\displaystyle \tan {z}={\frac {\tan {x}+i\tanh {y}}{1-i\tan {x}\tanh {y}}}}

global function CTan(AComplex z)
-- z = Real(x) + Imaginary(y)
-- f(x,y) = cos(2x) + cosh(2y)
-- tan(z) = (sin(2x)/f(x,y)) + (sinh(2y)/f(x,y))i
  atom x, y, f
  sequence r
  x = z[AREAL] * 2
  y = z[AIMAG] * 2
  f = CosAtom(x) + CoshAtom(y)
  r = {0, 0}
  r[AREAL] = DivAtom(SinAtom(x), f)
  r[AIMAG] = DivAtom(SinhAtom(y), f)
  return r
end function

-- cot(z) = −((1 + i * cot(x) * coth(y)) / (cot(x) − i * coth(y)))
-- {\displaystyle \cot {z}=-{\frac {1+i\cot {x}\coth {y}}{\cot {x}-i\coth {y}}}}

global function CCot(AComplex z)
  atom a, b
  sequence r, n, d
  a = z[AREAL]
  b = z[AIMAG]
  a = CotAtom(a)
  b = CothAtom(b)
  n = {0, 0}
  n[AREAL] = 1
  n[AIMAG] = a * b
  d = {0, 0}
  d[AREAL] = a
  d[AIMAG] = - (b)
  r = - CDiv(n, d)
  return r
end function

-- Hyperbolic Trig functions
-- more information on this site:
-- https://mathworld.wolfram.com/HyperbolicFunctions.html

-- sinh(z) = sinh(x) * cos(y) + i * cosh(x) * sin(y)
-- {\displaystyle \sinh {z}=\sinh {x}\cos {y}+i\cosh {x}\sin {y}}

global function CSinh(AComplex z)
-- Sinus hyperbolic (sinh)
-- sinh(z) = (sinh(x) * cos(y)) - (cosh(x) * sin(y))i
    sequence r
    r = {0, 0}
    r[AREAL] = SinhAtom(z[AREAL]) * CosAtom(z[AIMAG])
    r[AIMAG] = CoshAtom(z[AREAL]) * SinAtom(z[AIMAG]) -- ? plus or minus
    return r
end function

-- cosh(z) = cosh(x) * cos(y) + i * sinh(x) * sin(y)
-- {\displaystyle \cosh {z}=\cosh {x}\cos {y}+i\sinh {x}\sin {y}}

global function CCosh(AComplex z)
-- Cosine hyperbolic (cosh)
-- cosh(z) = (cosh(x) * cos(y)) - (sinh(x) * sin(y))i
    sequence r
    r = {0, 0}
    r[AREAL] = CoshAtom(z[AREAL]) * CosAtom(z[AIMAG])
    r[AIMAG] = SinhAtom(z[AREAL]) * SinAtom(z[AIMAG]) -- ? plus or minus
    return r
end function

-- tanh(z) = (tanh(x) + i * tan(y)) / (1 + i * tanh(x) * tan(y))
// {\displaystyle \tanh {z}={\frac {\tanh {x}+i\tan {y}}{1+i\tanh {x}\tan {y}}}}

global function CTanh(AComplex z)
  atom a, b
  sequence r, n, d
  a = z[AREAL]
  b = z[AIMAG]
  a = TanhAtom(a)
  b = TanAtom(b)
  n = {0, 0}
  n[AREAL] = a
  n[AIMAG] = b
  d = {0, 0}
  d[AREAL] = 1
  d[AIMAG] = a * b
  r = CDiv(n, d)
  return r
end function

-- coth(z) = (1 − i * coth(x) * cot(y)) / (coth(x) − i * cot(y))
-- {\displaystyle \coth {z}={\frac {1-i\coth {x}\cot {y}}{\coth {x}-i\cot {y}}}}

global function CCoth(AComplex z)
  atom a, b
  sequence r, n, d
  a = z[AREAL]
  b = z[AIMAG]
  a = CothAtom(a)
  b = CotAtom(b)
  n = {0, 0}
  n[AREAL] = 1
  n[AIMAG] = - (a * b)
  d = {0, 0}
  d[AREAL] = a
  d[AIMAG] = - (b)
  r = CDiv(n, d)
  return r
end function


global function CArcTan(AComplex z)
-- Given: arctan(x + iy), z = x + iy
-- (1/2) * i * log(1 - i(x + iy)) - (1/2) * i * log(1 + i(x + iy))
-- (1/2) * i * (log(-ix + y + 1) - log(ix - y + 1))
-- (1/2) * i * (log(1 - i * z) - log(1 + i * z))

  sequence a, b, r
  a = CMult({0, 1}, z)
  b = CAdd(CLn(CAdd({1, 0}, -(a))), - CLn(CAdd({1, 0}, a)))
  r = CMult({0, 0.5}, b)
  return r
end function

global function CQuadraticEquation(AComplex a, AComplex b, AComplex c)
    -- The quadratic equation produces two answers (the answers may be the same)
    -- ax^2 + bx + c
    -- f(a,b,c) = (-b +-sqrt(b*b - 4*a*c)) / (2*a)
    -- answer[0] = ((-b + sqrt(b*b - 4*a*c)) / (2*a))
    -- answer[1] = ((-b - sqrt(b*b - 4*a*c)) / (2*a))
    --
    sequence r, s, t
    t = CMult(b, b) - 4 * CMult(a, c))
    s = CSqrt(t)
    r = {0, 0}
    r[1] = CDiv(-b + s, 2 * a)
    r[2] = CDiv(-b - s, 2 * a)
    return r
end function

-- end of file.
