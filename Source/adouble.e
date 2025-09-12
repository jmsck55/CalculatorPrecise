-- Copyright (c) 2025 James Cook
-- Class file for doubles under 64-bit Euphoria

namespace double


include std/convert.e
include std/stack.e
include std/console.e

include amath.e
include acomplex.e

function encode(atom a)
-- encode a 64-bit double in as little memory as possible, as a 63-bit integer (without its smallest bit)
ifdef BITS64 then
atom b
integer f
sequence s
s = atom_to_float64(a)
b = s[$]
for i = length(s) - 1 to 1 by -1 do
  b *= #100
  b = or_bits(b, s[i])
end for
f = floor(b / 2)
return f
elsedef
return a
end ifdef
end function

function decode(atom f)
ifdef BITS64 then
f *= 2
return float64_to_atom(int_to_bytes(f, 8))
elsedef
return f
end ifdef
end function

public function format(atom a)
  return decode(encode(a))
end function

stack doubles = new()

public procedure store(atom val)
  push(doubles, encode(val))
end procedure

public function recall()
  return decode(pop(doubles))
end function

atom n, d, ans
n = 0
d = 0
ans = 0

public procedure display(atom arg, sequence p = "%.18g\n")
  printf(1, p, {arg})
end procedure

public procedure prompt(sequence p = ": ")
  puts(1, p)
end procedure



public function negate()
  ans = -(ans)
end function


public procedure start_calc()
sequence s

while 1 do

puts(1, "function (qnsr-/)") -- some functionality, not done yet.
prompt()
s = prompt_string()

switch s do
case "q" then
  exit
case "n" then
  ans = format(prompt_number("number: ", {}))
case "s" then
  store(ans)
case "r" then
  if is_empty(doubles) then
    puts(1, "nothing in storage.\n")
  else
    ans = recall()
  end if
case "-" then
  negate()
case "/" then
  d = format(prompt_number("number: ", {}))
  ans = format(DivAtom(ans, d))
case else
end switch
display(ans)

end while

end procedure
