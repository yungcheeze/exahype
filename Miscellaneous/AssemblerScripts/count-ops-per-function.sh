# just a small script composed by SvenK to determine
# the cost of NCP's pdematrixb_ vs the pdencp_

# First, disassemble. Do this with a release mode binary!

objdump -Sd ExaHyPE-Z4 > Z4-assembly.txt

# typical sizes are for 20MB binary -> 200MB disassembled txt.

# Use grep to find the starting lines of functions:

grep -inP '^[0-9a-z]+ <.+>:$' Z4-assembly.txt

# To determine the length of functions (numbers of assembler instructions),
# count the lines between each function. Grep this for the functions
# you are interested in:

awk '/^[0-9a-z]+ <.+>:$/{p=$2} {a[p]++} END{for(i in a) print i" "a[i]-1}' Z4-assembly.txt | grep -i pde

# This gives you helpful output like:

# <pdematrixb_>: 17993
# <pdesource_>: 5450
# <pdencp_>: 2894
# <pdeflux_>: 133

# So you know well what's expensive and what not.

