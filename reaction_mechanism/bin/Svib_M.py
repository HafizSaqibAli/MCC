#! /usr/bin/python
# coding: utf-8
# Hafiz Saqib
from __future__ import division  # Python 3 compatibility

import math

file_name = "f1/WMEigenValues.txt"

out_file_name = "Svib_M.txt"

atoms = set (["O", "C", "H", "B", "Cl", "F", ""])
molecules = [[]]

#Constants
Kb = 1.38064852e-23     # J/K
T = 298.15               # K
hbar = 6.6260755E-34    # Js
c = 30000000000         # m/s
Na = 6.022e+23          # mol
R = 8.3145              # J/K/mol

my_file = open(file_name, "r")
out_file = open(out_file_name, "w")

i = 0
for line in my_file:
    if line !='\n':
        molecules[i].append(line.strip('\n'))
    else:
        molecules.append([])
        i += 1
#print (len(molecules))

for n in range(0, len(molecules)-1):
    eigenvalues = []
    frequencies = []
    fraction = []
    wavenumbers = []
    sums = []
    for sl in molecules[n]:
        lines = sl.split('\n')
        for line in lines:
            if line[0:1] not in atoms:
                floats = (float(line))
                eigenvalues.append(floats)
            else:
                continue
    eigenvalues.sort(reverse = True)
    for value in eigenvalues:
        w = value / (Kb * T)
        w = math.sqrt(w)
        coeff = (1 / (2 * math.pi))
        w = coeff * w
        frequencies.append(w)
    for freq in frequencies:
        interm = (hbar * freq) / (Kb * T)
        fraction.append(interm)
    for value in fraction:
        mid = math.exp(value) - 1
        mid1 = math.log(1 - math.exp(-value))
        mid2 = (value / mid) - mid1
        sums.append(mid2)
    Srot = sum(sums[0:3])
    Strans = sum(sums[3:6])
    Srot_M = Srot * R
    Strans_M = Strans * R
    out_file.write("For molecules {0}, Srot at M_level in J/K/mol is {1}\n".format(molecules[n][0], Srot_M))
    out_file.write("For molecules {0}, Strans at M_level in J/K/mol is {1}\n".format(molecules[n][0], Strans_M))
my_file.close()
out_file.close()
