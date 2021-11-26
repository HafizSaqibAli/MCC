#!/usr/bin/env python3

import math

file_name = "WMEigenValues.txt"
out_file_name = "WMEigenValues_out.txt"


atoms = set (["O", "C", "N", "H", "S", "P", "B", ""])
molecules = [[]]

#constants values
Kb = 1.38064852e-23
T = 300.0
hbar = 1.0545718e-34
c = 30000000000
Na = 6.022e+23

my_file = open(file_name, "r")
out_file = open(out_file_name, "w")

i = 0
for line in my_file:
    if line !='\n':
        molecules[i].append(line.strip('\n'))
    else:
        molecules.append([])
        i += 1
for n in range(len(molecules)-1):
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
	for value in eigenvalues:
		w = value/(Kb*T)
		w = math.sqrt(w)
		frequencies.append(w)
		w = w/(2*math.pi)/c
		wavenumbers.append(w)
	for freq in frequencies:
		interim = (hbar*freq)/(Kb*T)
		fraction.append(interim)
	for value in fraction:
		mid = (math.exp(value)-1)
		mid2 = (math.log(1-math.exp(-value)))
		mid3 = value/mid-mid2
		sums.append(mid3)
	S_qm = sum(sums)
	S_qm = Na * Kb * S_qm
	out_file.write("For molecules {0}, qm entropy is {1}\n".format(molecules[n][0], S_qm))

my_file.close()
out_file.close()

