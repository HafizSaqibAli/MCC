#!/bin/csh
set mols = "rc-0.10 rc-0.20 rc-0.30 rc-0.40 rc-0.50 rc-0.60 rc-0.70 rc-0.80 rc-0.90 rc-1.00 rc-1.10 rc-1.20 rc-1.30 rc-1.40 rc-1.50 rc0.00 rc0.10 rc0.20 rc0.30 rc0.40 rc0.50 rc0.60 rc0.70 rc0.80 rc0.90 rc1.00 rc1.10 rc1.20 rc1.30 rc1.40 rc1.50"

set mols = "rc-0.10 rc-0.20 rc-0.30 rc-0.40 rc-0.50 rc-0.60 rc-0.70 rc-0.80 rc-0.90 rc-1.00 rc-1.10 rc-1.20 rc-1.30 rc-1.40 rc-1.50 rc0.00 rc0.10 rc0.20 rc0.30 rc0.40 rc0.50 rc0.60 rc0.70 rc0.80 rc0.90 rc1.00 rc1.10 rc1.20 rc1.30 rc1.40 rc1.50"

foreach  mol ($mols)

	mv $mol/mdn.crd $mol/trajectory.crd
	mv $mol/mdn.frc $mol/forces.frc
	#mv $mol/mol.prmtop $mol/molecules.top
	#mv $mol/mden $mol/md.ene
end
