#!/bin/csh

#set mols = "cyclohexane methanol propanol acetic_acid benzene ethane toluene acetone butane ethanol methane propane argon water"

set mols = "argon water methanol ethanol propanol methane ethane propane butane cyclohexane benzene toluene acetone acetic_acid"

rm -f results.txt components.txt
foreach  mol ($mols)
        cd $mol
#        echo $mol
#        ln -s mol.top molecules.top
#        ln -s md.frc forces.frc
#        ln -s md.mdcrd trajectory.crd
#	#Run the executable
#        ../hierarchical/vib_ent.pl $mol EigenValues.txt WMEigenValues.txt ../../temperatures.dat
        ../../bin/vib_ent.pl $mol EigenValues.txt WMEigenValues.txt ../../temperatures.dat internal.dat
#        ../../bin/simple_vib_ent.pl $mol EigenValues.txt WMEigenValues.txt ../../temperatures.dat
        ../../bin/wm_config.pl $mol WMShells.txt ../../symmetry.dat

# Expt data
     set sexp = `grep "s $mol" ../../sexp.dat | awk '{print $3}' `
# 2PT data
     set s2pt = `grep "s $mol" ../../s2pt.dat | awk '{print $3}' `
     set s2pttrans = `grep "s $mol" ../../s2pt.dat | awk '{print $6*4.184}' `
     set s2ptrot = `grep "s $mol" ../../s2pt.dat | awk '{print $5*4.184}' `
     set s2ptint = `grep "s $mol" ../../s2pt.dat | awk '{print $4*4.184}' `
# Cell data
     set svib = `tail -1 $mol.svib | awk '{print $1}'`
     set strans = `tail -1 $mol.svib | awk '{print $2}'`
     set srot = `tail -1 $mol.svib | awk '{print $3}'`
     set sint = `tail -1 $mol.svib | awk '{print $4}'`
     set sorient = `tail -1 $mol.sconfig`
     set sconf = 0
     if ( -f dih.dat ) then
        set sconf = `tail -1 dih.dat`
     endif
     set sconfig = `echo $sorient $sconf | awk '{print $1 + $2}'`
     set stotal = `echo $svib $sconfig  | awk '{print $1 + $2}'`
     echo $mol $stotal $sexp $s2pt | awk '{printf "%-12s%5.0f %s %5.0f\n", $1, $2, $3, $4*4.184}' >> ../results.txt
#     echo $mol $strans $srot $sint $sorient $sconf | awk '{printf "%-12s%6.1f%6.1f%6.1f%6.1f%6.1f\n", $1, $2, $3, $4, $5, $6}' >> ../components.txt
     echo $mol $strans $s2pttrans $srot $s2ptrot $sint $s2ptint $sorient $sconf $svib $sconfig | awk '{printf "%-12s%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f\n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}' >> ../components.txt
        cd ..
end
