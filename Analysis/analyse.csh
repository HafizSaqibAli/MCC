#!/bin/csh

set mols = "formamide formamide_g ammonia 1,4-dioxane 2-butaoxyethanol 2-methyl-2-propanol acetic-acid acetone acetophenone aniline benzene benzyl-alcohol benzyldehyde butane butanol carbon-dioxide chloroform cyclohexane diazene dichloromethane diethanolamine diethyl-ether dimethylformamide dmso ethane ethanol ethyl-acetate ethylamine ethylene-glycol formic-acid furan hexanol hydrazine hydrogen-peroxide methane methanethiol methanol methylamine m-xylene n-hexane n-methylacetamide octanol o-xylene pentane pentanol piperidine propane propanol p-xylene pyridine styrene tetrahydrofuran toluene triethylamine"
set mols = "butane"


foreach  mol ($mols)
     cd $mol/
#     echo $mol
#        ln -s mol.top molecules.top
#        ln -s md.frc forces.frc
#        ln -s md.mdcrd trajectory.crd

#	#Run the executable
#     ~/liquids/MCC/Exe 100 1

# Analyse
     set sconf = 0
     set ndih = 0
     if ( -f dih.dat ) then
        set sconf = `tail -1 dih.dat`
        set ndih = `wc dih.dat | awk '{print (0.25*($1-1))}'`
        echo $mol "Number of dihedrals $ndih"
     endif

     /home/hafiz/liquids/MCC/vib_ent.pl $mol f100/WMEigenValues.txt f100/FEigenValues.txt f100/TEigenValues.txt ../../../temperatures.dat $ndih
     /home/hafiz/liquids/bin/wm_config.pl $mol f1/WMShells.txt ../../../symmetry.dat

# Expt data
    # set sexp = `grep "s $mol" ../../../sexp.dat | awk '{print $3}' `
# Cell data
     set svib = `tail -1 $mol.svib | awk '{print $1}'`
     set strans = `tail -1 $mol.svib | awk '{print $2}'`
     set srot = `tail -1 $mol.svib | awk '{print $3}'`
     set sfint = `tail -1 $mol.svib | awk '{print $4}'`
     set stint = `tail -1 $mol.svib | awk '{print $5}'`
     set sorient = `tail -1 $mol.sconfig`
     set sconfig = `echo $sorient $sconf | awk '{print $1 + $2}'`
     set stotal = `echo $svib $sconfig  | awk '{print $1 + $2}'`
     echo $mol $stotal $sexp $strans $srot $sfint $stint  $sorient $sconf $svib $sconfig | awk '{printf "%12s %6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f  \n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}' 
#     echo $length $strans $s2pttrans $srot $s2ptrot $sint $s2ptint $sorient $sconf $svib $sconfig | awk '{printf "%5d %6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.  1f%6.1f\n", 10*$1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}' 

#tail -1 $mol.svib
        cd ../
end
