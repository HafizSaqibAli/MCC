#!/bin/csh

set mols = "1,4-dioxane 2-butaoxyethanol 2-methyl-2-propanol acetic-acid acetone acetophenone ammonia aniline benzene benzyl-alcohol benzyldehyde butane butanol carbon-dioxide chloroform cyclohexane diazene dichloromethane diethanolamine diethyl-ether dimethylformamide dmso ethane ethanol ethyl-acetate ethylamine ethylene-glycol formamide formic-acid furan hexanol hydrazine hydrogen-peroxide methane methanethiol methanol methylamine m-xylene n-hexane n-methylacetamide octanol o-xylene pentane pentanol piperidine propane propanol p-xylene pyridine styrene tetrahydrofuran toluene triethylamine"

foreach  mol ($mols)
        cd $mol
        echo $mol
        mkdir ../../graphs/$mol
#        ../../bin/format_nc.pl f1/WMShells.txt ../../graphs/$mol/ncoord.txt
 #       ../../bin/format_frequencies.pl $mol f1/TEigenValues.txt f1/TWMEigenValues.txt ../../temperatures.dat ../../graphs/$mol/freq.txt
  #      ~/liquids/bin/format_matrix.pl f1/FCMatrix.txt ../../graphs/$mol/fcmatrix.txt
        ~/liquids/bin/format_matrix.pl f1/TCMatrix.txt ../../graphs/$mol/tcmatrix.txt
   #     ~/liquids/bin/reorder_eigenvectors.pl f1/TEigenVectors.txt f1/TEigenValues.txt ../../graphs/$mol/teigenvectors.txt
        cd ..
end
