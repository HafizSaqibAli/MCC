#!/bin/csh

set time = "10000"
set freqs = "1 3 10 32 100 250 500 1000 10000"
set freqs = "100"
set mols = "water formamide formamide_g ammonia 1,4-dioxane 2-butaoxyethanol 2-methyl-2-propanol acetic-acid acetone acetophenone ammonia aniline benzene benzyl-alcohol benzyldehyde butane butanol carbon-dioxide chloroform cyclohexane diazene dichloromethane diethanolamine diethyl-ether dimethylformamide dmso ethane ethanol ethyl-acetate ethylamine ethylene-glycol formic-acid furan hexanol hydrazine hydrogen-peroxide methane methanethiol methanol methylamine m-xylene n-hexane n-methylacetamide octanol o-xylene pentane pentanol piperidine propane propanol p-xylene pyridine styrene tetrahydrofuran toluene triethylamine"
set mols = "100000f1ns"
#set mols = "10000f1ns"

foreach  mol ($mols)
        cd $mol
        echo $mol
        foreach freq ($freqs)
	   echo $freq
	   ../../MCC/Exe $time $freq

	   mkdir f${freq}
	   mv *.txt f${freq}
	end

        cd ..
end
