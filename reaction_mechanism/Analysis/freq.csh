#!/bin/csh

set time = "1000"
set freqs = "1 3 10 32 100 250 500 1000"
set freqs = "1"

set mols = " rc-0.10 rc-0.20 rc-0.30 rc-0.40 rc-0.50 rc-0.60 rc-0.70 rc-0.80 rc-0.90 rc-1.00 rc-1.10 rc-1.20 rc-1.30 rc-1.40 rc-1.50 rc0.00 rc0.10 rc0.20 rc0.30 rc0.40 rc0.50 rc0.60 rc0.70 rc0.80 rc0.90 rc1.00 rc1.10 rc1.20 rc1.30 rc1.40 rc1.50"
set mols = "rc-0.90 rc-0.80"

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
