#!/bin/csh

#set mols = "formamide"

set mols = "formamide"
#set mols = "formamide"
 
set freqs = "1 3 10 32 100 250 500 1000"
set freqs = "100"

foreach  mol ($mols)
    cd $mol
    #ln -s 1000 f1
    echo $mol

    rm -f simfreq.dat freqcomponents.dat

    foreach freq ($freqs)

#if ( -d f$freq ) then
       cd f{$freq}

# Run the entropy calculations

        ../../../bin/vib_ent.pl $mol EigenValues.txt WMEigenValues.txt ../../../temperatures.dat ../../../internal.dat
        ../../../bin/wm_config.pl $mol WMShells.txt ../../../symmetry.dat

# Expt data
     set sexp = `grep "s $mol" ../../../sexp.dat | awk '{print $3}' `
# 2PT data
     set s2pt = `grep "s $mol" ../../../s2pt.dat | awk '{print $3}' `
     set s2pttrans = `grep "s $mol" ../../../s2pt.dat | awk '{print $6*4.184}' `
     set s2ptrot = `grep "s $mol" ../../../s2pt.dat | awk '{print $5*4.184}' `
     set s2ptint = `grep "s $mol" ../../../s2pt.dat | awk '{print $4*4.184}' `
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

     echo $freq $stotal $sexp | awk '{printf "%7.3f %7.2f %s \n", 100/$1, $2, $3}' >> ../simfreq.dat
#     echo $mol $strans $srot $sint $sorient $sconf | awk '{printf "%-12s%6.1f%6.1f%6.1f%6.1f%6.1f\n", $1, $2, $3, $4, $5, $6}' >> ../components.dat
     echo $freq $strans $s2pttrans $srot $s2ptrot $sint $s2ptint $sorient $sconf $svib $sconfig | awk '{printf "%7.3f %6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f\n", 100/$1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}' >> ../freqcomponents.dat
        cd ..
      endif
     end
     cd ..
end
