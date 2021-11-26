#!/bin/csh

set mols = "butane butanol propanol 2-butaoxyethanol diethanolamine diethyl-ether ethyl-acetate hexanol triethyamine octanol pentanol pentane n-hexane"
set mols = "triethylamine octanol pentanol pentane n-hexane"

foreach  mol ($mols)
        cd $mol
#        echo $mol
     ../../bin/calc_dih.pl $mol trajectory.crd
     set sconf = `tail -1 dih.dat`
#set stotal = `echo $svib $sconfig | awk '{print $1 + $2}'`
     echo $mol $sconf
        cd ..
end
