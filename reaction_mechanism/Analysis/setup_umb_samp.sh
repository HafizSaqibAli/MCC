#!/bin/bash

# Prepare for Umbrella sampling, WITHOUT spawning:
#  - create directory for each rc
#  - create md inputfile(s) for each rc (requires presence of 'template' $md_files)
#  - create restraint (.RST) files for each rc
#  - create job submission script that runs series of jobs (run_umb_samp.sh)
#  Hafiz Saqib Ali, 28/10/2019
#
##### THE FOLLOWING IS THE ESSENTIAL USER INPUT! #####
#
# Necessary input files
prmtop="mol.prmtop"
restart="mol.inpcrd"
md_files="md.in"
# Reaction coordinate values, restraint, definition
start_rc=-1.50
end_rc=1.50
step=0.10	# NB: if start_rc > end_rc, make step<0 
kumb=300	# Force constant, in ..
# Define expr OR rc_iat for rc definition 


# NB rc_expr is currently not always handled correctly by AMBER/sander. Use rc_iat instead.
#rc_line="restraint = \"coordinate(distance(6,1),2.0,distance(1,2),-2.81)\""
rc_line="iat=1,2,6,1, rstwt=1,-1,"
#
################### END USER INPUT ###################
#

#### Start preparing directories and files
qsub_file="run_umb_samp.sh"

echo "Preparing the following reaction coordinate values:"

# Write the job submission script, except for jobs
# Now PBS options are hard-coded, so users may need to change them
# Working dir
workdir=`pwd`
jobname=`basename $workdir`
jobname=$jobname"_wt"
# Guess time needed on 1 node / 8 procs
nodes=1
procs=1
#for ((i=$start; i<=$end; i++)); do hours=`echo "$hours + $rc_hrs" | bc`; done
#hrs=`printf "%.0f" $hours`

# Print the first part of the job script
#
printf "#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V

#$ -pe smp.pe 16

module load apps/intel-15.0/amber/16
#

# ---------------------------------------------------------------
# EXECUTE a number of amber jobs
" > $qsub_file


#### Start THE loop 
# In case 'bc' is not present, use hard-coded r1 and r4)
for i in `seq $start_rc $step $end_rc`; do
       rc=`printf "%3.2f" $i`
       r1=`echo "scale=2; $rc-10" | bc`
       r4=`echo "scale=2; $rc+10" | bc`
       #r1=-10.0
       #r4=10.0 
       printf "\t%s\n" $rc
# Continue
	mkdir rc$rc
	cd rc$rc 
	# Write .RST file 
	printf "# reaction coordinate d(C1-Cl2) - d(O6-C1)

&rst
$rc_line
r1=$r1, r2=$rc, r3=$rc, r4=$r4,
rk2=$kumb, rk3=$kumb,
/
 " > rc$rc.RST	
	# Write md inputfiles, using 'template' md files
	temp_count=0
	for template in $md_files; do
		temp_count=`expr $temp_count + 1`
		#file=`basename $template | sed -e "s,\.i,_rc$rc.i,g"`
		mdname=`basename $template | sed -e "s,\.i,,g"`
		sed -e "s/ifqnt=1,/ifqnt=1, nmropt=1,/g" ../$template > $mdname.i
		printf " &wt 
  type=\'DUMPFREQ\', istep1=50, 
 &end
 &wt 
  type=\'END\' 
 &end
  DISANG=rc$rc.RST
  DUMPAVE=rc${rc}_$temp_count.tra" >> $mdname.i
	        # Add execute lines to job submission script
        	printf "cd $workdir/rc$rc\n" >> ../$qsub_file
        	#printf "mpirun -i $mdname.i -o $mdname.log -r ../$restart\n" >> ../$qsub_file
		printf "mpirun -np 8 sander.MPI -O -i $mdname.i -o $mdname.log -p ../$prmtop -c ../$restart -frc $mdname.frc -x $mdname.crd -r $mdname.rst\n" >> ../$qsub_file 
		restart="rc$rc/$mdname.rst"
	done
	cd ..
done

echo "All done."
echo "   Carefully check the contents of the reaction coordinate directories created."
echo "   Check and, if neccesary, alter the job submission file: run_umb_samp.sh"
echo "Have fun umbrella sampling!"
