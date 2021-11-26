#!/usr/bin/perl -w

# Calculates the dihedral distribution for multiple dihedrals of > four-atom systems
# Takes into account correlations between dihedrals. Amazing!
# Prints out eigenvalues and the entropy using the PDL eigenvalue method
#use POSIX;
use PDL;
use PDL::MatrixOps;

$nmol = 500; # number of molecules. Be sure to set this.

$molecule = $ARGV[0] or die "Enter molecule \n";
$crdfile = $ARGV[1] or die "Enter trajectory file\n";
$rad2deg = 180/3.141592654;

&define_dihedrals( );

$ndih = $#heavyatoms + 1; # Number of dihedrals
$nconf = 3; # Number of conformations per dihedral

#$n_bins = 3;
#$b_range = 360/$n_bins;

#&calc_dih( \@initcoords, \@start_dihedrals, \@dih_time, \%legal_dihedral );
#@dihedrals = @prev_dihedrals = @start_dihedrals;
#&print_dih_status( "dih0.dat", \@start_dihedrals, \@prev_dihedrals, \@dihedrals );

print "Reading $crdfile for $molecule ...\n";
&read_crd( "$crdfile", \@coords );
#    &print_pdb( "whole".$ifile.".pdb", \@pdb, \@coords );
&calc_dih( \@coords, \@heavyatoms, \@dihedrals, \@dihedral_single, \@dihedral_corr_matrix );
#&print_dih_status( "data/d".$ifile.".dat", \@start_dihedrals, \@prev_dihedrals, \@dihedrals );

&print_correlation_matrix( \@dihedral_single, \@dihedral_corr_matrix );
&print_dihedrals( \@dihedrals );


#---------------------------------------------------------------------------
# SUBROUTINES

# Define the sets of four atoms for each molecule
sub define_dihedrals {
    if ($molecule eq "butanol") {
# CCCC
    $id = 0;
    for $ia ("1", "5", "8", "11") { 
	$heavyatoms[0]{$ia} = ++$id;
    }
    $id = 0;
# HCCC
    for $ia ("2", "1", "5", "8") { 
	$heavyatoms[1]{$ia} = ++$id;
    }
    $natoms = 14; 
} elsif ($molecule eq "propanol") {
# CCCO
    $id = 0;
    for $ia ("1", "2", "3", "4") { 
	$heavyatoms[0]{$ia} = ++$id;
    }
    $natoms = 12; 
} elsif ($molecule eq "butane") {
# CCCC
    $id = 0;
    for $ia ("1", "2", "3", "4") { 
	$heavyatoms[0]{$ia} = ++$id;
    }
    $natoms = 14;
 } elsif ($molecule eq "2-butaoxyethanol") {
# CCCC
    $id = 0;
    for $ia ("1", "2", "3", "4") { 
	$heavyatoms[0]{$ia} = ++$id;
    }
# CCCC
    $id = 0;
    for $ia ("2", "3", "4", "5") { 
	$heavyatoms[1]{$ia} = ++$id;
    }
# CCCC
    $id = 0;
    for $ia ("3", "4", "5", "6") { 
	$heavyatoms[2]{$ia} = ++$id;
    }
# CCCO
    $id = 0;
    for $ia ("4", "5", "6", "7") { 
	$heavyatoms[3]{$ia} = ++$id;
    }
# CCCC
    $id = 0;
    for $ia ("5", "6", "7", "8") { 
	$heavyatoms[4]{$ia} = ++$id;
    }
    $natoms = 22;
} elsif ($molecule eq "diethanolamine") {
# OCCN
    $id = 0;
    for $ia ("1", "2", "3", "4") { 
	$heavyatoms[0]{$ia} = ++$id;
    }
# CCNC
    $id = 0;
    for $ia ("2", "3", "4", "5") { 
	$heavyatoms[1]{$ia} = ++$id;
    }
# CNCC
    $id = 0;
    for $ia ("3", "4", "5", "6") { 
	$heavyatoms[2]{$ia} = ++$id;
    }
# NCCO
    $id = 0;
    for $ia ("4", "5", "6", "7") { 
	$heavyatoms[3]{$ia} = ++$id;
    }
    $natoms = 18;
} elsif ($molecule eq "diethyl-ether") {
# CCOC
    $id = 0;
    for $ia ("1", "2", "3", "4") { 
	$heavyatoms[0]{$ia} = ++$id;
    }
# CCOC
    $id = 0;
    for $ia ("2", "3", "4", "5") { 
	$heavyatoms[1]{$ia} = ++$id;
    }
    $natoms = 15;
} elsif ($molecule eq "ethyl-acetate") {
# CCOC
    $id = 0;
    for $ia ("1", "2", "3", "4") { 
	$heavyatoms[0]{$ia} = ++$id;
    }
# COCC 
    $id = 0;
    for $ia ("6", "2", "3", "4") { 
	$heavyatoms[1]{$ia} = ++$id;
    }
# OCOC
    $id = 0;
    for $ia ("2", "3", "4", "5") { 
	$heavyatoms[2]{$ia} = ++$id;
    }
    $natoms = 14;
} elsif ($molecule eq "hexanol") {
# CCCC
    $id = 0;
    for $ia ("1", "2", "3", "4") { 
	$heavyatoms[0]{$ia} = ++$id;
    }
# CCCC
    $id = 0;
    for $ia ("2", "3", "4", "5") { 
	$heavyatoms[1]{$ia} = ++$id;
    }
# CCCC
    $id = 0;
    for $ia ("3", "4", "5", "6") { 
	$heavyatoms[2]{$ia} = ++$id;
    }
# CCCO
    $id = 0;
    for $ia ("4", "5", "6", "7") { 
	$heavyatoms[3]{$ia} = ++$id;
    }
    $natoms = 21;
} elsif ($molecule eq "triethylamine") {
# CCCC
    $id = 0;
    for $ia ("7", "6", "1", "4") { 
	$heavyatoms[0]{$ia} = ++$id;
    }
# CCCC
    $id = 0;
    for $ia ("3", "2", "1", "6") { 
	$heavyatoms[1]{$ia} = ++$id;
    }
# CCCC
    $id = 0;
    for $ia ("5", "4", "1", "2") { 
	$heavyatoms[2]{$ia} = ++$id;
    }
    $natoms = 22;
} elsif ($molecule eq "octanol") {
# CCCC
    $id = 0;
    for $ia ("1", "2", "3", "4") { 
	$heavyatoms[0]{$ia} = ++$id;
    }
# CCCC
    $id = 0;
    for $ia ("2", "3", "4", "5") { 
	$heavyatoms[1]{$ia} = ++$id;
    }
# CCCC
    $id = 0;
    for $ia ("3", "4", "5", "6") { 
	$heavyatoms[2]{$ia} = ++$id;
    }
# CCCC
    $id = 0;
    for $ia ("4", "5", "6", "7") { 
	$heavyatoms[3]{$ia} = ++$id;
    }
# CCCC
    $id = 0;
    for $ia ("5", "6", "7", "8") { 
	$heavyatoms[4]{$ia} = ++$id;
    }
# CCCO
    $id = 0;
    for $ia ("6", "7", "8", "9") { 
	$heavyatoms[5]{$ia} = ++$id;
    }
    $natoms = 27;
} elsif ($molecule eq "pentanol") {
# CCCC
    $id = 0;
    for $ia ("1", "2", "3", "4") { 
	$heavyatoms[0]{$ia} = ++$id;
    }
# CCCC 
    $id = 0;
    for $ia ("2", "3", "4", "5") { 
	$heavyatoms[1]{$ia} = ++$id;
    }
# CCCO
    $id = 0;
    for $ia ("3", "4", "5", "6") { 
	$heavyatoms[2]{$ia} = ++$id;
    }
    $natoms = 18;
} elsif ($molecule eq "hexane") {
# CCCC
    $id = 0;
    for $ia ("1", "2", "3", "4") {
        $heavyatoms[0]{$ia} = ++$id;
    }
# CCCC
    $id = 0;
    for $ia ("2", "3", "4", "5") {
        $heavyatoms[1]{$ia} = ++$id;
    }
# CCCC
    $id = 0;
    for $ia ("3", "4", "5", "6") {
        $heavyatoms[2]{$ia} = ++$id;
    }
    $natoms = 20;
} elsif ($molecule eq "pentane") {
# CCOC
    $id = 0;
    for $ia ("1", "2", "3", "4") { 
	$heavyatoms[0]{$ia} = ++$id;
    }
# CCOC
    $id = 0;
    for $ia ("2", "3", "4", "5") { 
	$heavyatoms[1]{$ia} = ++$id;
    }
    $natoms = 17;
} else {
    die "No/Invalid molecule entered\n";
}

}

sub read_crd {
    my ($file, $cref) = @_;
    open FILE, "$file" or die "Can't find $file\n" ;
    $_ = <FILE>;
    $ic = 0;
    LINE: while (<FILE>) {
	@line = (split);
	foreach $num (@line) {
	    $$cref[$ic++] = $num;
#	    $ic > $n_crds and last LINE;
	}
    }
    close FILE;
    print "Read trajectory with $ic atoms\n";
}

sub calc_dih {
    my ($cref, $haref, $dref, $dsingref, $dcmref) = @_;
    # Loop over all coordinates
    $moltot = 0;
      $ia = 1; # atom count
      $ixyz = 0; # xyz count
      $imol = 0; # molecule count
      $ic = 0; # coordinate count
      for  $idih (0 .. $#$haref) {
	  $id[$idih] = 0; # index of atom in dihedral
      }
# Loop over all trajectory frames
      IC: while ($ic <= $#$cref) {
# Loop over all dihedrals
        for $idih (0 .. $#$haref) {
	    if ($$haref[$idih]{$ia}) { # Save atom if heavy
		$p[$idih][$$haref[$idih]{$ia}-1][$ixyz] = $$cref[$ic];
#	        printf "%8.3f , %6d at $ic $ia $ixyz\n", $$cref[$ic],  $$haref[$idih]{$ia}-1;
		($ixyz == 2) and $id[$idih]++; # increment atom count
	    }
	}
	$ixyz++; # increment xyz counter
	if ($ixyz == 3) { # Increment if new atom
	    $ixyz = 0;
	    $ia++;
	}
	if ($ia > $natoms) { # Reached end of molecule, calc dih
	    $ia = 1; # Initialise counters for next molecule
	    for  $idih (0 .. $#$haref) {
		$id[$idih] = 0; # index of atom in dihedral
	    }
	    $imol++;
	    undef %dstatus; # Hash saving each dihedral
	    # Calculate dihedral angles for each dihedral
	    for $idih (0 .. $#$haref) {
#		print "$ic $imol $idih\n";
#		for $ia2 (0 .. 3) {
#		    for $ixyz2 (0 .. 2) {
#			printf "%8.3f", $p[$idih][$ia2][$ixyz2];
#		    }
#		    printf "\n";
#		}
		for $ix (0 .. 2) {
		    $vim1[$ix] = $p[$idih][1][$ix] - $p[$idih][0][$ix];
		    $vi[$ix] = $p[$idih][2][$ix] - $p[$idih][1][$ix];
		    $vip1[$ix] = $p[$idih][3][$ix] - $p[$idih][2][$ix];
		}
# First normal
		$nim1[0] = $vim1[1]*$vi[2] - $vi[1]*$vim1[2];
		$nim1[1] = $vim1[2]*$vi[0] - $vi[2]*$vim1[0];
		$nim1[2] = $vim1[0]*$vi[1] - $vi[0]*$vim1[1];
		$rmag = 1/sqrt($nim1[0]*$nim1[0] + $nim1[1]*$nim1[1] + $nim1[2]*$nim1[2]);
		for $ix (0 .. 2) {
		    $nim1[$ix] *= $rmag;
		}
# Second normal
		$ni[0] = $vi[1]*$vip1[2] - $vip1[1]*$vi[2];
		$ni[1] = $vi[2]*$vip1[0] - $vip1[2]*$vi[0];
		$ni[2] = $vi[0]*$vip1[1] - $vip1[0]*$vi[1];
		$rmag = 1/sqrt($ni[0]*$ni[0] + $ni[1]*$ni[1] + $ni[2]*$ni[2]);
		for $ix (0 .. 2) {
		    $ni[$ix] *= $rmag;
		}
		$dotprod = ($ni[0]*$nim1[0] + $ni[1]*$nim1[1] + $ni[2]*$nim1[2]);
# Calculate cross product
		$cross[0] = $ni[1]*$nim1[2] - $nim1[1]*$ni[2];
		$cross[1] = $ni[2]*$nim1[0] - $nim1[2]*$ni[0];
		$cross[2] = $ni[0]*$nim1[1] - $nim1[0]*$ni[1];

		$sign = $cross[0]*$vi[0] + $cross[1]*$vi[1] + $cross[2]*$vi[2];
		$phi = $rad2deg*acos($dotprod);
		$sign < 0 and $phi *= -1;
#	$phi += 180;
#        ($phi > 360 or $phi < -180) and printf "Phi = $phi ($ic) outside normal range\n";
#	$phi >= 180 and $phi -= 360;
# Create dihedral histogram
		if (($phi >= 120) or ($phi < -120)) {
		    $$dref[$idih][0]++;
		    $dstatus{$idih} = 0;
		} elsif (($phi >= 0) and ($phi < 120)) {
		    $$dref[$idih][1]++;
		    $dstatus{$idih} = 1;
		} elsif (($phi >= -120) and ($phi < 0)) {
		    $$dref[$idih][2]++;
		    $dstatus{$idih} = 2;
		}
#		(($idih == 0) and ($dstatus{$idih} == 2)) and print "X";
#		printf "$idih %8d%6d phi= %10.3f %3s\n", $ic, $imol, $phi, $dstatus{$idih};
	    }
 # Reset molecule count
	    if ($imol == $nmol) {
		$imol = 0;
#		printf "box = %8.3f%8.3f%8.3f\n", $$cref[$ic+1], $$cref[$ic+2], $$cref[$ic+3];
		$ic+=3; # skip box
	    }
# Calculate correlation matrix
	    for $idih (0 .. $#$dref) { 
		$istatus = $dstatus{$idih};
		$i = $idih*$nconf + $istatus;
		$$dsingref[$i] += 1;
                for $jdih (0 .. $#$dref) { 
		    $jstatus = $dstatus{$jdih};
		    $j = $jdih*$nconf + $jstatus;
                    $$dcmref[$i][$j] += 1;
		}
	    }
	    $moltot ++;
#	    $imol > 0 and last IC;
	}
	$ic++;
      }
}

sub print_correlation_matrix {
    my ($dsingref, $dcmref) = @_;
    ($moltot) ? ($invmoltot = 1/$moltot) : ($invmoltot = 0);
    printf "molecules $moltot and normalisation $invmoltot\n";
# Print dihedral probabilities
    for $i (0 .. $#$dsingref) { 
#	$dsingav = 0.726*$moltot;
	if ($$dsingref[$i]) {            
	    $sigma[$i] = $$dsingref[$i]*$invmoltot- ($$dsingref[$i]*$invmoltot)**2;
	    printf "%8d %10.4f %10.4f %10.4f\n", $$dsingref[$i], $$dsingref[$i]*$invmoltot, ($$dsingref[$i]*$invmoltot)**2, $sigma[$i];
	} else {
	    printf "%8d %10.4f\n", 0, 0;
	}
    }
    printf "\n";
# N dihedrals with M conformations gives and (N+M)x(N+M) correlation matrix
    $matrixtotal = 0;
    for $i (0 .. $#$dcmref) { 
	for $j (0 .. $#$dcmref) {
	    if (($$dcmref[$i][$j])) { # and ($i==$j)) {
		if ($i==$j) {
		    printf " %7d", $$dcmref[$i][$j]- $$dsingref[$i]*$$dsingref[$j]*$invmoltot;
#		    $$dcmref[$i][$j] = ($$dcmref[$i][$j]*$invmoltot) - ($$dsingref[$i]*$$dsingref[$j]*$invmoltot*$invmoltot);
		    $$dcmref[$i][$j] = $$dcmref[$i][$j]*$invmoltot*($$dcmref[$i][$j]*$invmoltot - ($$dsingref[$i]*$$dsingref[$j]*$invmoltot*$invmoltot))/sqrt($sigma[$i]*$sigma[$j]);
		} else {
		    printf " %7d", $$dcmref[$i][$j] -  $$dsingref[$i]*$$dsingref[$j]*$invmoltot;
#		    $$dcmref[$i][$j] = ($$dcmref[$i][$j]*$invmoltot) - ($$dsingref[$i]*$$dsingref[$j]*$invmoltot*$invmoltot);
		    $$dcmref[$i][$j] = $$dcmref[$i][$j]*$invmoltot*($$dcmref[$i][$j]*$invmoltot - ($$dsingref[$i]*$$dsingref[$j]*$invmoltot*$invmoltot))/sqrt($sigma[$i]*$sigma[$j]);
		}
	    } else {
		printf " %7d", 0;
		$$dcmref[$i][$j] = 0;# - 0*$$dsingref[$i]*$$dsingref[$j]*$invmoltot*$invmoltot;
	    }
	    $matrixtotal += $$dcmref[$i][$j];
	}
	    printf "\n";
    }
# Print normalised matrix
    ($matrixtotal) ? ($matnorm = $ndih/$matrixtotal) : ($matnorm = 0);
#    $matnorm = 1;
    printf "Matrix norm = $matnorm\n";
    printf "[";
    for $i (0 .. $#$dcmref) { 
	for $j (0 .. $#$dcmref) {
	    $$dcmref[$i][$j] *= $matnorm;
	    printf "%7.3f", $$dcmref[$i][$j];#*$matnorm;
	}
	printf ";\n";
    }
    printf "]\n";
}

# Obsolete except for a track
sub partial_print_correlation_matrix {
    my ($dsingref, $dcmref) = @_;
    ($moltot) ? ($invmoltot = 1/$moltot) : ($invmoltot = 0);
    for $i (0 .. $#$dsingref) { 
	if ($$dsingref[$i]) {
	    printf "%10.4f\n", $$dsingref[$i]*$invmoltot;
	} else {
	    printf "%10.4f\n", 0;
	}
    }
    printf "\n";
# N dihedrals with M conformations gives and (N+M)x(N+M) correlation matrix
    printf "[";
    for $i (0 .. $#$dcmref) { 
	foreach $j (0 .. $#$dcmref) {
            if (($$dcmref[$i][$j])) { # and ($i==$j)) {
		if ($i==$j) {
		    printf "%6.2f", $$dcmref[$i][$j]*$invmoltot - 0*$$dsingref[$i]*$$dsingref[$j]*$invmoltot*$invmoltot;
		} else {
		    printf "%6.2f", $$dcmref[$i][$j]*$invmoltot - 0*$$dsingref[$i]*$$dsingref[$j]*$invmoltot*$invmoltot;
		}
	    } else {
		printf "%6.2f", 0 - 0*$$dsingref[$i]*$$dsingref[$j]*$invmoltot*$invmoltot;
	    }
	}
	printf ";\n";
    }
    printf "]\n";
}

sub print_dihedrals {
    my ($dref) = @_;
    $outfile = "corr_dih.dat";
    open PHI, ">$outfile";
    $R = 8.3145;
    $stot = 0; # entropy total
# Calculate entropy for independent dihedrals
    printf PHI "Uncorrelated entropy \n";
    for $idih (0 .. $#$dref) { 
	$s = 0;
      $tot = 0; # normalisation
      foreach $dih (0 .. $#{$$dref[$idih]}) {
	$$dref[$idih][$dih] and $tot += $$dref[$idih][$dih];
      }
      ($tot)? ($norm = 1/$tot) : ($norm = 0);
      foreach $dih (0 .. $#{$$dref[$idih]}) {
	if ($$dref[$idih][$dih]) {
	    $prob = $$dref[$idih][$dih]*$norm;
	    printf PHI "%10.4f%10d\n", $prob, $$dref[$idih][$dih];
	    $s += -$prob*log($prob);
	} else {
	    printf PHI "%10.4f%10d\n", 0, 0;
	}
      }
      printf PHI "%8.3f\n", $R*$s;
	$stot += $s;
    }
    printf PHI "%8.3f\n", $R*$stot;
    printf PHI "\n";
# Diagonalise the matrix and get eigenvalues
    my $matrix = pdl(\@dihedral_corr_matrix);
#($ev,$e) = eigens_sym $a;        
    $eigenvalues = eigens_sym($matrix);
#print $eigenvalues;
    $stot = 0;
    $neigs = nelem($eigenvalues)-1;
    printf PHI "Correlated entropy \n";
    for $i (0 .. $neigs) {
	$eig = $eigenvalues->index($i);
	if ($eig) {
	    $stot += $eig*log($eig);
	}
	printf PHI "%10.5f %8.3f\n",  $eig, -$R*$stot;
    }
    printf PHI "%8.3f\n", -$R*$stot;
    close PHI;
}

