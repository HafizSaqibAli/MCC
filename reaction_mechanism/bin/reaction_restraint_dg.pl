#!/usr/bin/perl -w

# 22 Jan Calculates the restraint dG for the umbrella reaction calculations
# of Hafiz Sn2 reactions.
# Arguments; trajectory file, c, x, o (in the topology file, starting at 1), reaction coordinate value (Ang)

#use POSIX;
use PDL;
#use PDL::MatrixOps;

$natoms = 4507 ; # number of atoms. Be sure to set this.
$k = 300 ; # 
#$requil = -1.5*0; # Default equilibrium umbrella distance
$temp = 298.15;
$RT = $temp*8.3145/4.184;
$e = 2.718281828;

$crdfile = $ARGV[0] or die "Enter trajectory file\n";
$catom = $ARGV[1] or die "Enter c atom\n";
$xatom = $ARGV[2] or die "Enter halide atom \n";
$oatom = $ARGV[3] or die "Enter o atom \n";
$requil = $ARGV[4] or die "Enter rc value\n";
($catom, $xatom, $oatom) = ($catom-1, $xatom-1, $oatom-1);
printf "Number of atoms = %d, force constant = %d, C = %d, X = %d, O = %d, rc = %6.2f\n", $natoms, $k, $ARGV[1], $ARGV[2], $ARGV[3], $requil;

print "Reading $crdfile ...\n";
&read_crd( "$crdfile", \@coords );
#    &print_pdb( "whole".$ifile.".pdb", \@pdb, \@coords );
&calc_energy( \@coords );


#---------------------------------------------------------------------------
# SUBROUTINES

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
    $nframes = $ic / (3*($natoms+1));
    print "Read trajectory with $ic points and $nframes frames\n";
}

sub calc_energy {
    my ($cref) = @_;
    $boltzsum = 0;
    # Loop over all coordinates
    for $iframe (0 .. ($nframes - 1)) {
	$baseref = $iframe*3*($natoms+1);
        ($c[0], $c[1], $c[2]) = ( $$cref[$baseref+$catom*3],  $$cref[$baseref+$catom*3+1], $$cref[$baseref+$catom*3+2]);
	($x[0], $x[1], $x[2]) = ( $$cref[$baseref+$xatom*3],  $$cref[$baseref+$xatom*3+1], $$cref[$baseref+$xatom*3+2]);
	($o[0], $o[1], $o[2]) = ( $$cref[$baseref+$oatom*3],  $$cref[$baseref+$oatom*3+1], $$cref[$baseref+$oatom*3+2]);
	($box[0], $box[1], $box[2]) = ( $$cref[$baseref+$natoms*3],  $$cref[$baseref+$natoms*3+1], $$cref[$baseref+$natoms*3+2]);
#	printf "coords: %4d:  c:%8.3f%8.3f%8.3f   x:%8.3f%8.3f%8.3f  o:%8.3f%8.3f%8.3f  box:%8.3f%8.3f%8.3f \n", $iframe, $c[0], $c[1], $c[2], $x[0], $x[1], $x[2], $o[0], $o[1], $o[2],  $box[0], $box[1], $box[2];
	($dx,$dy,$dz) = ($c[0]-$x[0], $c[1]-$x[1], $c[2]-$x[2]);
#	printf "dx1:%8.3f%8.3f%8.3f ", $dx,$dy,$dz;
	$dx > 0.5*$box[0] and $dx -= $box[0];
	$dx < -0.5*$box[0] and $dx += $box[0];
	$dy > 0.5*$box[1] and $dy -= $box[1];
	$dy < -0.5*$box[1] and $dy += $box[1];
	$dz > 0.5*$box[2] and $dz -= $box[2];
	$dz < -0.5*$box[2] and $dz += $box[2];
#	printf " dx2:%8.3f%8.3f%8.3f\n", $dx,$dy,$dz;
	$cx = sqrt($dx**2 + $dy**2 + $dz**2);

	($dx,$dy,$dz) = ($c[0]-$o[0], $c[1]-$o[1], $c[2]-$o[2]);
	$dx > 0.5*$box[0] and $dx -= $box[0];
	$dx < -0.5*$box[0] and $dx += $box[0];
	$dy > 0.5*$box[1] and $dy -= $box[1];
	$dy < -0.5*$box[1] and $dy += $box[1];
	$dz > 0.5*$box[2] and $dz -= $box[2];
	$dz < -0.5*$box[2] and $dz += $box[2];
	$co = sqrt($dx**2 + $dy**2 + $dz**2);
#       $cx = sqrt(($c[0]-$x[0])**2 + ($c[1]-$x[1])**2 + ($c[2]-$x[2])**2);
#	$co = sqrt(($c[0]-$o[0])**2 + ($c[1]-$o[1])**2 + ($c[2]-$o[2])**2);
	# Calculate energy
	$umbrellaenergy = $k*($cx-$co-$requil)**2;
	$boltzfactor = $e**(-$umbrellaenergy/$RT);
	printf "%4d:  cx:%8.3f    co:%8.3f   diff:%8.3f   rc:%8.3f   u:%8.3f   e:%12.6f\n", $iframe, $cx, $co, $cx-$co, $cx-$co-$requil, $umbrellaenergy, $boltzfactor;
	$boltzsum += $boltzfactor;
    }
    $boltzsum /= $nframes;
    printf "exp = %12.6f dg = %12.3f\n", $boltzsum, -$RT*log($boltzsum);
}
