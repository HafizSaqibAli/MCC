#!/usr/bin/perl -w

# Format matrix to make easier to read and shuffle torques/fragments

$infile = "$ARGV[0]" or die "no eigfile\n";
$outfile = "$ARGV[1]" or die "no outfile\n";
  
&read_data( $infile, \@data, \%ic2newref );
&print_entropy( $outfile, \@data, \%ic2newref ); 

#--------------------------------------------------------------------------

sub read_data {
    my ($file, $dref, $ic2newref) = @_;
    open FILE, "$file" or die "Can't find $file\n";
    $_ = <FILE>; # Skip first line
    @atomlist = (split);
    $ic = 0;
    $total_force = 0;
    $total_torque = 0;
    for $ia (0 .. $#atomlist) {
	$ntorque = 3;
	$atom = $atomlist[$ia];
#	printf "$ia %s\n", $atom;
	(($atom eq "OH,") or ($atom eq "CH,")) and $ntorque = 2;
	(($atom eq "C,") or ($atom eq "O,") or ($atom eq "Ar,")) and $ntorque = 0;
        for $ixyz (0 .. 2) {
	    $$ic2newref{$ic} = $ia*3 + $ixyz;
	    printf "$ic f goes to %2d\n", $$ic2newref{$ic};
	    $ic++;
	}
        for $ixyz (0 .. ($ntorque-1)) {
	    $$ic2newref{$ic} = ($#atomlist+1)*3 + $total_torque + $ixyz;
	    printf "$ic t goes to %2d\n", $$ic2newref{$ic};
	    $ic++;
	}
	$total_force+= 3;
	$total_torque+= $ntorque;
    }
    printf "tot f: %5d  tot t: %5d\n", $total_force, $total_torque;
 
# Read in matrix
    RST: while (<FILE>) {
	chomp;
	$_ or last; #Quit on first empty line 
	push @$dref, [split];
    }
    close FILE;
}

# Print out the pair data as an unwrapped matrix
sub print_entropy {
    my ($file, $dref, $ic2newref) = @_;
# Renumber matrix using the mapping
    for $i  (0 .. $#$dref) {
	for $j (0 .. $#{$$dref[$i]}) {
            $inew = $$ic2newref{$i};
            $jnew = $$ic2newref{$j};
	    $atommatrix[$inew][$jnew] = $$dref[$i][$j];
#	    printf "%3d%3d%3d%3d %12.0f\n", $i, $j, $inew, $jnew, $$dref[$i][$j];
	}
    }
# print out new matrix
    open FILE, ">$file";
    for $i (0 .. $#atommatrix) {
	for $j (0 .. $#{$atommatrix[$i]}) {
	    printf FILE "%6d", $atommatrix[$i][$j]*1e-4;
	}
	printf FILE "\n";
    }
    close FILE;
}

