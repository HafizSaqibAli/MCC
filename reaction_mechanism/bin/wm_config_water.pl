#!/usr/bin/perl -w

# Determines water rotational probabilities based on the WMShells.txt file
# ... and possibly Shells.txt with atomc-specific information

#$molecule = "$ARGV[0]" or die "enter molecule\n";
$shellfile = "$ARGV[0]" or die "no shellfile\n";
#$symfile = "$ARGV[2]" or die "no symfile\n";
$outfile = "/home/hafiz/reaction_mechanism/CH3Br/test/rotopo_w_TS2.dat";

&define_constants;
&read_shells( $shellfile, \@shells );
#&read_symmetry( $symfile, \%symmetry );
#$sym = $symmetry{$molecule};
#print "$molecule symmetry $sym\n";
$sym = "wat"; # enforce no symmetry#

&print_entropy( $outfile, \@shells ); #Prints out the shell entropy

#--------------------------------------------------------------------------

sub define_constants {
    $pi = 3.141592654;
    $R = 8.3145;
}

sub read_shells {
    my ($file, $sref) = @_;
    open FILE, "$file" or die "Can't find $file\n";
#    print "Reading $file\n";
    RST: while (<FILE>) {
	chomp;
	@line = (split);
	($line[0] eq "1") and push @$sref, [split]; # Only take water as centre
#	printf FILE "%3d%20d\n", $$;
    }
    close FILE;
}

sub read_symmetry {
    my ($file, $symref) = @_;
    open FILE, "$file" or die "Can't find $file\n";
#    print "Reading $file\n";
    RST: while (<FILE>) {
	chomp;
	$_ or last; #Quit on first empty line 
	@line = (split);
	$$symref{$line[0]} = $line[1];
    }
    close FILE;
}

# Print out the pair data as an unwrapped matrix
sub print_entropy {
    my ($file, $sref) = @_;
    open FILE, ">$file";
    $poptot = 0;
    for $i  (0 .. $#$sref) {
	$poptot += $$sref[$i][2];
    }
    ($poptot) ? ($norm = 1/$poptot): ($norm = 0);
    $stot = 0;
    for $i  (0 .. $#$sref) {
	$frag = $$sref[$i][0];
	$shell = $$sref[$i][1];
	$pop = $$sref[$i][2];
	$prob = $pop*$norm;
	my $count0 = () = $shell =~ /0/g;
	#	$count0 > 0 and next;
	my $count1 = () = $shell =~ /1/g;
	$count0 > 0 and $count1++; # Add 1 to replace solute
	$count = $count1;
        $count or $count = 1;	# Set count to 1 if zero
	if ($sym eq "2D") {
            $s = $R*log($count);
#	    print "2D model\n";
	} elsif ($sym eq "2D2") {
            $s = $R*(log($count) - log(2));
#	    print "2D2 model\n";
        } elsif ($sym eq "2D1") {
            $s = $R*(log($count));
	} elsif ($sym eq "1D") {
	    $s = 0;
        } elsif ($sym eq "wat") {
            $s = $R*(1.5*log($count) + 0.5*log($pi) - log(2) - log(4));
	} else {
	    $s = $R*(1.5*log($count) + 0.5*log($pi) - log($sym));
#	    print "3/2 model\n";
	}
	$s < 0 and $s = 0;
	printf FILE "%3s%20s%3d%8d%8.3f%8.3f $sym\n", $frag, $shell, $count, $pop, $prob, $s;
	$stot += $s*$prob;
    }
    printf "%10.3f\n" , $stot;
    printf FILE "%10.3f\n" , $stot;
    close FILE;
}

