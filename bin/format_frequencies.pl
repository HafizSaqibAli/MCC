#!/usr/bin/perl -w

# Prints out the wavenumbers of the eigenvalues as sharp delta functions

$molecule = "$ARGV[0]";
$eigfile = "$ARGV[1]" or die "no eigfile\n";
$eigfilewm = "$ARGV[2]" or die "no eigfile\n";
$tempfile = "$ARGV[3]" or die "no temperature file\n";
$outfile = "$ARGV[4]" or die "no output file\n";

&define_constants;
&read_eigs( $eigfile, \@eigs );
&read_eigs( $eigfilewm, \@wmeigs );
&read_temperature( $tempfile, \%temps );
$t = $temps{$molecule} or die "No temperature for molecule $molecule\n";

&print_entropy( \@eigs, \@wmeigs, $outfile ); #Prints out the eigenvectors

#--------------------------------------------------------------------------

sub define_constants {
#    $h = 6.626e-34;
    $kb = 1.3806485e-23;
    $pi = 3.141592654;
#    $R = 8.3145;
    $c = 299792458;
}

sub read_eigs {
    my ($file, $eref) = @_;
    open FILE, "$file" or die "Can't find $file\n";
#    print "Reading $file\n";
    $_ = <FILE>; # Skip first line
    RST: while (<FILE>) {
	chomp;
	$_ or last; #Quit on first empty line 
	@line = (split);
	push @$eref, $line[0]*1.0;
#        push @$dref, $_;
    }
    close FILE;
}

sub read_temperature {
    my ($file, $tref) = @_;
    open FILE, "$file" or die "Can't find $file\n";
#    print "Reading $file\n";
    RST: while (<FILE>) {
	chomp;
	$_ or last; #Quit on first empty line 
	@line = (split);
	$$tref{$line[0]} = $line[1];
    }
    close FILE;
}

# Print out the pair data as an unwrapped matrix
sub print_entropy {
    my ($eref, $wmeref, $file) = @_;
    $kbT = $kb*$t;
# Sort eigenvalues
    @sorted_eigs = sort {$a<=>$b} @$eref;
    $binwidth = 50;
    for $i  (0 .. $#sorted_eigs) {
	$nu = sqrt($sorted_eigs[$i]/$kbT);
	$wavenumber = $nu/(2*$pi*$c*100);
	$bin = int($wavenumber/$binwidth + 0.5);
	$i >=6 and $hist[$bin]++;
	$i >=6 and printf "%12.2f%4d\n", $wavenumber, $bin;
    }
# Sort eigenvalues
    @wmsorted_eigs = sort {$a<=>$b} @$wmeref;
# Whole molecule entropy
    for $i  (0 .. $#wmsorted_eigs) {
# A check for argon to only read 3 eigenvalues
#	(($molecule eq "argon") and ($i <= 2)) and next;
# Calculate frequency
	$nu = sqrt($wmsorted_eigs[$i]/$kbT);
	$wavenumber = $nu/(2*$pi*$c*100);
	$bin = int($wavenumber/$binwidth + 0.5);
	$hist[$bin]++;
	printf "%12.2f%4d\n", $wavenumber, $bin;
    }
    open FILE, ">$file";
    for $i (0 .. $#hist) {
	if ($hist[$i]) {
	    printf FILE "%4d%4d\n", $i*$binwidth-1, 0;
	    printf FILE "%4d%4d\n", $i*$binwidth, $hist[$i];
	    printf FILE "%4d%4d\n", $i*$binwidth+1, 0;
	} else {
#	    printf FILE "%4d%4d\n", $i*$binwidth, 0;
	}
    }
    close FILE;
}
