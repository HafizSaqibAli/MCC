#!/usr/bin/perl -w

# Calculates vibrational entropy using the 6N-6 modes of fragment (highest frequency) and 6 of whole molecule

$molecule = "$ARGV[0]";
$eigfile = "$ARGV[1]" or die "no eigfile\n";
$eigfilewm = "$ARGV[2]" or die "no eigfile\n";
$tempfile = "$ARGV[3]" or die "no temperature file\n";
$outfile = $molecule.".svib";

&define_constants;
&read_eigs( $eigfile, \@eigs );
&read_eigs( $eigfilewm, \@wmeigs );
&read_temperature( $tempfile, \%temps );
$t = $temps{$molecule} or die "No temperature for molecule $molecule\n";

&print_entropy( \@eigs, \@wmeigs, $outfile ); #Prints out the eigenvectors

#--------------------------------------------------------------------------

sub define_constants {
    $h = 6.626e-34;
    $kb = 1.3806485e-23;
    $pi = 3.141592654;
    $R = 8.3145;
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
    open FILE, ">$file";
    $kbT = $kb*$t;
# Sort eigenvalues
    @sorted_eigs = sort {$a<=>$b} @$eref;
    $stot = 0;
    for $i  (0 .. $#sorted_eigs) {
	$nu = sqrt($sorted_eigs[$i]/$kbT); # BIG TEST TO SEE IF NO FORCE HALVING
	$hnukbT = $h*$nu/($kbT*2*$pi);
	$s[$i] = $hnukbT/(exp($hnukbT)-1) - log(1-exp(-$hnukbT));
	$s[$i] *= $R;
	$scl = -$R*log($hnukbT) + $R;
	$i >=6 and $stot += $s[$i];
#	printf "%12d %10.3e%10.3e%8.3f\n", $sorted_eigs[$i],$kbT, $nu, $hnukbT;
	printf FILE "%13d%10.3f%10.3f\n", $sorted_eigs[$i], $s[$i], $scl;
    }
# Sort eigenvalues
    @wmsorted_eigs = sort {$a<=>$b} @$wmeref;
# Whole molecule entropy
    for $i  (0 .. $#wmsorted_eigs) {
# A check for argon to only read 3 eigenvalues
#	(($molecule eq "argon") and ($i <= 2)) and next;
# Calculate frequency
	$nu = sqrt($wmsorted_eigs[$i]/$kbT);
	$hnukbT = $h*$nu/($kbT*2*$pi);
	$s[$i] = $hnukbT/(exp($hnukbT)-1) - log(1-exp(-$hnukbT));
	$s[$i] *= $R;
	$scl = -$R*log($hnukbT) + $R;
	$stot += $s[$i];
#	printf "%12d %10.3e%10.3e%8.3f\n", $wmsorted_eigs[$i],$kbT, $nu, $hnukbT;
	printf FILE "%13d%10.3f%10.3f\n", $wmsorted_eigs[$i], $s[$i], $scl;
    }
    printf FILE "%10.3f\n" , $stot;
    close FILE;
}

