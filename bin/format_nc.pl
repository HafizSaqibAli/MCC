#!/usr/bin/perl -w

# Print coordination shell file to have number of elements instead of 0

$infile = "$ARGV[0]" or die "enter infile\n";
$outfile = "$ARGV[1]" or die "enter outfile\n"; 

&read_shells( $infile, \%shells );
&print_shells( $outfile, \%shells ); 

#--------------------------------------------------------------------------

sub read_shells {
    my ($file, $sref) = @_;
    open FILE, "$file" or die "Can't find $file\n";
#    print "Reading $file\n";
    RST: while (<FILE>) {
	chomp;
	@line = (split);
	$shell = $line[1];
	my $count = () = $shell =~ /0/g;
	$$sref{$count} = $line[2];
#	printf FILE "%3d%20d\n", $$;
    }
    close FILE;
}

sub print_shells {
    my ($file, $sref) = @_;
# Get the normalisation
    $poptot = 0;
    foreach $i  (keys %$sref) {
	$poptot += $$sref{$i};
    }
    ($poptot) ? ($norm = 1/$poptot): ($norm = 0);
# Print shell count for range 0 to 20
    open FILE, ">$file";
    $shellmax = 20;
    for $i  (0 .. $shellmax) {
	if ($$sref{$i}) {
	    $pop = $$sref{$i};
	    $prob = $pop*$norm;
	    printf FILE "%3d %10.5f %d\n", $i, $prob, $pop;
	} else {
	    printf FILE "%3d %10.5f %d\n", $i, 0, 0;
	}
    }
    close FILE;
}

