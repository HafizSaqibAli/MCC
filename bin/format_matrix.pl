#!/usr/bin/perl -w

# Format matrix to make easier to read

$infile = "$ARGV[0]" or die "no eigfile\n";
$outfile = "$ARGV[1]" or die "no outfile\n";
  
&read_data( $infile, \@data );
&print_entropy( $outfile, \@data ); 

#--------------------------------------------------------------------------

sub read_data {
    my ($file, $dref) = @_;
    open FILE, "$file" or die "Can't find $file\n";
    $_ = <FILE>; # Skip first line
    RST: while (<FILE>) {
	chomp;
	$_ or last; #Quit on first empty line 
	push @$dref, [split];
    }
    close FILE;
}

# Print out the pair data as an unwrapped matrix
sub print_entropy {
    my ($file, $dref) = @_;
    open FILE, ">$file";
    for $i  (0 .. $#$dref) {
	for $j (0 .. $#{$$dref[$i]}) {
	    printf FILE "%6d", $$dref[$i][$j]*1e-4;
	}
	printf FILE "\n";
    }
    close FILE;
}

