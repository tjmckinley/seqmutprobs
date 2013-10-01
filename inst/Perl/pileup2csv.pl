#!/usr/bin/perl

print "NAME\tPOS\tCOV\tREF\tA\tC\tG\tT\tD\n";

my %nuc;

while(<>) {
	@list = split /\t/, $_;
	my $pileupnucleotides = $list[4];
	$nuc{"A"} = 0;
    $nuc{"C"} = 0;
    $nuc{"G"} = 0;
    $nuc{"T"} = 0;
    $nuc{"*"} = 0;
	while (length $pileupnucleotides) {
		my @matches = ($pileupnucleotides  =~ m/\.|,|\*|\+[0-9]+[ACGTNacgtn]+|-[0-9]+[ACGTNacgtn]+|[ACGTNacgtn]/g);
		my %ins = ();
		$pileupnucleotides = "";
		foreach $m (@matches) {
			if ($m eq "." || $m eq ",") {
				$nuc{$list[2]}++;
			} elsif ($m =~ m/\+/) {
				#print "insertion here: " . $m . "\n";
				my @number = ($m =~ m/[0-9]+/g);
				my @sequence = ($m =~ m/[ACGTNacgtn]+/g);
				$pileupnucleotides .= substr $sequence[0], int($number[0]);
			} elsif ($m =~ m/\-/) {
				#print "deletion here: " . $m . "\n";
				my @number = ($m =~ m/[0-9]+/g);
				my @sequence = ($m =~ m/[ACGTNacgtn]+/g);
				$pileupnucleotides .= substr $sequence[0], int($number[0]);
			} else {
				$nuc{uc($m)}++;
			}
		}
	}
	print $list[0] . "\t" . $list[1] . "\t" . $list[3] . "\t" . $list[2] . "\t" . $nuc{"A"} . "\t" . $nuc{"C"} . "\t" . $nuc{"G"} . "\t" . $nuc{"T"} . "\t" . $nuc{"*"} . "\n" ;
	
}
