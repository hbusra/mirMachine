#!/bin/perl -w


my ($mirnaquery, $minlen, $maxlen, $rpm_cutoff) = @ARGV or die "Please specify the name of a file containing miRNAs (in fasta format), minimum length for a miRNA, maximum length of a miRNA, and the RPM threshold\n";
my $total_count = 0;
open ( MIRNAS, $mirnaquery ) or die "Could not open $mirnaquery.  Is it in the right folder?\n";
my (%QuerymiRNAs);
while ($line = <MIRNAS>) {
	chomp $line;
	$line =~ s/\r//;
	if ($line =~ /^>(\S+)\s*/) {
		my $idkey = $1;
	}
	else {
		my $seqval = $2;
		$QuerymiRNAs{$seqval}++;
		$total_counts++;
	}
}
close MIRNAS;


open ( NONREDUNDANT, ">".$mirnaquery.".nr.fasta" ) or die "Could not open an output file!\n";
open ( FILTERED, ">".$mirnaquery.".filtered.fasta" ) or die "Could not open an output file!\n";

my $count = 0;
foreach my $seq(sort keys %QuerymiRNAs) {
	$rpm = $QuerymiRNAs{$seq}*1000000/$total_counts;
	print NONREDUNDANT ">read$count_R$rpm\n$seq\n";
	if ($rpm >= $rpm_cutoff && length($seq)>=$minlen && length($seq)<=$maxlen) { print FILTERED ">read$count_R$rpm\n$seq\n"; }
	$count++;
}
close NONREDUNDANT;
close FILTERED;
