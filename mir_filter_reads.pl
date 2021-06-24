#!/usr/bin/perl

use warnings;
use strict;


my ($mirnaquery, $minlen, $maxlen, $rpm_cutoff) = @ARGV or die "Please specify the name of a file containing miRNAs (in fasta format), minimum length for a miRNA, maximum length of a miRNA, and the RPM threshold\n";
my $total_counts = 0;
open ( MIRNAS, $mirnaquery ) or die "Could not open $mirnaquery.  Is it in the right folder?\n";
my (%QuerymiRNAs);
while (my $line = <MIRNAS>) {
	chomp $line;
	if ($line =~ /^>/){
		next;
	}
	else{
		$QuerymiRNAs{$line}++;
		$total_counts++;
	}
}
close MIRNAS;


open ( NONREDUNDANT, ">".$mirnaquery.".nr.fasta" ) or die "Could not open an output file!\n";
open ( FILTERED, ">".$mirnaquery.".filtered.fasta" ) or die "Could not open an output file!\n";

my $count = 0;
my $rpm;
foreach my $seq(sort keys %QuerymiRNAs) {
	$rpm = $QuerymiRNAs{$seq}*1000000/$total_counts;
	if ($rpm >= $rpm_cutoff && length($seq)>=$minlen && length($seq)<=$maxlen) { 
		$rpm = sprintf('%d', $rpm); 
		print FILTERED ">read$count"."_R$rpm\n$seq\n";
	}
	print NONREDUNDANT ">read$count"."_R$rpm\n$seq\n";
	
	
	$count++;
}
close NONREDUNDANT;
close FILTERED;

