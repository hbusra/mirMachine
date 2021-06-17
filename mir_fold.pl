#!/usr/bin/perl -w

use List::Util qw[min max];
my $goodcount = 0;
my $badcount = 0;
my $suspectcount = 0;
my $loopno = 0;

if ( $#ARGV != 4 )
{
	print "Correct usage is:\n";
	print "perl mir_fold.pl [miRNA query file] [HitTable] [BLASTdb] [long/short] [true/false]\n";
}

my ( $mirnaquery, $infile, $blastdb, $run, $sRNAseq ) = @ARGV or die "Please specify names of a file containing known miRNAs (in fasta format), BLAST hit table (generated using '-outfmt 6' when running BLAST) and blast database from which the hit table was generated.\n";

# Check for presence of UNAFold and maximum number of folds

print "---  miRNaL_fold.pl v1.0.0  ---\n";
system ("RNAfold -V > unachk.tmp");
open ( UNAREPORT, "unachk.tmp") or die "RNAfold is not available.  Please check that it is installed and in your PATH.\n";
while ($firstline = <UNAREPORT>)
{
	chomp $firstline;
#	$firstline = substr ($firstline, 21);
	print "RNAfold version $firstline detected\n";
	last;
}
unlink ("unachk.tmp");
print "miRNAs that produced a large number of hits with mir_find\n may be matching repeat elements instead of real miRNAs.\nTo save computing time, you can limit the number of hits \nto be folded for any miRNA.\n";
#print "Please enter the maximum no. of hits to be folded (press Return for no limit): ";
my $maxhits = <STDIN>;
chomp $maxhits;
if ( $maxhits )
{
	if ( $maxhits =~ m/\D+/ ) { die "That wasn't a number.  Try again!\n" }
	print "Proceeding to fold a maximum of $maxhits hits per miRNA.\n";
}
else
{
	print "Folding all hits.\n";
}

# Convert fasta file into a table of miRNAs and
# Populate hash table of miRNAs from newly generated table

open ( MIRNAS, $mirnaquery ) or die "Could not open $mirnaquery.  Is it in the right folder?\n";
open ( STARTABLE, ">".$mirnaquery.".tbl" ) or die "Could not open an output file!\n";

print STARTABLE "# miRNA ID\tsequence\n";
my (%QuerymiRNAs);
while ($line = <MIRNAS>) {
	chomp $line;
	$line =~ s/\r//;
	if ($line =~ /^>(\S+)\s*/) {
		$idkey = $1;
	}
	else {
		$QuerymiRNAs{ "$idkey" } = "$line";
		print STARTABLE "$idkey\t$line\n";
	}
}
close MIRNAS;
close STARTABLE;

my %AllmiRNAs = reverse %QuerymiRNAs;

unlink ( "$mirnaquery.tbl" );
print "miRNA data analysed.  Preparing to fold RNA sequences.\n";

# Reading in data from BLAST hit table and carrying out initial folds:

unless ( -e $infile && -f $infile && -r $infile )
{
	die "$infile cannot be accessed.  Does it exist?\n";
}

open ( IN, $infile ) or die "I don't have permission to open $infile!\n";
open ( OUT, ">".$infile.".seqtable" ) or die "Can't open an output file!\n";
open ( OUTTABLE, ">".$infile.".tbl") or die "Can't open an output file!\n";
open ( DISCARD, ">".$infile.".rejects" ) or die "Can't open an output file!\n";
open ( SUSPECT, ">".$infile.".suspect.tbl" ) or die "Can't open an output file!\n";
open ( SUSPECTOUT, ">".$infile.".suspect.seqtable" ) or die "Can't open an output file!\n";
open ( LOGFILE, ">".$infile.".log" ) or die "Can't open an output file!\n";
print DISCARD "Matches rejected:\n";
print DISCARD "Seq ID\tmiRNA\tStart\tEnd\tReason for rejection\n";
print OUTTABLE "#Structure       \tMatched                \tConserved\tMatch\tMature miRNA\tComplement\n";
print OUTTABLE "#Filename        \tSequence ID\tstart\tend\tmiRNA    \tlength\tstart\tend\tstart\tend\n";
print SUSPECT "#Structure       \tMatched                \tConserved\tMatch\tMature miRNA\tComplement\n";
print SUSPECT "#Filename        \tSequence ID\tstart\tend\tmiRNA    \tlength\tstart\tend\tstart\tend\n";

my $samecount = 0;
my $samename = "";
while ($line = <IN>) {
	chomp $line;
	$line =~ s/\r//;
	next if $line =~ /^#/;
	next if $line eq "";
#	print ".\n";

# Get the values from the table, check whether the maximum folds criterion has been exceeded;
# if not, retrieve the sequence from the BLAST database and write it to a fasta file

	my ( $qid, $sid, $percent, $allength, $mismatch, $gaps, $qstart, $qend, $sstart, $send, $evalue, $bitscore ) = split /\t/, $line;
	$loopno++;
	my $uniqueid = $loopno."_".$sid;
	if ( $maxhits ) {
		if ( $qid ne $samename ) {
			$samename = $qid;
			$samecount = 0;
		}
		$samecount++;
		if ( $samecount > $maxhits ) {
			print LOGFILE "Maximum hits for $qid exceeded, ignoring hit $uniqueid\n";
			next;
		}
	}
	my $qlength = length $QuerymiRNAs{ $qid };	
	system ("blastdbcmd -db $blastdb -entry $sid -outfmt %f -out $uniqueid.fsa" );
	print LOGFILE "Testing hit for $qid on sequence $uniqueid\n";

# Reverse complement the sequence if it is on the negative strand; otherwise just convert Ts to Us.  Either way, clip the sequence if it is too long.
# Also, truncate the header line if it is longer than 130 characters (otherwise UNAFold doesn't generate structure files)
	my $mirstart = $sstart;
	my $mirend = $send;
	if ( $sstart > $send )
	{
		print LOGFILE "This hit is on the negative strand.  Initiating reverse complement.\n";		
		my $tempseq = "";
		my $defline = "";
		my $clipstart = 0;		
		open ( NEGSTRAND, $uniqueid.".fsa" );
		while ( $line = <NEGSTRAND> )
		{
			chomp $line;
			if ( $line =~ m/^>/)
			{ 
				$defline = $line;
				if ( length($defline) > 130 )
				{
					$defline = substr ($defline, 0, 130);
				}
			}
			else 
			{
				$tempseq = $tempseq.$line;
			}
		}
		close NEGSTRAND;
		$tempseq = scalar reverse ("$tempseq");
		$tempseq =~ tr/[UT]/a/;
		$tempseq =~ tr/A/u/;
		$tempseq =~ tr/C/g/;
		$tempseq =~ tr/G/c/;
		$tempseq =~ tr/acgu/ACGU/;
		$seqlength = length ($tempseq);
		$sstart = 1 + $seqlength - $sstart;
		$send = 1 + $seqlength - $send;
		if ($seqlength > 750)
		{
			print LOGFILE "Long sequence, clipping 300bp either side of miRNA\n";
			$clipstart = $sstart - 300;
			if ($clipstart < 0)
			{
				$clipstart = 0;
			}
			$tempseq = substr $tempseq, $clipstart, 725;
			$sstart = $sstart - $clipstart;
			$send = $send -$clipstart;
		}
		open  ( POSSTRAND, ">".$uniqueid.".fsa" );		
		print POSSTRAND "$defline\n$tempseq\n";
		close POSSTRAND;			
		print LOGFILE "Reverse complement stats (id, length, hit start and end): $sid, $seqlength, $sstart, $send\n";
	}
	else
	{
		my $tempseq = "";
		my $defline = "";
		open ( DNASTRAND, $uniqueid.".fsa" );
		while ( $line = <DNASTRAND> )
		{
			chomp $line;
			if ( $line =~ m/^>/)
			{ 
				$defline = $line;
				if ( length($defline) > 130 )
				{
					$defline = substr ($defline, 0, 130);
				}
			}
			else 
			{
				$tempseq = $tempseq.$line;
			}
		}
		close DNASTRAND;
		$tempseq =~ tr/T/U/;
		if (length ($tempseq) > 750)
		{
			print LOGFILE "Long sequence, clipping 300bp either side of miRNA\n";
			$clipstart = $sstart - 300;
			if ($clipstart < 0)
			{
				$clipstart = 0;
			}
			$tempseq = substr $tempseq, $clipstart, 725;
			$sstart = $sstart - $clipstart;
			$send = $send - $clipstart;
		}
		open  ( RNASTRAND, ">".$uniqueid.".fsa" );		
		print RNASTRAND "$defline\n$tempseq\n";
		close RNASTRAND;
	}

# Adjust ends of putative mature miRNA if it is shorter than query miRNA, and note 'extended' miRNA duplex end	
	$sstart = $sstart - $qstart + 1;
	$send = $send + $qlength - $qend;
	my $exstart = $sstart - 2; 

# Run RNAfold on the fasta file
	print LOGFILE "Running RNAfold on $sid\n";
	system ("RNAfold --noPS -i $uniqueid.fsa > $uniqueid.fsa.st");

# Check that RNAfold worked
	unless (-e $uniqueid.".fsa.st") {
		print "RNAfold gave no result.  Probably a low complexity sequence. Files deleted.\n";
		print LOGFILE "RNAfold gave no result.  Probably a low complexity sequence. Files deleted.\n";
		print DISCARD "$sid\t$qid\t$sstart\t$send\tNo result from RNAFold\n";
		unlink ( <$uniqueid.*> );
		$badcount++;
		next;
	}
	
# Examine the st file for unpaired bases in the miRNA sequence
	open ( FOLD, $uniqueid.".fsa.st" ) or die "where is the st file?";
	my $linecount = 0;
	my $ssflag = 0;
	my $revstart = 0;
	my $revstartbackup = 0;
	my $revend = 0;
	my $revendbackup = 0;
	my $revlength = 0;
	$seq = "";
	$mfe=0; 
	$nucleotides="";
	while ($line = <FOLD>) {
		chomp $line;
		$linecount++;
		if ( $linecount==2 ) {
			$nucleotides = $line;
		}
		next if $linecount < 3;
		( $seq, $mfe ) = split / /, $line;
		$mfe = substr $mfe, 1, -1;
		$mfe = $mfe *1;
	}
	close FOLD;

	my @spl = split(//, $seq);
	my @stack;
	my %data;
	$linecount = 0;
	foreach my $i (@spl) { 
		$linecount++;
		if ( $i eq "." ) {
			$data{$linecount} = 0 ;
		} elsif ( $i eq "(" ) {
			push(@stack, $linecount);
		} elsif ( $i eq ")" ) {
			$value = pop(@stack);
			$data{$linecount} = $value;
			$data{$value} = $linecount;
		}
	}
	$GCcount = () = $nucleotides =~ /G|C/g;

	$linecount=0;
	foreach my $i (@spl) { 
		$linecount++;
		next if $linecount < ($exstart);
		next if $linecount > ($send);
		$ssflag++ if $data{$linecount} == 0;
		if ( $linecount == $exstart && $data{$linecount} ne 0) {
			$revend = $data{$linecount} if $revend == 0;
		}
		if ( $linecount == $sstart && $data{$linecount} ne 0) {
			$revend = $data{$linecount} +2;
		}
		if ( $linecount == $sstart - 1 && $data{$linecount} ne 0) {
			$revend = $data{$linecount} +1 if $revend == 0;
		}
		if ( $linecount == $send - 2 && $data{$linecount} ne 0) {
			$revstart = $data{$linecount};
		}
		if ( $linecount == $send - 1 && $data{$linecount} ne 0) {
			$revstart = $data{$linecount} + 1 if $revstart == 0;
		}
		if ( $linecount == $send && $data{$linecount} ne 0) {
			$revstart = $data{$linecount} + 2 if $revstart == 0;
		}
	}
	$revlength = $revend - $revstart + 1;
	if ( $ssflag > 4 )
	{
		print LOGFILE "Too many unpaired bases in the miRNA region of $sid. Files deleted.\n";
		print DISCARD "$sid\t$qid\t$sstart\t$send\t$ssflag unpaired bases in miRNA\n";
		unlink ( <$uniqueid.*> );
		$badcount++;
	}
	elsif ( $revend < 0 or $revstart < 0)
	{
		print LOGFILE "One end of the miRNA region was not base-paired, so could not locate the end of the miRNA complementary sequence.  This hit will be placed in the 'suspect' table.\n";
		if ( $revend < 0 )
		{
			$revend = $revstart + $qlength - 1;
		}
		elsif ( $revstart < 0 )
		{
			$revstart = $revend - $qlength + 1;
		}
		print SUSPECT "$uniqueid\t$sid\t$mirstart\t$mirend\t$qid\t$qlength\t$sstart\t$send\t$revstart\t$revend\n";
		open ( FASTA, $uniqueid.".fsa");
		while ($line = <FASTA>)
		{	
			chomp $line;		
			if ( $line =~ /^>/ )
			{
				print SUSPECTOUT "$uniqueid\t";
			}
			else
			{
				print SUSPECTOUT "$line";
			}
		}
		print SUSPECTOUT "\n";
		close FASTA;
		unlink ( <$uniqueid.fsa*> );
		$suspectcount++;
	}	
	elsif ( abs($revlength - $qlength) > 3 )
	{
		print LOGFILE "The miRNA complementary sequence of $sid contains breaks or a large loop. Files deleted. $revlength vs $qlength\n";
		print DISCARD "$sid\t$qid\t$sstart\t$send\tmiRNA complementary sequence is broken \n";
		unlink ( <$uniqueid.*> );
		$badcount++;
	}
	elsif ( $ssflag == 0 )
	{
		print LOGFILE "The putative miRNA sequence of $sid is perfectly base-paired, so it is more likely to be an inverted repeat or siRNA.  This hit will be placed in the 'suspect' table.\n";
		print SUSPECT "$uniqueid\t$sid\t$mirstart\t$mirend\t$qid\t$qlength\t$sstart\t$send\t$revstart\t$revend\n";
		open ( FASTA, $uniqueid.".fsa");
		while ($line = <FASTA>)
		{	
			chomp $line;		
			if ( $line =~ /^>/ )
			{
				print SUSPECTOUT "$uniqueid\t";
			}
			else
			{
				print SUSPECTOUT "$line";
			}
		}
		print SUSPECTOUT "\n";
		close FASTA;
		unlink ( <$uniqueid.fsa*> );
		unlink ( "$uniqueid.fsa_1.ss" );
		unlink ( "$uniqueid.fsa_1.pdf" );
		$suspectcount++;
	}	
	else
	{
		print LOGFILE "The secondary structure of $uniqueid passes initial analysis, writing to output table and fasta file.\n";		
		print OUTTABLE "$uniqueid\t$sid\t$mirstart\t$mirend\t$qid\t$qlength\t$sstart\t$send\t$revstart\t$revend\n";		
		open ( FASTA, $uniqueid.".fsa");
		while ($line = <FASTA>)
		{	
			chomp $line;		
			if ( $line =~ /^>/ )
			{
				print OUT "$uniqueid\t";
			}
			else
			{
				print OUT "$line";
			}
		}
		print OUT "\n";
		close FASTA;
		unlink ( <$uniqueid.fsa*> );
		$goodcount++;
	}
}
close IN;
close OUT;
close OUTTABLE;
#close SUSPECT;
close SUSPECTOUT;

print LOGFILE "From the input table, $goodcount sequence(s) gave folds that could contain a miRNA.\n$badcount sequence(s) were rejected after folding.\n$suspectcount sequence(s) are more likely to be repeats or siRNAs.\n\n";
close LOGFILE;

# Conduct further analysis on hairpin regions of all good hits

open ( GOODTABLE, $infile.".tbl" );
open ( RESULTS, ">".$infile.".hairpins.tbl" ) or die "Can't open an output file!\n";
open ( HAIRPINS, ">".$infile.".hairpins.fsa" ) or die "Can't open an output file!\n";
open ( LOGFILE2, ">".$infile.".hairpins.log" ) or die "Can't open an output file!\n";
print RESULTS "Unique\tNew miRNA\t\t\tConserved miRNA\t\t\tSequence\tMature\tMature\tmiRNA*\tmiRNA*\tmiRNA*\tHairpin\tPre-miRNA stats\n";
print RESULTS "Hit ID\tID\tSequence\tLength\tID\tSequence\tMismatch\tID\tStart\tEnd\tStart\tEnd\tSequence\tlocation\tlength\tMFE\tGC%\tMFEI\tstart\tsequence\n";
my $goodhairpins = 0;

while ( $line = <GOODTABLE> ) {
	chomp $line;
	next if $line =~ /^#/;
	my ( $uniqueid, $sseqid, $mirstart, $mirend, $mirnaid, $length, $sstart, $send, $revstart, $revend ) = split /\t/, $line;
	my ( $armflag, $seq ) = "";
	my ( $hairpinstart, $hairpinend, $hairpinlength ) = 0;
	print LOGFILE2 "Analysing hit $uniqueid for $mirnaid\n";	
	
# Determine which part of the sequence corresponds to putative hairpin and retrieve it from the sequence table
	if ( $sstart > $revstart ) {
		$armflag = "3'";
		$hairpinstart = $revstart - 20;
		$hairpinend = $send + 20;
	} elsif ( $sstart < $revstart ) {
		$armflag = "5'";
		$hairpinstart = $sstart - 20;
		$hairpinend = $revend + 20;
	} else {
		print LOGFILE2 "The miRNA co-ordinates for $uniqueid don't make sense!  Skipping it.\n";
		next;
	}
	
	open ( GOODSEQS, $infile.".seqtable" );
	while ( $line = <GOODSEQS> ) {
		my ( $fseqid, $fseq) = split /\t/, $line;
		if  ( $fseqid =~ m/$uniqueid/ ) {
			$seq = $fseq;	
			last;
		}
	}
	close GOODSEQS;
	if ( $hairpinend - $hairpinstart < (2*$length) + 40 ) {
		print LOGFILE2 "The miRNA region of $sseqid goes round the head of the hairpin; discarded.\n";
		print DISCARD "$sseqid\t$mirnaid\t$sstart\t$send\tmiRNA goes round head of hairpin\n";
		next;
	}
 	if ( $hairpinstart < 1 )
	{
		$hairpinstart = 1;
	}
	if ( $hairpinend > length $seq )
	{
		$hairpinend = length $seq;
	}
	$hairpinlength = $hairpinend - $hairpinstart + 1;
	
	$hairpinseq = substr $seq, ($hairpinstart-1), $hairpinlength;
	$matseq = substr $seq, ($sstart-1), $length;	
	$matstart = $sstart - $hairpinstart+1;
	$matend = $send - $hairpinstart+1;
	
	print LOGFILE2 "Hairpin stats (start, end, length):\n";
	print LOGFILE2 "$hairpinstart\t$hairpinend\t$hairpinlength\n";
	print LOGFILE2 "Running RNAfold on the hairpin sequence... ";
	open ( PFASTA, ">".$uniqueid.".hairpin.fsa" ) or die "Can't open an output file!\n";
	print PFASTA ">$sseqid\t$mirnaid\n";
	print PFASTA "$hairpinseq\n";
	close PFASTA;


# Run RNAfold on the fasta file
	system ("RNAfold --noPS -i $uniqueid.hairpin.fsa > $uniqueid.hairpin.fsa.st");

# Check that RNAfold worked
	unless (-e $uniqueid.".hairpin.fsa.st") {
		print "RNAfold gave no result.  Probably a low complexity sequence. Files deleted.\n";
		print LOGFILE2 "RNAfold gave no result.  Probably a low complexity sequence. Files deleted.\n";
		unlink ( <$uniqueid.*> );
		next;
	}

# Examine the st file for unpaired bases in the miRNA sequence
	open ( HAIRPINFOLD, $uniqueid.".hairpin.fsa.st" ) or die "where is the st file?";

	$linecount = 0;
	$seq = "";
	$mfe=0; 
	$nucleotides="";
	while ($line = <HAIRPINFOLD>) {
		chomp $line;
		$linecount++;
		if ( $linecount==2 ) {
			$nucleotides = $line;
		}
		next if $linecount < 3;
		( $seq, $mfe ) = split / /, $line;
		$mfe = substr $mfe, 1, -1;
		$mfe = $mfe * 1;
	}
	close HAIRPINFOLD;

	@spl = split(//, $seq);
	@stack = ();
	%data=();
	$linecount = 0;
	foreach my $i (@spl) { 
		$linecount++;
		if ( $i eq "." ) {
			$data{$linecount} = 0 ;
		} elsif ( $i eq "(" ) {
			push(@stack, $linecount);
		} elsif ( $i eq ")" ) {
			$value = pop(@stack);
			$data{$linecount} = $value;
			$data{$value} = $linecount;
		}
	}
	$GCcount = () = $nucleotides =~ /G|C/g;

	my $starstart = 0;
	my $starend = 0;
	my $starendbackup = 0;
	my $starstartbackup = 0;
	my $starseq = "Not defined";
	my $ssflag = 0;
	my $dicerflag = 0;
	my $ssflag2 = 0;
	$linecount=0;
	foreach my $i (@spl) { 
		$linecount++;
		if ($linecount == $matstart -2 && $data{$linecount} ne 0) {
			$starend = $data{$linecount} if $starend == 0;
		} elsif ($linecount == $matstart -1 && $data{$linecount} ne 0) {
			$starend = $data{$linecount} + 1 if $starend == 0;
		} elsif ($linecount == $matstart && $data{$linecount} ne 0) {
			$starend = $data{$linecount} + 2;
		} elsif ($linecount == $matend && $data{$linecount} ne 0) {
			$starstart = $data{$linecount} + 2 if $starstart == 0;
		} elsif ($linecount == $matend -1 && $data{$linecount} ne 0) {
			$starstart = $data{$linecount} + 1 if $starstart == 0;
		} elsif ($linecount == $matend-2 && $data{$linecount} ne 0) {
			$starstart = $data{$linecount};
		}
	}

	$start = min($starstart, $matstart, $starend, $matend);
	$end = max($starstart, $matstart, $starend, $matend);
	$multiloopflag=0;
	$headflag = 0;
	$check = 0;
	$headsize = 0;
	my $k = 0;
	$linecount=0;
	foreach my $i (@spl) { 
		$linecount++;
		$temp = $data{$linecount} if $linecount == $start;
		if ($linecount >= $matstart && $linecount <= $matend) {
			$ssflag++ if $data{$linecount} == 0;
		}	
		if ($linecount >= $starstart && $linecount <= $starend) {
			$ssflag2++ if $data{$linecount} == 0;
		}
		if ($linecount == $starstart or $linecount == $matstart) {
			$dicerflag = 1 if $data{$linecount} == 0;
		}
		if ($linecount > $start && $linecount <= $end) {
			if ($data{$linecount} < $temp) {
				$temp = $data{$linecount} if $data{$linecount} ne 0;
			} else {
				$multiloopflag =1;
			}
			
			if ($data{$linecount} == 0) {
				$headsize++;
				$check = $linecount -1 if $headsize == 1;
			} else {
				$headsize = 0;
				if ($check == $data{$linecount}) {
					for ($k = $linecount - $headsize; $k < $linecount; $k++){
						$headflag =1 if (($k > $matstart && $k <= $matend) or ($k > $starstart && $k <= $starend));
					}
				}
			}
		}
	}
	if ($starend > 0 && $starstart > 0) {
		$starseq = substr $hairpinseq, ($starstart-1), ($starend-$starstart+1);
	}
	my $tempseq = $starseq;
	$tempseq =~ tr/[UT]/T/;

	if ( $ssflag > 4 ) {
		print LOGFILE2 "Too many unpaired bases in the miRNA region of $uniqueid. Files deleted.\n";
		print DISCARD "$sseqid\t$mirnaid\t$sstart\t$send\t$ssflag unpaired bases in miRNA\n";
		unlink ( <$uniqueid.*> );
		next;
	} elsif ( $ssflag2 > 6 ) {
		print LOGFILE2 "Too many unpaired bases in the miRNA* region of $uniqueid. Files deleted.\n";
		print DISCARD "$sseqid\t$mirnaid\t$sstart\t$send\t$ssflag2 unpaired bases in miRNA*\n";
		unlink ( <$uniqueid.*> );
		next;
	} elsif ( $dicerflag > 0 ) {
		print LOGFILE2 "Dicer-like enzyme cut regions are not accesible in the miRNA region of $uniqueid. Files deleted.\n";
		print DISCARD "$sseqid\t$mirnaid\t$sstart\t$send\tmismatches at the start of miRNA or miRNA* sequences\n";
		unlink ( <$uniqueid.*> );
		next;
	} elsif ( $multiloopflag > 0 ) {
		print LOGFILE2 "Multiple loops occur in the miRNA region of $uniqueid. Files deleted.\n";
		print DISCARD "$sseqid\t$mirnaid\t$sstart\t$send\tmultiple loops between the miRNA and miRNA* sequences\n";
		unlink ( <$uniqueid.*> );
		next;
	} elsif ( $headflag > 0 ) {
		print LOGFILE2 "The miRNA region of $uniqueid is in the head region. Files deleted.\n";
		print DISCARD "$sseqid\t$mirnaid\t$sstart\t$send\tmiRNA is in the head region\n";
		unlink ( <$uniqueid.*> );
		next;
	}

# Process data and print results files
	my $mirnafam = substr ($mirnaid, 3);
	my $conseq = $QuerymiRNAs{"$mirnaid"};
	my $GCcomp = 100*$GCcount/$hairpinlength;
	my $amfe = (0-100*$mfe/$hairpinlength);
	if ( $GCcomp == 0 ) {
		print LOGFILE2 "GC content appears to be 0.  That's weird.\n";
		print DISCARD "$sseqid\t$mirnaid\t$sstart\t$send\tGC content outside acceptable range\n";
		unlink ( <$uniqueid.*> );
		next;
	}
	my $mfei = $amfe/$GCcomp;
	my $mism = 0;
	for ($loopcount = 0; $loopcount < $length; $loopcount++) {
		$mism++ if substr ($conseq, $loopcount, 1) ne substr ($matseq, $loopcount, 1);
	}
	if ( $GCcomp < 24 or $GCcomp > 71 ) {
		print LOGFILE2 "The hairpin for $uniqueid has too low or high GC content.  Deleting files.\n";
		print DISCARD "$sseqid\t$mirnaid\t$sstart\t$send\tGC content outside acceptable range\n";
		unlink ( <$uniqueid.*> );
	} elsif ( $mfei < 0.67 ) {
		print LOGFILE2 "The hairpin for $uniqueid has MFEI < 0.67.  Deleting files.\n";
		print DISCARD "$sseqid\t$mirnaid\t$sstart\t$send\tMFEI outside acceptable range\n";
		unlink ( <$uniqueid.*> );
	} elsif ($sRNAseq eq "true" && ! exists $AllmiRNAs {$tempseq}) {
 		print LOGFILE2 "The miRNA star expression evidence could not be found; $uniqueid moved to suspect tbl.\n";
		print SUSPECT "$uniqueid\t$sseqid\t$mirstart\t$mirend\t$mirnaid\t$length\t$sstart\t$send\t$revstart\t$revend\n";
 		unlink ( <$uniqueid.*> );
 		next;
	} else {
		print LOGFILE2 "The hairpin passes analysis, printing to results table\n";
		print RESULTS "$uniqueid\tcandidate$mirnafam\t$matseq\t$length\t$mirnaid\t$conseq\t$mism\t$sseqid\t$mirstart\t$mirend\t$matstart\t$matend\t$starstart\t$starend\t$starseq\t$armflag\t$hairpinlength\t$mfe\t$GCcomp\t$mfei\t$hairpinstart\t$hairpinseq\n";
		print HAIRPINS ">$sseqid\t$mirnaid\n$hairpinseq\n";
		unlink ( <$uniqueid.hairpin.fsa*> );
		$goodhairpins++
	}

}
close SUSPECT;
close RESULTS;
close HAIRPINS;
close GOODTABLE;
close DISCARD;
close LOGFILE2;
unlink ( "$infile.seqtable" );

open ( RESULTS, $infile.".hairpins.tbl" ) or die "Can't open an output file!\n";
open ( OUTTBL, ">".$infile.".hairpins.tbl.out.tbl" ) or die "Can't open an output file!\n";
my $k=0;
my @blacklist = ('Unique', 'Hit', 'Hit ID');
my (%stock, %stocknames);
while ( $line = <RESULTS> )
{
	chomp $line;
	my @lineelems = split(/\t/, $line);

	if ($line ne "" && substr($line, 0, 2) ne "\t" && !(grep( /^$lineelems[0]$/, @blacklist ))){
		my $n = scalar(@lineelems);
		my $temp = $line;
		if ($n < 22){
			$temp .= $lines[$lineID+1];
			my @tmparr = split(/\t/, $lines[$lineID+1]);
			$n += scalar(@tmparr);
			if ($n < 22){
				$temp .= $lines[$lineID+2];
			}
		}
		my @store = split(/\t/, $temp);
		$k+=1;
		my $loc = substr($store[15], 0, 1);
		my $mir="";
		if ($sRNAseq eq "true") {
			$mir = $store[4].= "-".$loc."p";
		} else {
			$mirID = (split(/\./, (split(/-/, $store[1]))[1]))[0];
			$mir = $mirID;
			$mir =~ s/[^\d]//g;
			$homolog_loc = (split(/-/, $store[4]))[-1];
			if (grep(/^$homolog_loc$/, ("3p", "5p")) && $loc ne substr($homolog_loc, 0, 1)){
				$mir = $mir;
			} else {
				$mir .= "-".$loc."p";
			}
			$mir="miR".$mir;
		}
		unless ( exists $stock{"$mir, $store[2], $store[21]"} ){ $stocknames{$mir}++; $stock{"$mir, $store[5], $store[21]"}++; }
		$mir.="-".$stocknames{$mir};
		print OUTTBL "$mir\t$store[2]\t$store[21]\t$k\t$store[4]\t$store[14]\t$store[7]\t$store[8]\t$store[9]\n";

	}
}

close RESULTS;
close OUTTBL;



if ( $run eq "long" ){
	# Carry out the same hairpin analysis on suspect hits
	open ( SUSPECTTABLE, $infile.".suspect.tbl" );
	open ( SUSPECTRESULTS, ">".$infile.".suspecthairpins.tbl" ) or die "Can't open an output file!\n";
	open ( SUSPECTHAIRPINS, ">".$infile.".suspecthairpins.fsa" ) or die "Can't open an output file!\n";
	open ( LOGFILE2, ">>".$infile.".hairpins.log" ) or die "Can't open an output file!\n";
	print LOGFILE2 "Moving on to examine suspect hits from initial folding analysis...\n";
	print SUSPECTRESULTS "Unique\tNew miRNA\t\t\tConserved miRNA\t\t\tSequence\tMature\tMature\tmiRNA*\tmiRNA*\tmiRNA*\tHairpin\tPre-miRNA stats\n";
	print SUSPECTRESULTS "Hit ID\tID\tSequence\tLength\tID\tSequence\tMismatch\tID\tStart\tEnd\tStart\tEnd\tSequence\tlocation\tlength\tMFE\tGC%\tMFEI\tstart\tsequence\n";
	open ( DISCARD, ">>".$infile.".rejects" ) or die "Can't open an output file!\n";
	my $suspecthairpins = 0;

	while ( $line = <SUSPECTTABLE> ) {
		chomp $line;
		next if $line =~ /^#/;
		my ( $uniqueid, $sseqid, $mirstart, $mirend, $mirnaid, $length, $sstart, $send, $revstart, $revend ) = split /\t/, $line;
		my ( $armflag, $seq ) = "";
		my ( $hairpinstart, $hairpinend, $hairpinlength ) = 0;
		print LOGFILE2 "Analysing hit $uniqueid for $mirnaid\n";	
	
	# Determine which part of the sequence corresponds to putative hairpin and retrieve it from the sequence table
		if ( $sstart > $revstart ) {
			$armflag = "3'";
			$hairpinstart = $revstart - 20;
			$hairpinend = $send + 20;
		} elsif ( $sstart < $revstart ) {
			$armflag = "5'";
			$hairpinstart = $sstart - 20;
			$hairpinend = $revend + 20;
		} else {
			print LOGFILE2 "The miRNA co-ordinates for $uniqueid don't make sense!  Skipping it.\n";
			next;
		}
	
		open ( SUSPECTSEQS, $infile.".suspect.seqtable" );
		while ( $line = <SUSPECTSEQS> ) {
			my ( $fseqid, $fseq) = split /\t/, $line;
			if  ( $fseqid =~ m/$uniqueid/ ) {
				$seq = $fseq;	
				last;
			}
		}
		close SUSPECTSEQS;
		if ( $hairpinend - $hairpinstart < (2*$length) + 40 ) {
			print LOGFILE2 "The miRNA region of $sseqid goes round the head of the hairpin; discarded.\n";
			print DISCARD "$sseqid\t$mirnaid\t$sstart\t$send\tmiRNA goes round head of hairpin\n";
			next;
		}
		if ( $hairpinstart < 1 )
		{
			$hairpinstart = 1;
		}
		if ( $hairpinend > length $seq )
		{
			$hairpinend = length $seq;
		}
		$hairpinlength = $hairpinend - $hairpinstart + 1;
	
		$hairpinseq = substr $seq, ($hairpinstart-1), $hairpinlength;
		$matseq = substr $seq, ($sstart-1), $length;	
		$matstart = $sstart - $hairpinstart+1;
		$matend = $send - $hairpinstart+1;
		print LOGFILE2 "Hairpin stats (start, end, length):\n";
		print LOGFILE2 "$hairpinstart\t$hairpinend\t$hairpinlength\n";
		print LOGFILE2 "Running UNAFold.pl on the hairpin sequence... ";
		open ( PFASTA, ">".$uniqueid.".hairpin.fsa" ) or die "Can't open an output file!\n";
		print PFASTA ">$sseqid\t$mirnaid\n";
		print PFASTA "$hairpinseq\n";
		close PFASTA;

	# Run RNAfold on the fasta file
		system ("RNAfold --noPS -i $uniqueid.hairpin.fsa > $uniqueid.hairpin.fsa.st");

	# Check that RNAfold worked
		unless (-e $uniqueid.".hairpin.fsa.st") {
			print "RNAfold gave no result.  Probably a low complexity sequence. Files deleted.\n";
			print LOGFILE2 "RNAfold gave no result.  Probably a low complexity sequence. Files deleted.\n";
			unlink ( <$uniqueid.*> );
			next;
		}

	# Examine the st file for unpaired bases in the miRNA sequence
		open ( HAIRPINFOLD, $uniqueid.".hairpin.fsa.st" ) or die "where is the st file?";

		$linecount = 0;
		$seq = "";
		$mfe=0; 
		$nucleotides="";
		while ($line = <HAIRPINFOLD>) {
			chomp $line;
			$linecount++;
			if ( $linecount==2 ) {
				$nucleotides = $line;
			}
			next if $linecount < 3;
			( $seq, $mfe ) = split / /, $line;
			$mfe = substr $mfe, 1, -1;
			$mfe = $mfe * 1;
		}
		close HAIRPINFOLD;

		@spl = split(//, $seq);
		@stack = ();
		%data=();
		$linecount = 0;
		foreach my $i (@spl) { 
			$linecount++;
			if ( $i eq "." ) {
				$data{$linecount} = 0 ;
			} elsif ( $i eq "(" ) {
				push(@stack, $linecount);
			} elsif ( $i eq ")" ) {
				$value = pop(@stack);
				$data{$linecount} = $value;
				$data{$value} = $linecount;
			}
		}
		$GCcount = () = $nucleotides =~ /G|C/g;

		my $starstart = 0;
		my $starend = 0;
		my $starendbackup = 0;
		my $starstartbackup = 0;
		my $starseq = "Not defined";
		my $ssflag = 0;
		my $dicerflag = 0;
		my $ssflag2 = 0;
		$linecount=0;
		foreach my $i (@spl) { 
			$linecount++;
			if ($linecount == $matstart -2 && $data{$linecount} ne 0) {
				$starend = $data{$linecount} if $starend == 0;
			} elsif ($linecount == $matstart -1 && $data{$linecount} ne 0) {
				$starend = $data{$linecount} + 1 if $starend == 0;
			} elsif ($linecount == $matstart && $data{$linecount} ne 0) {
				$starend = $data{$linecount} + 2;
			} elsif ($linecount == $matend && $data{$linecount} ne 0) {
				$starstart = $data{$linecount} + 2 if $starstart == 0;
			} elsif ($linecount == $matend -1 && $data{$linecount} ne 0) {
				$starstart = $data{$linecount} + 1 if $starstart == 0;
			} elsif ($linecount == $matend-2 && $data{$linecount} ne 0) {
				$starstart = $data{$linecount};
			}
		}

		$start = min($starstart, $matstart, $starend, $matend);
		$end = max($starstart, $matstart, $starend, $matend);
		$multiloopflag=0;
		$headflag = 0;
		$check = 0;
		$headsize = 0;
		my $k = 0;
		$linecount=0;
		foreach my $i (@spl) { 
			$linecount++;
			$temp = $data{$linecount} if $linecount == $start;
			if ($linecount >= $matstart && $linecount <= $matend) {
				$ssflag++ if $data{$linecount} == 0;
			}	
			if ($linecount >= $starstart && $linecount <= $starend) {
				$ssflag2++ if $data{$linecount} == 0;
			}
			if ($linecount == $starstart or $linecount == $matstart) {
				$dicerflag = 1 if $data{$linecount} == 0;
			}
			if ($linecount > $start && $linecount <= $end) {
				if ($data{$linecount} < $temp) {
					$temp = $data{$linecount} if $data{$linecount} ne 0;
				} else {
					$multiloopflag =1;
				}
			
				if ($data{$linecount} == 0) {
					$headsize++;
					$check = $linecount -1 if $headsize == 1;
				} else {
					$headsize = 0;
					if ($check == $data{$linecount}) {
						for ($k = $linecount - $headsize; $k < $linecount; $k++){
							$headflag =1 if (($k > $matstart && $k <= $matend) or ($k > $starstart && $k <= $starend));
						}
					}
				}
			}
		}
		if ($starend > 0 && $starstart > 0) {
			$starseq = substr $hairpinseq, ($starstart-1), ($starend-$starstart+1);
		}
		if ( $ssflag > 4 ) {
			print LOGFILE2 "Too many unpaired bases in the miRNA region of $uniqueid. Files deleted.\n";
			print DISCARD "$sseqid\t$mirnaid\t$sstart\t$send\t$ssflag unpaired bases in miRNA\n";
			unlink ( <$uniqueid.*> );
			next;
		} elsif ( $ssflag2 > 6 ) {
			print LOGFILE2 "Too many unpaired bases in the miRNA* region of $uniqueid. Files deleted.\n";
			print DISCARD "$sseqid\t$mirnaid\t$sstart\t$send\t$ssflag2 unpaired bases in miRNA*\n";
			unlink ( <$uniqueid.*> );
			next;
		} elsif ( $dicerflag > 0 ) {
			print LOGFILE2 "Dicer-like enzyme cut regions are not accesible in the miRNA region of $uniqueid. Files deleted.\n";
			print DISCARD "$sseqid\t$mirnaid\t$sstart\t$send\tmismatches at the start of miRNA or miRNA* sequences\n";
			unlink ( <$uniqueid.*> );
			next;
		} elsif ( $multiloopflag > 0 ) {
			print LOGFILE2 "Multiple loops occur in the miRNA region of $uniqueid. Files deleted.\n";
			print DISCARD "$sseqid\t$mirnaid\t$sstart\t$send\tmultiple loops between the miRNA and miRNA* sequences\n";
			unlink ( <$uniqueid.*> );
			next;
		} elsif ( $headflag > 0 ) {
			print LOGFILE2 "The miRNA region of $uniqueid is in the head region. Files deleted.\n";
			print DISCARD "$sseqid\t$mirnaid\t$sstart\t$send\tmiRNA is in the head region\n";
			unlink ( <$uniqueid.*> );
			next;
		}

	# Process data and print results files
		my $mirnafam = substr ($mirnaid, 3);
		my $conseq = $QuerymiRNAs{"$mirnaid"};
		my $GCcomp = 100*$GCcount/$hairpinlength;
		my $amfe = (0-100*$mfe/$hairpinlength);
		if ( $GCcomp == 0 ) {
			print LOGFILE2 "GC content appears to be 0.  That's weird.\n";
			print DISCARD "$sseqid\t$mirnaid\t$sstart\t$send\tGC content outside acceptable range\n";
			unlink ( <$uniqueid.*> );
			next;
		}
		my $mfei = $amfe/$GCcomp;
		my $mism = 0;
		for ($loopcount = 0; $loopcount < $length; $loopcount++) {
			$mism++ if substr ($conseq, $loopcount, 1) ne substr ($matseq, $loopcount, 1);
		}
		if ( $GCcomp < 24 or $GCcomp > 71 ) {
			print LOGFILE2 "The hairpin for $uniqueid has too low or high GC content.  Deleting files.\n";
			print DISCARD "$sseqid\t$mirnaid\t$sstart\t$send\tGC content outside acceptable range\n";
			unlink ( <$uniqueid.*> );
		} elsif ( $mfei < 0.67 ) {
			print LOGFILE2 "The hairpin for $uniqueid has MFEI < 0.67.  Deleting files.\n";
			print DISCARD "$sseqid\t$mirnaid\t$sstart\t$send\tMFEI outside acceptable range\n";
			unlink ( <$uniqueid.*> );
		} else {
			print LOGFILE2 "The hairpin passes analysis, printing to results table\n";
			print SUSPECTRESULTS "$uniqueid\tcandidate$mirnafam\t$matseq\t$length\t$mirnaid\t$conseq\t$mism\t$sseqid\t$mirstart\t$mirend\t$matstart\t$matend\t$starstart\t$starend\t$starseq\t$armflag\t$hairpinlength\t$mfe\t$GCcomp\t$mfei\t$hairpinstart\t$hairpinseq\n";
			print SUSPECTHAIRPINS ">$sseqid\t$mirnaid\n$hairpinseq\n";
			unlink ( <$uniqueid.hairpin.fsa*> );
			$goodhairpins++
		}
	}
	close SUSPECTRESULTS;
	close SUSPECTHAIRPINS;
	close SUSPECTTABLE;
	close DISCARD;
	close LOGFILE2;

	open ( SUSPECTRESULTS, $infile.".suspecthairpins.tbl" ) or die "Can't open an output file!\n";
	open ( OUTTBL, ">".$infile.".suspecthairpins.tbl.out.tbl" ) or die "Can't open an output file!\n";
	my $k=0;
	my @blacklist = ('Unique', 'Hit', 'Hit ID');
	my (%stock, %stocknames);
	while ( $line = <SUSPECTRESULTS> )
	{
		chomp $line;
		my @lineelems = split(/\t/, $line);

		if ($line ne "" && substr($line, 0, 2) ne "\t" && !(grep( /^$lineelems[0]$/, @blacklist ))){
			my $n = scalar(@lineelems);
			my $temp = $line;
			if ($n < 22){
				$temp .= $lines[$lineID+1];
				my @tmparr = split(/\t/, $lines[$lineID+1]);
				$n += scalar(@tmparr);
				if ($n < 22){
					$temp .= $lines[$lineID+2];
				}
			}
			my @store = split(/\t/, $temp);
			my $loc = substr($store[15], 0, 1);
			my $mir="";
			if ($sRNAseq eq "true") {
				$mir = $store[4].= "-".$loc."p";
			} else {
				$mirID = (split(/\./, (split(/-/, $store[1]))[1]))[0];
				$mir = $mirID;
				$mir =~ s/[^\d]//g;
				$homolog_loc = (split(/-/, $store[4]))[-1];
				if (grep(/^$homolog_loc$/, ("3p", "5p")) && $loc ne substr($homolog_loc, 0, 1)){
					$mir = $mir;
				} else {
					$mir .= "-".$loc."p";
				}
				$mir="miR".$mir;
			}
			unless ( exists $stock{"$mir, $store[5], $store[21]"} ){
				$stocknames{$mir}++; 
				$stock{"$mir, $store[2], $store[21]"}++;
			}
			$mir.="-".$stocknames{$mir};
			print OUTTBL "$mir\t$store[2]\t$store[21]\t$k\t$store[4]\t$store[14]\t$store[7]\t$store[8]\t$store[9]\n";

		}
	}

	close SUSPECTRESULTS;
	close OUTTBL;
}
unlink ( "$infile.suspect.seqtable" );

