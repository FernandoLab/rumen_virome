#!/usr/bin/perl -w

use strict;
use Getopt::Long;

#Command line parameters:
my $rrna = "";

#Setup the command line options using Getopt:Long
my $commandline = GetOptions("rrna:s", \$rrna);
			
if (!$commandline || $rrna eq "" ) {
	print STDERR "\nUsage: $0 {-rrna } \n\n";
	exit;
}

open (my $FILE, "$rrna") or die "Can't open input file";

my $count_all = 0;
my $count_prok = 0;
my $count_line = 0;

while (my $line = readline($FILE)) {
	chomp $line;
	$count_line++;
	if ($count_line == 1) {
		next;
	}
	my @split = split /\t/, $line;
	$count_all++;
	if ($split[7] =~ /23S_rRNA/) {
		$count_prok++;
	}
	if ($split[7] =~ /16S_rRNA/) {
		$count_prok++;
	}
}

print "\nTotal number of rRNA: $count_all\n";
print "Total number of prokaryptic 16S and 23S rRNA: $count_prok\n\n";