#!/usr/bin/perl -w

use strict;
use Getopt::Long;

#Command line parameters:
my $distance = "";
my $config = "";
my $diet1 = "";
my $diet2 = "";

#Setup the command line options using Getopt:Long
my $commandline = GetOptions("distance:s", \$distance,
			"diet1:s", \$diet1,
			"diet2:s", \$diet2,
			"config:s", \$config);

if (!$commandline || $distance eq "" || $config eq "" || $diet1 eq "" || $diet2 eq "" ) {
	print STDERR "\nUsage: $0 -distance -config -diet1 -diet2  \n\n";
	exit;
}

open (my $CONFIG, "$config") or die "Can't open config file!";

my $config_line_number = 0;
my %group_hash;
my $col = 0;

while (my $line = readline($CONFIG)) {
	chomp $line;
	$config_line_number ++;
	if ($config_line_number == 4) {
		my @line_split = split /\s/, $line;
		foreach my $sample (@line_split) {
			$col++;
			$group_hash{$col} = $sample;
		}
	}
}

open (my $DIST, "$distance") or die "Can't open $distance file!" ;
open (my $OUTPUT, ">distance_$diet1"."_$diet2".".txt") or die "Can't open output file!" ;

my $line_number = -1;
my %line_hash;
my @dist;

while (my $line = readline($DIST)) {
	chomp $line;
	$line_number++;
	if ($line_number == 0) {
		my @sample_split = split /\t/, $line;
		my $count = -1;
		foreach my $sample (@sample_split) {
			$count++;
			$line_hash{$count} = $sample;
		}
		next;
	}
	my @line_split = split /\t/, $line;
	my $sample_count = 0;
	my $columns = -1;
	my $id = "";
	my $cur_diet = $group_hash{$line_number};
	if (!($cur_diet eq $diet1)) {
		next;
	}
	
	foreach my $sample (@line_split) {
		$columns ++;
		$sample_count++;
		if ($columns == 0) { #skip the column with sample id information
			$id = "$sample";
			next;
		}
		
		if ($columns == $sample_count) { #skip the column matching the current sample
			next;
		}
			
		my $col_diet = $group_hash{$columns};
		my $col_sample = $line_hash{$columns};
			
		if ($cur_diet eq $diet1) {
			if ($col_diet eq $diet2) {
				if (!($id eq $col_sample)) {
					print $OUTPUT "$id\t$cur_diet\t$col_sample\t$col_diet\t$sample\n";
					push @dist, $sample;
				}
			}	
		}
	}
	
}

my $count = 0;
my $total = 0;
foreach my $ele (@dist) {
	$count++;
	$total = $total + $ele;
}

my $average = ($total / $count );

print "For $diet1 and $diet2 the average distance between samples is $average\n\n";
		
