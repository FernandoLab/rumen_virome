#!/usr/bin/perl -w

use strict;
use Getopt::Long;

#Command line parameters:
my $phast_blast = "";
my $orf_ko_assign = "";
my $knum = "";

#Setup the command line options using Getopt:Long
my $commandline = GetOptions("phast_blast:s", \$phast_blast,
			"orf_ko_assign:s", \$orf_ko_assign,
			"knum:s", \$knum);
			
if (!$commandline || $phast_blast eq "" || $orf_ko_assign eq "" || $knum eq "" ) {
	print STDERR "\nUsage: $0 {-knum -orf_ko_assign -phast_blast  } \n\n";
	exit;
}

#inputs and outputs
open (my $PHAST, "$phast_blast") or die "Can't open input file";
open (my $ORF_KO, "$orf_ko_assign") or die "Can't open input file";
open (my $KNUM, "$knum") or die "Can't open input file";
open (my $OUT, ">orf_ko_phast_hits.txt") or die "Can't open output file";
open (my $OUT2, ">ko_w_phast_hits_list.txt") or die "Can't open output file";

my %knum_hash;
my %contig_hash;
my %knum_check_hash;
my %knum_out_hash;
my $top_hit_count = 0;
my $current = "current";
my %orf_out;
my %contig_out;

while (my $line = readline($KNUM)) {
	chomp $line;
	++$knum_hash{$line};
}

while (my $line = readline($ORF_KO)) {
	chomp $line;
	my @split = split /\t/, $line;
	if (exists $knum_hash{$split[2]}) {
		my @split2 = split /\_length/, $split[0];
		push(@{$contig_hash{$split2[0]}}, $line);
		#$contig_hash{$split2[0]} = "$line";
		++$knum_check_hash{$split[2]};
	}
}



while (my $line = readline($PHAST)) {
	chomp $line;
	my @blast_split = split /\t/, $line;
	if ($current eq $blast_split[0]) {
		next;
	}
	$top_hit_count++;
	my @split = split /\_length/, $blast_split[0];
	if (exists $contig_hash{$split[0]}) {
		foreach my $nextORF (@{$contig_hash{$split[0]}}) {
			my @nextORF_split = split /\t/, $nextORF;
			print $OUT "$nextORF\t$blast_split[0]\t $blast_split[1]\n";
			++$knum_out_hash{$nextORF_split[2]};
			++$orf_out{$nextORF_split[0]};
			++$contig_out{$split[0]};
		}
	}
	$current = $blast_split[0];
}

my $count_check = 0;
my $count_check2 = 0;
my $count_check3 = 0;
my $count_check4 = 0;

while ( my ( $key1, $value1 ) = each %knum_check_hash ) { #should equal the number of k nums input
	$count_check++;
}

while ( my ( $key1, $value1 ) = each %knum_out_hash ) { 
	$count_check2++;
	print $OUT2 "$key1\n";
}

while ( my ( $key1, $value1 ) = each %orf_out ) { 
	$count_check3++;
}

while ( my ( $key1, $value1 ) = each %contig_out ) { 
	$count_check4++;
}

print "\nNumber of unique differential k numbers found with ORFs: $count_check\n";
print "Total number of PHAST top hits: $top_hit_count\n";
print "Contigs with phast and KO hit: $count_check4 \n";
print "ORFs on a contig with phast hit and ORF has a KO: $count_check3\n";
print "K numbers printed that were on a contig with PHAST hit and KO: $count_check2\n";