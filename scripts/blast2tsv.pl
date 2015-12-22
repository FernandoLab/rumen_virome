#!/usr/bin/perl -w

use strict;
use Getopt::Long;

#Command line parameters:
my $blast_file = "";
my $bowtie_dir = "";
my $gene_ko = "";

#Setup the command line options using Getopt:Long
my $commandline = GetOptions("blast_file:s", \$blast_file,
			"gene_ko:s", \$gene_ko,
			"bowtie_dir:s", \$bowtie_dir);

if (!$commandline || $blast_file eq "" || $bowtie_dir eq "" || $gene_ko eq "" ) {
	print STDERR "\nUsage: $0 {-blast_file (tab delim output) -bowtie_dir (has bowtie output in .psl format for each sample) -gene_ko (kegg/genes/ko_genes.list)  } \n\n";
	exit;
}

#inputs and outputs
open (my $TOPHIT, ">blast_top_hit.txt") or die "Can't open output TOP HIT file";
open (my $BLASTFILE, "$blast_file") or die "Can't open input BLAST file";
opendir (my $BOWTIEDIR,"$bowtie_dir") or die "Can't open input BOWTIE directory";
open (my $KO, "$gene_ko") or die "Can't open input KO MAPPING file!";
open (my $ORFASSIGN, ">orf_ko_assignment.txt") or die "Can't open output ORF to KO text file";

#set variables
my $top_hit_count = 0;
my $current = "current";
my %blast_hash;
my %orf_hash;
my %ko_hash;
my %normalize_hash;
my $total_length = 0;
my $orf_assign_print = 0;

while (my $line = readline($BLASTFILE)) {
	chomp $line;
	my @blast_split = split /\t/, $line;
	#will likely have >1 hit per entry so use first entry only
	if ($current eq $blast_split[0]) {
		next;
	}
	my @orf_split = split /\s/, $blast_split[0];
	print $TOPHIT "$orf_split[0]\t$blast_split[1]\n";
	$blast_hash{$orf_split[0]} = $blast_split[1]; #hash key is orf ID and value is gene ID
	$top_hit_count++; #keep track of the number of ORFs with a hit
	$current = $blast_split[0];
}

print "\nNumber of top hits from BLAST file: $top_hit_count\n\n";

#create a hash with the gene ID (nr) as the key and K number as value
while (my $line = readline($KO)) {
	chomp $line;
	my @split_tab = split /\t/, $line; #k number is 0 and genes are 1
	my @split_colon = split /:/, $split_tab[0]; #want split_colon[1] for k number
	$ko_hash{$split_tab[1]} = $split_colon[1]; #gene ID is key and K number is value
}

#Should have a directory with bowtie output for each sample to be analyzed
while (my $file = readdir($BOWTIEDIR)) {
	if ($file =~ /.psl/) { 
    	my @sample_id = split /.psl/, $file; #split in front of the period to get sample ID
    	
    	#create output to store ORF abundances for current sample
		open (my $TSV, ">$sample_id[0]".".orf_abundance.txt") or die "Can't open OUTPUT file";
		#create output to store raw gene abundances for current sample
		open (my $TSV2, ">$sample_id[0]".".kegg_gene_raw_abundances.txt") or die "Can't open OUTPUT file";
		#create output to store raw k number abundances for current sample
		open (my $TSV3, ">$sample_id[0]".".ko_raw_abundances.txt") or die "Can't open OUTPUT file";
		#create output to store ORFs with no k number hits for current sample
		open (my $KONOTHIT, ">$sample_id[0]".".no_ko_hit.txt") or die "Can't open input file";
		
		print $TSV "ID\t$sample_id[0]\n"; #print header with "ID\tSAMPLEID
		print $TSV2 "ID\t$sample_id[0]\n"; #print header with "ID\tSAMPLEID
		print $TSV3 "ID\t$sample_id[0]\n"; #print header with "ID\tSAMPLEID
    	
    	#set variables
    	my %orf_abund;
    	my %gene_hash;
    	my %ko_sample_hash;
    	my $alignment_count = 0;
    	my $orf_hit = 0;
    	my $orf_count = 0;
    	my $total_sum = 0;
    	my $total_sum_ko = 0;
    	my $nr_orf = 0;
    	my $nr_ko = 0;
    	my $gene_no_hit = 0;
    	my $sum_no_hit = 0;
    	my $gene_ko = 0;
    	my $sample_length = 0;
    	my $interest_sum = 0;
    	my $sample_interest_count = 0;
    	my $read_count = 0;
    	
    	open (my $MAPFILE, "$bowtie_dir/$file") or die "Can't open input bowtie file!";
		
		while (my $line = readline($MAPFILE)) {
			chomp $line;
			my $check_line = substr ($line, 0, 1);
			if ($check_line eq "@") { #skip lines with @ in the bowtie output, not needed
				next;
			}
			my @bowtie_split = split /\t/, $line;
			my $seq_length = length($bowtie_split[9]);
			$sample_length = $sample_length + $seq_length; 
			$total_length = $total_length + $seq_length;
			$read_count++;
			if ($bowtie_split[2] eq "*") { # * is no hit
				next;
			} else {
				$alignment_count++;
				++$orf_abund{$bowtie_split[2]}; #use hash to count ORF abund
			}
		}
		while ( my ( $orf, $value ) = each %orf_abund ) {
			$orf_count++;
			print $TSV "$orf\t$orf_abund{$orf}\n"; #print the abundance of each ORF
			#blast_hash key is orf ID and value is gene
			if (exists $blast_hash{$orf}) { #should always, but used as a check
				$orf_hit++;
				my $gene_key = $blast_hash{$orf}; #value is ORF
				push (@{$gene_hash{$gene_key}}, $orf_abund{$orf}); #creates nr hash with gene as key and orf abundance as value
				if (exists $ko_hash{$gene_key}) { #gene ID is key and K number is value
					#$orf_assign_print++;
					$orf_hash{$orf} = "$gene_key\t$ko_hash{$gene_key}"; #ORF ID as key and value is gene \t k number
					#print $ORFASSIGN "$orf\t$gene_key\t$ko_hash{$gene_key}\n";
				}	
			}
		}
	
		foreach my $key (sort keys %gene_hash) { #key is gene ID and value is ORF adundace
			my $sum = 0;
			$nr_orf++;
			foreach my $number (@{$gene_hash{$key}}) {
				$sum = $sum + $number;
			}	
			print $TSV2 "$key\t$sum\n";
			$total_sum = $total_sum + $sum;
			if (exists $ko_hash{$key}) {
				$gene_ko++;
				my $ko_key = $ko_hash{$key};
				push (@{$ko_sample_hash{$ko_key}}, $sum);
			} else {
				$gene_no_hit++;
				$sum_no_hit = $sum_no_hit + $sum;
				print $KONOTHIT "$key\t$sum\n";
			}
		}
		
		foreach my $key (sort keys %ko_sample_hash) {
			my $sum = 0;
			$nr_ko++;
			foreach my $number (@{$ko_sample_hash{$key}}) {
				$sum = $sum + $number;
			}	
			print $TSV3 "$key\t$sum\n";
			$total_sum_ko = $total_sum_ko + $sum;
			
		}
		$normalize_hash{$sample_id[0]} = $read_count;
		
		
		print "Alignment count for $sample_id[0]: $alignment_count\n";
		print "Total Read count $sample_id[0]: $read_count\n";
		print "Redundant gene count with aligned reads for $sample_id[0]: $orf_hit\n";
		print "Non-Redundant gene count with aligned reads for $sample_id[0]: $nr_orf\n";
		print "Total abundance in gene hits from aligned reads for $sample_id[0]: $total_sum\n";
		print "Non-Redundant ko count with aligned reads for $sample_id[0]: $nr_ko\n";
		print "Number of nr genes with ko hits: $gene_ko\n";
		print "Total abundance in ko hits from aligned reads for $sample_id[0]: $total_sum_ko\n";
		print "Number of nr genes without ko hits: $gene_no_hit\n";
		print "Total abundance of reads from genes without ko hits: $sum_no_hit\n\n";
	
	}
	
}


foreach my $key (sort keys %orf_hash) {
	$orf_assign_print++;
	print $ORFASSIGN "$key\t$orf_hash{$key}\n";
}

my $avg_length_sample = $total_length / 15;

my @files = glob ("*.ko_raw_abundances.txt");
foreach my $file (@files) {
	my @sample_id_raw = split /.ko/, $file;
	open (my $RAW, "$file") or die "Can't open input bowtie file!";
	open (my $TSV4, ">$sample_id_raw[0]".".ko_corrected_abundances.txt") or die "Can't open OUTPUT file";
	print $TSV4 "ID\t$sample_id_raw[0]\n";
	my $count = 0;
	while (my $line = readline($RAW)) {
		chomp $line;
		$count ++;
		if ($count == 1) {
			next;
		}
		my @raw_split = split /\t/, $line;
		my $norm1 = ($raw_split[1] / $normalize_hash{$sample_id_raw[0]});
		print $TSV4 "$raw_split[0]\t$norm1\n";
	}
	
}		

print "\nTotal nucleotides for all samples combined: $total_length\n";
print "Average nucleotides per samples: $avg_length_sample\n";
print "Number of ORFs with assignments printed to orf_assignment.txt: $orf_assign_print\n";

close ($TOPHIT) or die "Can't close top hit output!";
close ($BLASTFILE) or die "Can't close input blast file!";
close ($BOWTIEDIR) or die "Can't close bowtie dir!";
close ($KO) or die "Can't close gene_ko mapping file!";
close ($ORFASSIGN) or die "Can't close ORF Assignment output file!";