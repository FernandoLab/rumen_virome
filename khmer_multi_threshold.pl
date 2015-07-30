#!/usr/bin/perl -w

use strict;
use Getopt::Long;

#Command line parameters:
my $khmer_multi_dir = "";
my $threshold = "";

#Setup the command line options using Getopt:Long
my $commandline = GetOptions("khmer_multi_dir:s", \$khmer_multi_dir,
			"threshold:s", \$threshold);

if (!$commandline || $khmer_multi_dir eq "" || $threshold eq "" ) {
	print STDERR "\nUsage: $0 -khmer_multi_dir -threshold \n\n";
	exit;
}

opendir (my $DIR,"$khmer_multi_dir/") or die "Can't open directory!" ;

my $file_count = 0;
my @thefiles;

while (my $file = readdir ($DIR)) { #read through the file names in directory
	if ($file =~ /.comb/) { #if file extension is '.comb'
      	$file_count++ ; #increment the file count by one
		push @thefiles, $file ; #push each file to this array
	}
}

my $file_number = 0;

foreach my $file (@thefiles) {
    open (my $FILE, "$khmer_multi_dir/$file") or die "Can't open $file from $khmer_multi_dir!" ;
    open (my $OUTPUT, ">$file".".keep.txt") or die "Can't open output file!" ;
    
    $file_number ++;
    my $keep_line = 0;
    my $input_seqs = 0;
    
    while (my $line = readline($FILE)) {
		chomp $line;
		$input_seqs ++;
		my @line_split = split /\s/, $line;
		my $sample_count = 0;
		my $columns = -1;
		foreach my $sample (@line_split) {
			$columns ++;
			if ($columns == 0) { #skip the column with seq id information
				next;
			}
			
			if ($columns == $file_number) { #skip the column matching the current sample
				next;
			}
			
			if ($sample > 1) {
				$sample_count ++;
			}
		}
		
		if ($sample_count >= $threshold) {
			$keep_line ++;
			print $OUTPUT "$line\n";
		}
	
	}
	
	print "$file had $input_seqs sequences and kept $keep_line sequences\n";
	
}

print "\nTotal number of files analyzed: $file_count\n\n";