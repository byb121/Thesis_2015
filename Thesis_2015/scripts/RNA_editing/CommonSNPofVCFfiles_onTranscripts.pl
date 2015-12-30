#!/usr/bin/perl
use strict;
use warnings;
use Cwd;

my ($folder) = @ARGV;

$folder =~ s/\/$//;
my @files;

my $dir = cwd();
my $output_file_name = "CommonSNPs.txt";
my $output_file = $dir."/".$output_file_name;
my %hash;

opendir (DIR, $folder) or die "Cann't find the directory";
while (my $file = readdir(DIR)) {
	next if ($file =~ m/^\./);
	if ($file =~ m/\.onTranscripts.vcf$/){ ## This line controls the name pattern of files that are read in
		my $input_file = $folder."/".$file;
		open INPUT, $input_file or die "cannot open file $input_file.\n";
		
		print "reading ".$input_file."\n";
		while (my $line =<INPUT>) {
			chomp $line;
			if( $line =~ m/^chr/ ){
				my @words = split("\t", $line);
				my $key = $words[0]."_".$words[1];
				if (exists $hash{$key}) {
					$hash{$key} += 1;
				} else {
					$hash{$key} = 1;
				}
			} 
		}
		close INPUT;
	}

}
closedir(DIR);

open OUTPUT, ">$output_file" or die "cannot open file $output_file.\n";

for my $key ( keys %hash ) {
	if ($hash{$key} == 16) {
		my @words = split("_", $key);
		print OUTPUT $words[0]."\t".$words[1]."\n";
	}
}

close OUTPUT;
exit;



