#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
#use DBI;
#use DBD::mysql;
use LWP::Simple;
use Getopt::Long;

################ Intro ##################
# The script will look for those variants that pass the filters.
# If they are exist in at least half of total sampels, 
# the script will look for their position in Gencode V10.
# 
# Genotype of left variants will be outputed into results file.
#########################################


print "\n";
print "#################################################################\n";
print "# Author is Yaobo. The script should be faster than SQL version #\n";
print "#################################################################\n";
print "\n";
print "\n";
print "\n";
print "\n";

my ($folder) = @ARGV;

my $common_heterozygous = Looking_4_common_hereozygous_in_VCF_file_folder($folder);

my $output_file_no_matched = "$common_heterozygous.dbSNP135.noSQL.filtered.txt";
my $output_file_matched = "$common_heterozygous.dbSNP135.noSQL.Matches.txt";
my $final_output = "$common_heterozygous.dbSNP135.genotype.filtered.txt";

my $SNP_file = "/users/a5907529/GenomeData/hg19/snp135.txt";

print "Variants data will be taken from: $common_heterozygous.\n";
print "Filtered result will be in: $output_file_no_matched.\n";
print "Matched result will be in: $output_file_matched.\n";
print "The script will take the first 3 colums as chromosome_name, position and Allele1/Allele2.\n";
print "Comparing to SNPs recorded in $SNP_file. Make sure that the file exists."."\n";

Comparing_2_DB_SNPfile($common_heterozygous, $SNP_file, $output_file_no_matched, $output_file_matched);	

my $filteredresult = DetermineIfRNAediting($output_file_no_matched);

print "#################################################\n";
print "########### output final results ################\n";
print "#################################################\n";

open FINAL, ">$final_output";
print FINAL @$filteredresult;
close FINAL;

exit;

sub Looking_4_common_hereozygous_in_VCF_file_folder {
	
	my ($folder) = @_;
	$folder =~ s/\/$//;
	my @files;
	my $dir = cwd();
	
	my %hash;
	
	my $depth_filter = 10; #locus with seqencing depth lower than the value will be dropped.
	my $alt_quality_filter = 13; # Phred-scaled p-value of the alternative base calling is wrong
	my $gq_quality_filter = 13; # Phred-scaled p-value of the genotype calling is wrong
	
	opendir (DIR, $folder) or die "Cann't find the directory";
	while (my $file = readdir(DIR)) {
		next if ($file =~ m/^\./);
		if ($file =~ m/\.flt.vcf$/){ ## This line controls the name pattern of files that are read in
			my $input_file = $folder."/".$file;
			open INPUT, $input_file or die "cannot open file $input_file.\n";
			
			print "reading ".$input_file."\n";
			while (my $line =<INPUT>) {
				chomp $line;
				if( $line =~ m/^chr/ ){
					my @words = split("\t", $line);
					my $key = $words[0]."_".$words[1];
					
					##filters: include Sequencing depth filter, ALT quaility filter and GenoType call quality filter
					my $depth = $words[5];
					my $alt_quality;
					my @vcf_info = split (";", $words[7]);
					foreach my $tag(@vcf_info) {
						if ($tag =~ m/^DP\=/) {
							$tag =~ s/^DP\=//;
							$alt_quality = $tag;
						}
					}
					my @vcf_sample = split (":", $words[9]);
					my $gq_quality = $vcf_sample[2];
					
					if($depth >= $depth_filter && $alt_quality >= $alt_quality_filter && $gq_quality >= $gq_quality_filter) {
						if (exists $hash{$key}) {
							$hash{$key} += 1;
						} else {
							$hash{$key} = 1;
						}	
					}
				} 
			}
			close INPUT;
		}
	}
	closedir(DIR);
	
	print "################################################################\n";
	print "#                     Asign SNP to genes                       #\n";
	print "#REF/Var could be biased, check the script it self for detail./#\n";
	print "################################################################\n";
	print "\n";
	print "######\n";
	print "Warning: Files will only have the ref and var bases extracted from the last vcf file processed. God blessed all samples have the same ref and var!\n";
	print "######\n";
	print "\n";
	
	my %annoted_hash;
	my $temp_output_file = $dir."/"."temptemptemp_file.txt";
	my $temp_output_file_2 = $dir."/"."temptemptemp_file_with_genes_Names.txt";
	
	open OUTPUT, ">$temp_output_file" or die "cannot open file $temp_output_file.\n";
	for my $key ( keys %hash ) {
		if ($hash{$key} >= 8) { 												#select those variants that exists in half of my samples
			my @words = split("_", $key);
			print OUTPUT $words[0]."\t".$words[1]."\t".$words[1]."\n";
		}
	}
	close OUTPUT;
	
	`/users/a5907529/biosofts/BEDTools-Version-2.12.0/bin/intersectBed -wo -a $temp_output_file -b ~/GenomeData/hg19/gencode.v10.annotation.gtf > $temp_output_file_2`;
	
	open INPUT, $temp_output_file_2 or die "cannot open file $temp_output_file_2.\n";
	while (my $line =<INPUT>) {
		chomp $line;
		if( $line =~ m/^chr/ ){
			my @words = split("\t", $line);
			my $gene_name;
			my $gene_id;
			my $key = $words[0]."_".$words[1];
			if (exists $words[11]) {
				my @attributes = split ("; ", $words[11]);
				#print "word 11 are::";
				#print @attributes;
				#print "\n";
				
				foreach my $attribute (@attributes) {
					if ($attribute =~ m/^gene\_id/){
						$attribute =~ s/gene\_id//;
						$attribute =~ s/\"//g;
						#print "subtitued gene_id is:";
						#print $attribute;
						#print "\n";
						$gene_id = $attribute;
					}
					
					if ($attribute =~ m/^gene\_name/){
						$attribute =~ s/gene\_name//;
						$attribute =~ s/\"//g;
						#print "subtitued gene_name is:";
						#print $attribute;
						#print "\n";
						$gene_name = $attribute;
					}
				}
				$annoted_hash{$key} = $gene_name."_".$gene_id;	
			}
		}
	}
	close INPUT;
	
	#`rm $temp_output_file`;
	#`rm $temp_output_file_2`;
	
	for my $key ( keys %hash ) {
		if ($hash{$key} >= 8 && ! exists $annoted_hash{$key} ) { 												#select those variants that exists in half of my samples
			$annoted_hash{$key} = "NA/Intergenic"."_"."NA/Intergenic";
		}
	}
	
	
	
	##Counting Genotype#
	
	my %genoType_hash;
	my %ref_var_hash; # this hash will only contain the ref and var from the last vcf file processed. God blessed all samples have the same ref and var!"
	my %homo_count;
	my %hetre_count;
	
	my $file_counts = 0;
	#my %export_2_R;
	
	opendir (DIR, $folder) or die "Cann't find the directory";
	while (my $file = readdir(DIR)) {
		next if ($file =~ m/^\./);
		if ($file =~ m/\.flt.vcf$/){ ## This line controls the name pattern of files that are read in
			$file_counts += 1;
			my $input_file = $folder."/".$file;
			open INPUT, $input_file or die "cannot open file $input_file.\n";
			
			print "reading ".$input_file."\n";
			while (my $line =<INPUT>) {
				chomp $line;
				if( $line =~ m/^chr/ ){
					my @words = split("\t", $line);
					my $key = $words[0]."_".$words[1];
					$ref_var_hash{$key} = $words[3]."/".$words[4];
					if (exists $annoted_hash{$key}) {
						my @vcf_info = split (";", $words[7]);
						my @vcf_sample = split (":", $words[9]);
						my $ref2var_coverage;
						foreach my $tag(@vcf_info) { # it requires that the vcf file info column has to have a DP4 tag
							if ($tag =~ m/^DP4/) {
								$tag =~ s/^DP4\=//;
								my @numbers = split (",", $tag);
								my $temp1 = $numbers[0] + $numbers[1]; #extract coverage of ref and revers_ref
								my $temp2 = $numbers[2] + $numbers[3]; #extract coverage of var and revers_var
								$ref2var_coverage = $temp1.":".$temp2; #recorde the ratio
								if (exists $genoType_hash{$key}) {
									$genoType_hash{$key} = $genoType_hash{$key}."\t".$vcf_sample[0]."_".$ref2var_coverage; #genotye (eg 0/1 or 1/1) and base ratio
								} else {
									$genoType_hash{$key} = $vcf_sample[0]."_".$ref2var_coverage;
								}
								
							} else {
								next;
							}
						}
						
						if (!exists $hetre_count{$key}) {# initial the count to allow value 0
							$hetre_count{$key} = 0;
						}
						if (!exists $homo_count{$key}) {
							$homo_count{$key} = 0;
						}
						
						if ($vcf_sample[0] =~ m/0\/1/) {# start to count
							$hetre_count{$key} += 1;
						} else {
							$homo_count{$key} += 1;						
						}
						
					} else {
						next;
					}
				} 
			}
			
			for my $key (keys %annoted_hash) { #find those variants that not recored in this vcf file but in annoted_hash, assign an value 0/0_1:0
				if(!exists $genoType_hash{$key}) {
					$genoType_hash{$key} = '0/0_1:0';	
				} else {
					my @genoTypes = split("\t", $genoType_hash{$key});
					if (scalar(@genoTypes) < $file_counts) {
						$genoType_hash{$key} = $genoType_hash{$key}."\t".'0/0_1:0';
					}
				}
			}
			
			close INPUT;
		}
	}
	closedir(DIR);
	
	
	##### To calculate the average imbalance
	###Choose SNP that are hetreozygous in at least 3 samples in OA and NOF group
	###Remove homozygous at the position
	###Return avg(Nr/(Nr+Nv))
	
	my %OA_group = ( 2 => 1, 3 => 1, 4 => 1, 5 => 1, 6 => 1, 7 => 1,  8 => 1, 9 => 1, 11 => 1, 12 => 1);
	my %NOF_group = ( 1 => 1, 10 => 1, 13 => 1, 14 => 1, 15 => 1, 16 => 1);
	my @output;
	my @intermidiate_genotype_output;
	my @output_hetereo_freq;
	
	my $total_count = 0;
	my $hetero_in_NOF_only = 0;
	my $hetero_in_OA_only = 0;
	
	for my $key ( keys %genoType_hash ) {
		my @words = split("_", $key);
		if ($words[0] =~ m/chrY/){
			next;
		} else {
			#print OUTPUT $words[0]."\t".$words[1]."\t".$genoType_hash{$key}."\t".$hetre_count{$key}."\t".$homo_count{$key}."\n";
			my $OA_hetero_count = 0;
			my $NOF_hetero_count = 0;
			
			
	
			my @sample_genoType = split("\t", $genoType_hash{$key});
			for(my $i=0;$i<=$#sample_genoType;$i++){
				#print $key." sample_genoType is ".$sample_genoType[$i]. "\n";
				my @temp = split("_", $sample_genoType[$i]);
				my @numbers = split(":", $temp[1]); 
				if($temp[0] =~ m/0\/1/) {
					if(exists($OA_group{$i+1})){
						$OA_hetero_count += 1;
					} elsif (exists($NOF_group{$i+1})) {
						$NOF_hetero_count += 1;
					} else {
						print "Fatal error, Grouping failed!"."\n";
						exit;
					}
				}
			}
			
			if($OA_hetero_count == 0 && $NOF_hetero_count != 0) {
				$hetero_in_NOF_only += 1;
			} elsif($OA_hetero_count != 0 && $NOF_hetero_count == 0) {
				$hetero_in_OA_only += 1;
			}
			
			my @words_2 = split("_", $annoted_hash{$key});
	
			if ($OA_hetero_count + $NOF_hetero_count >= 8) {
				
				$total_count += 1;
				my $average_OA;
				my $average_NOF;
				
				my @words_2 = split("_", $annoted_hash{$key});
				my $total_hetereo = $OA_hetero_count+$NOF_hetero_count;
				#next line is used to check results
				#my $output_string = $words[0]."\t".$words[1]."\t".$ref_var_hash{$key}."\t".$words_2[0]."\t".$words_2[1]."\t".$OA_hetero_count."\t".$NOF_hetero_count."\t".$total_hetereo."\t".$genoType_hash{$key}."\n";
				my $output_string = $words[0]."\t".$words[1]."\t".$ref_var_hash{$key}."\t".$words_2[0]."\t".$words_2[1]."\t".$OA_hetero_count."\t".$NOF_hetero_count."\t".$total_hetereo."\n";
			
				push (@output, $output_string);
				push @intermidiate_genotype_output, $words[0]."\t".$words[1]."\t".$genoType_hash{$key}."\n";
			}
		}
	}
	
	print "###########################################\n";
	print "########### output results ################\n";
	print "###########################################\n";
	
		
	my $imbalance_output_file = $dir."/"."Common_Heterozygous_with_info_filtered.txt";
	open OUTPUT, ">$imbalance_output_file" or die "cannot open file $imbalance_output_file.\n";
	print OUTPUT "#Total qulified common heterozygous are: $total_count"."\n";
	print OUTPUT "#Heterozygous in NOF only: $hetero_in_NOF_only"."\n";
	print OUTPUT "#Heterozygous in OA only: $hetero_in_OA_only"."\n";
	print OUTPUT @output; 
	close OUTPUT;
	
	print "Total qulified common heterozygous are: $total_count"."\n";
	print "Heterozygous in NOF only: $hetero_in_NOF_only"."\n";
	print "Heterozygous in OA only: $hetero_in_OA_only"."\n";
	print "Result is in: $imbalance_output_file\n";
	
	my $intermidiate_output_file = $dir."/"."Common_Heterozygous_genotypes_intermidiatefile.txt";
	open OUTPUT, ">$intermidiate_output_file" or die "cannot open file $intermidiate_output_file.\n";
	print OUTPUT @intermidiate_genotype_output; 
	close OUTPUT;
	
	print "\n";
	print "\n";
	print "\n";
	print "Intermidiate file for genotype of each sample is: \n";
	print "it can be deleted if no desired\n";
	print "\n";
	
	return $imbalance_output_file;
}


	
sub Comparing_2_DB_SNPfile {
	my ($input_file, $SNP_file, $output_file_no_matched, $output_file_matched)=@_;
	
	print "\n";
	print "Input data is from: $input_file.\n";
	print "Variants in the file will be compared to SNP recorded in $SNP_file.\n";
	print "\n";
	
	my %query_variants;
	
	open VARIANTS, $input_file or die "Cannot open $input_file"; 
	while (my $line = <VARIANTS>){  
		chomp $line;
		if ($line =~ m/^chr/) {
			my @SplitLine=split(/\t/, $line);
			#if($het_hom=~/heterozygous/ and $SplitLine[8]>85){next A_loop;}
			#if($het_hom=~/homozygous/ and $SplitLine[8]<=85){next A_loop;}
			my $Chr=$SplitLine[0];
			my $Pos=$SplitLine[1];
			my $Mut=$SplitLine[2];
			#my $ComplimentMutant=$SplitLine[5];
			#$ComplimentMutant=~tr/ACGT/TGCA/;
			if (exists $query_variants{$Chr."_".$Pos}) {
				next;
				print "Warning: Ignored a duplicated entry found for ".$Chr."_".$Pos." in query file: $input_file.\n"
			} else {
				$query_variants{$Chr."_".$Pos} = $Mut;	
			}
		}
	}
	close VARIANTS;
	
	my %dbsnp_matches;
	open SNP_DB, $SNP_file or die "Cannot open $SNP_file"; 
	while (my $line = <SNP_DB>){  
		chomp $line;
		my @SplitLine=split(/\t/, $line);
		my $Chr=$SplitLine[1];
		my $Pos=$SplitLine[2]+1;
		my $Mut=$SplitLine[9];
		my $rs=$SplitLine[4];
		my $ref_base=$SplitLine[8];
		my $strand=$SplitLine[6];
		 
		#my $ComplimentMutant=$SplitLine[5];
		#$ComplimentMutant=~tr/ACGT/TGCA/;
		
		if (exists $query_variants{$Chr."_".$Pos}) {
			$dbsnp_matches{$Chr."_".$Pos} = $rs."\t".$Mut."\t".$ref_base."\t".$strand;
		} else {
			next;
		}
	}
	close SNP_DB;
	
	my @matched;
	my @filtered;
	
	push @matched, "chrom"."\t"."position"."\t"."genotype"."\t"."mapped gene"."\t"."ensembl ID"."\t"."heterozygous in OA"."\t"."heterozygous in NOF"."\t"."total heterozygous"."\t";
	push @matched, "rs"."\t"."recorded mutation"."\t"."ucsc ref base"."\t"."strand"."\n";
	push @filtered, "chrom"."\t"."position"."\t"."genotype"."\t"."mapped gene"."\t"."ensembl ID"."\t"."heterozygous in OA"."\t"."heterozygous in NOF"."\t"."total heterozygous"."\n";
	
	open VARIANTS, $input_file or die "Cannot open $input_file"; 
	while (my $line = <VARIANTS>){  
		chomp $line;
		if ($line =~ m/^chr/) {
			my @SplitLine=split(/\t/, $line);
			my $Chr=$SplitLine[0];
			my $Pos=$SplitLine[1];
			#my $Mut=$SplitLine[2];
			#my $ComplimentMutant=$SplitLine[5];
			#$ComplimentMutant=~tr/ACGT/TGCA/;
			if (exists $dbsnp_matches{$Chr."_".$Pos}) {
				push @matched, $line."\t".$dbsnp_matches{$Chr."_".$Pos}."\n";
			} else {
				push @filtered, $line."\n";
			}
		}
	}
	close VARIANTS;
	
	print "###########################################\n";
	print "########### output results ################\n";
	print "###########################################\n";
	print "$output_file_matched\n";
	print "$output_file_no_matched\n";
	
	open MATCHED, ">$output_file_matched";
	print MATCHED @matched;
	close MATCHED;
	
	open FILTERED, ">$output_file_no_matched";
	print FILTERED @filtered;
	close FILTERED;
	
	#return($output_file_matched, $output_file_no_matched);
	
}

sub DetermineIfRNAediting {
	
	my ($input_file)=@_;
		
	print "\n";
	print "Input data is from: $input_file.\n";
	print "Variants in the file will be filtered for RNA-eidting.\n";
	print "Possible RNA-editing events is define in the script as following:.\n";
	print "A-to-I: A/T, A/C, T/A, T/G.\n";
	print "C-to-U: C/T, G/A.\n";
	print "\n";
	
	my @filtered;
	
	open VARIANTS, $input_file or die "Cannot open $input_file"; 
	while (my $line = <VARIANTS>){  
		chomp $line;
		if ($line =~ m/^chr/) {
			my @SplitLine=split(/\t/, $line);
			my $Chr=$SplitLine[0];
			my $Pos=$SplitLine[1];
			my $Mut=$SplitLine[2];
			my $ComplimentMutant=$Mut;
			$ComplimentMutant=~tr/ACGT/TGCA/;
			
			#Determine if the Mut type is an possible RNA-editing event.
			if ($Mut =~ m/^[Aa]\/[Tt]|^[Aa]\/[Cc]|^[Tt]\/[Aa]|^[Tt]\/[Gg]|^[Cc]\/[Tt]|^[Gg]\/[Aa]/ | $ComplimentMutant =~ m/^[Aa]\/[Tt]|^[Aa]\/[Cc]|^[Tt]\/[Aa]|^[Tt]\/[Gg]|^[Cc]\/[Tt]|^[Gg]\/[Aa]/) {
				push @filtered, $line."\n";
			}
		}
	}
	close VARIANTS;
	
	return \@filtered;
}




