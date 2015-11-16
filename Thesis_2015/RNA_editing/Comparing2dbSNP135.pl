#!/usr/bin/perl
use strict;
use warnings;
use DBI;
use DBD::mysql;
use LWP::Simple;
use Getopt::Long;


print "\n";
print "#######################################################\n";
print "#Author is not me, the script was adapted from Helen's#\n";
print "#######################################################\n";
print "\n";
print "\n";
my($input_file) = @ARGV;
my $output_file = "$input_file.dbSNP135.Matches.txt";
print "SNP data will be taken from: $input_file\n";
print "Output data will be in: $output_file\n";
print "The script will take the first 3 colums as chromosome_name, position and Allele1/Allele2.\n";

#CALL: /users/nhrg/AnalysisScripts/VarScan2PolyPhenWithCCDS_hg19.pl --geno "homozygous"  --pat "patientname" --lane "lane_no." --file "snps_on_target.txt"
#my $patient="Patient";
#my $lane="Lane";
#my $vs_file="vs_file"; #probableSNPs file

#my $het_hom="all"; #SNP type
#my $Results=GetOptions("geno=s"=>\$het_hom, "pat=s"=>\$patient, "lane=s"=>\$lane, "file=s"=>\$vs_file);

my $CurrentDir=`pwd`;
chomp $CurrentDir;

#my $dbSNPmatchfile = MySQLforVariants($CurrentDir, $vs_file, $patient, $lane, $het_hom);
my $dbSNPmatchfile = MySQLforVariants($input_file, $output_file);
                        
#my $dbSNPmatchfile = "AlldbSNPmatches_Lane4_CuB_new_filter_all.txt";                                    

#FindRareVariants($CurrentDir, $dbSNPmatchfile, $vs_file, $patient, $lane, $het_hom);

#my $plhh=$patient."_".$lane."_".$het_hom;
#`rm ReformattedVariantFile_$plhh.txt`;

exit;

#################### SUBROUTINES ##########################

sub MySQLforVariants{
	#my ($CurrentDir, $probableSNPs, $patient, $lane, $het_hom)=@_;
	#my $IDplh=$patient."_".$lane."_".$het_hom;
	my ($input_file, $output_file)=@_;
	my @OutputFile;
	
	open ALLSNPS, $input_file or die "Cannot open $input_file"; 
	#my $header=<ALLSNPS>;
	
	A_loop: while(<ALLSNPS>){
		chomp $_;
		my @SplitLine=split(/\t/, $_);
		#if($het_hom=~/heterozygous/ and $SplitLine[8]>85){next A_loop;}
		#if($het_hom=~/homozygous/ and $SplitLine[8]<=85){next A_loop;}
		my $Chr=$SplitLine[0];
		my $Pos=$SplitLine[1];
		my $Mut=$SplitLine[2];
		#my $ComplimentMutant=$SplitLine[5];
		#$ComplimentMutant=~tr/ACGT/TGCA/;
		push @OutputFile, $Chr.",".$Pos.",".$Mut."\n";#only forward strand mut listed in file, no (rev strand) complement listed
	}
	close ALLSNPS;
	
	my $ReformattedOutput = "Comparing2dbSNP135.temp.txt";
	open REFORMAT, ">$ReformattedOutput";
	print REFORMAT @OutputFile;
	close REFORMAT;
	
	print "intermidiate file is generated with reformated SNPs info from $input_file.\n";
	
	
	##MySQL comparison.
	print"Connecting the database....\n";
	my $ds="DBI:mysql:Data:headnode";
	my $user="nhrg";
	my $dbh=DBI->connect($ds,$user)||die "Cant connect to MySQL";
	print"Connection is made.\n";
	
	#####Generate a table containing all the exons which should be covered by our baits.
	print"Querying the databse......\n";
	my $AllMatchingSNPsTable = $dbh->prepare("CREATE TABLE AllMatchingSNPs (col1 CHAR(5), col2 INT(15), col3 CHAR(25), col4 CHAR(15), col5 CHAR(10), col6 CHAR(54), col7 CHAR(10))");
	$AllMatchingSNPsTable -> execute();
	my $LocatedSNPsTable = $dbh -> prepare("CREATE TABLE LocatedSNPs (col1 CHAR(5), col2 INT(15), col3 CHAR(25))");
	$LocatedSNPsTable -> execute();
	my $LoadLocatedSNPsTable = $dbh -> prepare("LOAD DATA LOCAL INFILE '$CurrentDir/$ReformattedOutput' INTO TABLE LocatedSNPs FIELDS TERMINATED BY ',' LINES TERMINATED BY '\n'");
	$LoadLocatedSNPsTable -> execute();
	my $Index1 = $dbh -> prepare("CREATE INDEX col_index1 ON LocatedSNPs (col1)");
	my $Index2 = $dbh -> prepare("CREATE INDEX col_index2 ON LocatedSNPs (col2)");
	my $Index3 = $dbh -> prepare("CREATE INDEX col_index3 ON LocatedSNPs (col3)");
	$Index1 -> execute();
	$Index2 -> execute();
	$Index3 -> execute();
	
	print"Comparing........";
	my $LoadMatchingSNPsTable1 = $dbh -> prepare(
	"INSERT INTO AllMatchingSNPs
	SELECT LocatedSNPs.col1,LocatedSNPs.col2,LocatedSNPs.col3,snp135.name,snp135.observed,snp135.alleleFreqs,snp135.alleles
	FROM snp135, LocatedSNPs WHERE snp135.chrom=LocatedSNPs.col1 AND snp135.chromStart=LocatedSNPs.col2"); #output all chr pos matches 
	$LoadMatchingSNPsTable1 -> execute();
	print"Done!\n";
	
	print"Generating result.......";
	my $AllMatchingSNPsStatement = "SELECT col1,col2,col3,col4,col5,col6,col7 FROM AllMatchingSNPs";
	my $AllData = $dbh -> prepare($AllMatchingSNPsStatement);
	$AllData -> execute();
	my $AllDataContent = $AllData -> fetchall_arrayref();
	print"Done!\n";
	
	print"Cleaning tables and intermidiate files......";
	#####Removing all the tables form MySQL.
	my $RemoveAllTable = $dbh -> prepare("DROP TABLE AllMatchingSNPs");
	$RemoveAllTable -> execute();
	my $RemoveLocatedTable = $dbh -> prepare("DROP TABLE LocatedSNPs");
	$RemoveLocatedTable -> execute();
	`rm $ReformattedOutput`;
	print"Done!\n";
	
	
	my @Temp;
	my %RemoveDuplicates;
	my @RemoveDuplicatesArray;
	my $DuplicateCounter = 0;
	foreach my $Line(@$AllDataContent){
		push @Temp, "@$Line\n";
	}
	foreach my $TempLine(@Temp){
		if(!exists $RemoveDuplicates{$TempLine}){ ##There will be some duplicates as in matching to dbSNP, some of the reverse compliments are the same-such as C/G (Leads to G/C,C/G,GC).
		$RemoveDuplicates{$TempLine} = $DuplicateCounter;
		push @RemoveDuplicatesArray, $TempLine;
		$DuplicateCounter++;
		}
	}
	
	print "###########################################\n";
	print "########### output results ################\n";
	print "###########################################\n";
	print "$output_file\n";
	my $AllMatchesOutput = $output_file;
	open ALL, ">$AllMatchesOutput";
	print ALL @RemoveDuplicatesArray;
	#print ALL @Temp;
	close ALL;
	return($AllMatchesOutput);
}

sub FindRareVariants{
	###PUT INTO CODE###
	###Are allele freqs defined? No=next line.
	####For each dbSNP match, does Var-allele or Complement of Var-allele match any of the observed alleles with freqs? No=next line. N.B sometimes observed but no freqs present
	###If more than one allele freq. due to complement, choose rarest allele freq, if rare (MAF<0.01) =next line.
	###Common SNP match found = add to %dbSNP.
	###################
	
	my ($CurrentDir, $dbSNPmatchfile, $probableSNPs, $patient, $lane, $het_hom)=@_;
	
	my %dbSNP;
	my @NotIndbSNP;
	my $common_dbSNPcount=0;
	open DBSNP, $dbSNPmatchfile or die "cannot open $dbSNPmatchfile";
	snploop: while (<DBSNP>) {
		my $ln=$_;
		chomp $ln;
		my @line = split(/\s+/, $ln);
		my $v = $line[2];
		my $cv = $line[2];
		$cv=~tr/ACGT/TGCA/;
		if(!defined $line[5]){next snploop;}#don't inlcude SNPs with undefined MAFs in the dbSNPs to filter list
		my @MAFs = split(/,/, $line[5]);
		my @Alleles = split(/,/,$line[6]);
		my $numAlle=@Alleles;
		my %Alls=(); 
		for(my $a=0; $a<$numAlle; $a++){
			$Alls{$Alleles[$a]}=$MAFs[$a];
		}
		if(!exists $Alls{$v} and !exists $Alls{$cv}){next snploop;}#variant allele not observed
		my $lowest=1;
		if(exists $Alls{$v}){$lowest=$Alls{$v};}
		if(exists $Alls{$cv} and $Alls{$cv}<$lowest){$lowest=$Alls{$cv};}
		if($lowest<=0.01){next snploop;}#don't include rare dbSNP alleles 
		if(!exists $dbSNP{$line[0]."_".$line[1]."_".$v}) {
			$dbSNP{$line[0]."_".$line[1]."_".$v}="indbSNP";
			$common_dbSNPcount++;
		}
	}
	print "$patient $lane $het_hom:\nCommon dbSNPs $common_dbSNPcount\n";
	close(DBSNP);
	
	open ALLSNPS, $probableSNPs or die "cannot open $probableSNPs";
	my $head=<ALLSNPS>;
	var_loop: while (<ALLSNPS>) {
		chomp $_;
		my @Line = split(/\t/, $_);
		if($het_hom=~/heterozygous/ and $Line[8]>85){next var_loop;}
		if($het_hom=~/homozygous/ and $Line[8]<85){next var_loop;}
		my $newstart = $Line[3]; #0-based in dbSNP matches and SnpsOnTarget!!!!
		if (!exists $dbSNP{$Line[1]."_".$newstart."_".$Line[5]}) {
			push @NotIndbSNP, $Line[1]."\t".$newstart."\t".$Line[4]."\/".$Line[5]."\t".$Line[8]."\n";
		}
	}
	my $NonMatchOutput="VariantsNotMatchingdbSNP_MAF-0.01_".$patient."_".$lane."_".$het_hom.".txt";
	open NONMATCH, ">$NonMatchOutput";
	print NONMATCH @NotIndbSNP;
	close NONMATCH;

}


        

