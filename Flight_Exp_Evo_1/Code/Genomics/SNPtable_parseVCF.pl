###############################
# SNPtable_parseVCF.pl
###############################
use strict;
use warnings;
my $NumberOfStrains = shift(@ARGV);
chomp $NumberOfStrains;
my $AcceptableMissing = shift(@ARGV);
chomp $AcceptableMissing;

#stdin is standard input and comes from the command line- piped from zcat
#zcat doesn't work in mac osx

#while (my $line = <STDIN>){

#read file in instead
my $filein = shift(@ARGV);
chomp $filein;
open(FH, "<$filein") or die "can't open file\n";

while (my $line = <FH>){
chomp $line;
if($line =~ /^#CHROM\t/){
	$line =~ s/#CHROM/CHROM/;
	my @blowline = split("\t",$line);
	print "Nmiss\t$blowline[0]\t$blowline[1]\t$blowline[3]\t$blowline[4]";
	for(my $i = 0; $i<$NumberOfStrains; $i++){
		my $temp = $blowline[9+$i];
		print "\tfreq_".$temp."\tN_".$temp;
		}
	print "\n";
	}

# So this is what a line of a VCF file looks like
# this is what I am attempting to parse below
#  X	2594943	.	T	G	2689.23	PASS	AC=2;AF=0.056;AN=36;BaseQRankSum=-3.082;DP=1380;Dels=0.00;FS=1.487;HaplotypeScore=1.3573;InbreedingCoeff=1.0000;MLEAC=2;MLEAF=0.056;MQ=58.74;MQ0=0;MQRankSum=-0.694;QD=33.62;ReadPosRankSum=-0.738;EFF=DOWNSTREAM(MODIFIER||||357|Csat|protein_coding|CODING|FBtr0070501||1),DOWNSTREAM(MODIFIER||||368|Csat|protein_coding|CODING|FBtr0273276||1),UPSTREAM(MODIFIER||||352|eIF2B-beta|protein_coding|CODING|FBtr0070479||1),UPSTREAM(MODIFIER||||504|CG2658|protein_coding|CODING|FBtr0070503||1),UPSTREAM(MODIFIER||||819|CG2658|protein_coding|CODING|FBtr0070502||1),UTR_3_PRIME(MODIFIER||||446|CG2662|protein_coding|CODING|FBtr0070478|3|1)	GT:AD:DP:GQ:PL	0/0:66,1:63:99:0,184,2335	0/0:76,0:76:99:0,214,2600	1/1:0,80:79:99:2741,214,0	0/0:71,0:71:99:0,202,2555	0/0:88,0:88:99:0,259,3369	0/0:68,0:68:99:0,205,2655	0/0:59,0:59:99:0,178,2305	0/0:71,0:71:99:0,208,2686	0/0:64,0:64:99:0,184,2316	0/0:62,0:62:99:0,187,2331	0/0:100,0:100:99:0,292,3651	0/0:64,0:62:99:0,187,2441	0/0:70,0:70:99:0,211,2662	0/0:107,0:105:99:0,316,4088	0/0:79,0:79:99:0,229,3093	0/0:76,0:76:99:0,226,3043	0/0:88,0:88:99:0,262,3370	0/0:89,0:88:99:0,259,3326

unless($line =~ /^##/ || $line =~ /^CHROM\t/){

	my $temp = "";
	my $missing = 0;
	my @ss = split("\t",$line);
	my @biallele = split(",",$ss[4]);
	#  ancestral versus derived both length 1
	if(scalar @biallele == 1){
		$temp = "$ss[0]\t$ss[1]\t$ss[3]\t$ss[4]";
		for(my $i=9; $i<(9+$NumberOfStrains); $i++){
			if($ss[$i] =~ /^[0-9\.]\/[0-9\.]$/){
				$temp = $temp."\tNA\t0";
				$missing++;
				}else{
				my @tt = split(":",$ss[$i]);
				my @qq = split(",",$tt[1]);
				my $A = $qq[0];
				my $B = $qq[1];
				my $TotalCount = $A + $B;
				if($TotalCount==0){
					$temp = $temp."\tNA\t0";
					$missing++;
					}else{
					my $RefFreq = $A/$TotalCount;
					$temp = $temp."\t".sprintf("%.2f",$RefFreq)."\t$TotalCount";
					}
				}
			} # over strains
		if($missing <= $AcceptableMissing){print "$missing\t$temp\n";}
		} # biallelic true SNP
	} # line of vcf file describing a SNP
}  # lines in vcf file
exit();


#To run this script 
#perl SNPtable_parseVCF.pl Nsamples Nmissing  “MyVCF.vcf” > MySNPtable.txt

#for example:
#perl SNPtable_parseVCF.pl 34 0  “MyVCF.vcf” > Dispersal_SNPtable.txt