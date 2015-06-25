#!/usr/bin/perl  

#perl EFRshift.pl -G -I -T 

#use strict;
use File::Basename;
use IO::File;
use Getopt::Std;

my %opts=();
getopts('G:I:T:m:x:z:e:p:v:b:h',\%opts);

$gene=$opts{G};
$intron=$opts{I};
$integenic=$opts{T};

&usage if ($opts{h});
&usage if (!defined $gene);
&usage if (!defined $intron);
&usage if (!defined $integenic);

print "... ... ...\nStart Gene file Processing\n... ... ...\n";
print "... ... ...\n";
my $name=$gene;
if($name=~/\.gtf/){$name=~s/\.gtf$//i;}
if($name=~/gene.*$/){$name=~s/gene.*$//i;};
$name=~s/\_$//i;

my($outName, $directories, $suffix) = fileparse($name);

print "Gene file name: $gene\n";
print "Gene output file name: $outName\n";

my $exon=$outName."_exon.gtf";
my $grFilter=$outName."_GeneRegionFilter.bed";
my $exonMerge=$outName."_GeneRegionFilter_merge.bed";
my $intronMerge=$outName."_introns_merge.bed";
my $intronInter=$outName."_introInter_sort.bed";
my $geneFreeRegion=$outName."_GeneFreeRegion.bed";
my $gfrm=$outName."_GFRM.bed";

#Get exon only file, sort and combine overlap region
print "... ... ...\nGenerate exon and exon_merge file\n... ... ...\n";
system("grep -P '\texon\t' $gene | sort -k1,1 -k4,4n -k5,5n > $exon");
print "Add exon region to geneRegionFilter: $exon\n\n";
system("grep -P '\texon\t' $gene | cut -f1,4,5 > $grFilter");

#---optional gene-exon remove region-------                                                                                                                                                                                                     #
my $mRNA=$opts{m};
if (defined($mRNA)){
    print "Add mRNA region to GeneRegionFilter: $mRNA\n";
    system("sed '/start/I d' $mRNA >> $grFilter");
}
my $xenoMrna=$opts{x};
if (defined($xenoMrna)){
    print "Add xenoMran region to GeneRegionFilter: $xenoMrna\n";
    system("sed '/start/I d' $xenoMrna >> $grFilter");
}
my $xenoRefGene=$opts{z};
if (defined($xenoRefGene)){
    print "Add xenoRefGene region to GeneRegionFilter: $xenoRefGene\n";
    system("sed '/start/I d' $xenoRefGene >> $grFilter");
}
my $ensGene=$opts{e};
if (defined($ensGene)){
    print "Add eneGene region to GeneRegionFilter: $ensGene\n";
    system("sed '/start/I d' $ensGene >> $grFilter");
}
my $pseudoGene=$opts{p};
if (defined($pseudoGene)){
    print "Add pseudoGene region to GeneRegionFilter: $pseudoGene\n";
    system("sed '/start/I d' $pseudoGene >> $grFilter");
}
my $vegaGene=$opts{v};
if (defined($vegaGene)){
    print "Add vegaGene region to GeneRegionFilter: $vegaGene\n";
    system("sed '/start/I d' $vegaGene >> $grFilter");
}
my $optBed=$opts{b};
if (defined($optBed)){
    print "Add optional bed region to GeneRegionFilter: $optBed\n";
    system("sed '/start/I d' $optBed >> $grFilter");
}
#
#---------------------   
system("cat $grFilter | sort -k1,1 -k2,2n -k3,3n | bedtools merge -i stdin -n > $exonMerge");

#Sort and merge overlap intron region
print "... ... ...\n\nGenerate intron merge file\n... ... ...\n";
system("sort -k1,1 -k2,2n -k3,3n $intron | bedtools merge -i stdin -n | gawk 'BEGIN {FS=OFS=\"\t\"} {\$4=\"intron\"; print}' > $intronMerge");

#combine meged intron and integenic file
print "... ... ...\nCombine intron and intergenic regions\n... ... ...\n";
system("cat $intronMerge $integenic | sort -k1,1 -k2,2n -k3,3n > $intronInter");

#remove the region overlapping exon regions from hg19_UCSC_introInter_sort.bed
print "... ... ...\nFind gene free region\n... ... ...\n";
system("bedtools subtract -a $intronInter -b $exonMerge > $geneFreeRegion");
system("awk 'BEGIN {FS=OFS=\"\t\"} {if(\$4==\"intron\"){\$2=\$2+100; \$3=\$3-100; \$5=\$3-\$2}else{\$2=\$2+1000; \$3=\$3-1000; \$5=\$3-\$2} {print \$1, \$2, \$3, \$4, \$5}}' $geneFreeRegion | sort -k1,1 -k2,2n -k3,3n > $gfrm");

###
###
#Comapre gene.gtf with $gfrm and shift the exon to the closest gene free region
print "... ... ...\n";
print "*******************************************************\n";
print "Comapre $exon with $gfrm and shift the exon to gene free regions\n";
print "*******************************************************\n";
print "... ... ...\n";
system("perl /home/UTHSCSA/zou/HiSeq/ChenLab/GeneIntronShift/exonFreeRegionShift.pl -EX $exon -FR $gfrm");

print "Gene exon shifting finished\n... ... ...\n"; 
print "... ... ...\n"; 
print "... ... ...\n"; 
print "... ... ...\n"; 
#**************************************************************#

sub usage {
    die qq(
Program: GEFRshift.pl
Usage: GEFRshift.pl <-G gene-GTF.gtf > <-I intronRegion.tsv> <-T integenicRegion.tsv>
       optional: -m mRNA.bed -x xenoMrna.bed -z xenoRefGene.bed -e ensGene.bed -p pseudoGene.bed -v vegaGene.bed -b customizedBed.bed for filtering 

    \n);
}

sub pr_usage{
    print STDERR "ERROR: $_[0]\n";
    exit;
}
