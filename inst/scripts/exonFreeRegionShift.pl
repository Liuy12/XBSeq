#!/usr/bin/perl  

#exonFreeRegionShift.pl <-EX exon-GTF file > <-FR gene free region>

use strict;

my $exon="";
my $freeRegion="";

use Getopt::Long;

GetOptions('EX=s' => \$exon,	
	   'FR=s' => \$freeRegion);

&usage if ($exon eq "");
&usage if ($freeRegion eq "");

my $output=$exon;
my $output_shift=$exon;
my $output_noShift=$exon;
$output =~ s/\.gtf/\_shift_All\.gtf/;
$output_shift =~ s/\.gtf/\_shift\.gtf/;
$output_noShift =~ s/\.gtf/\_noShift\.gtf/;


print "... ... ...\nOutput: $output\n... ... ...\n";
print "... ... ...\nStart Processing\n... ... ...\n";
print "... ... ...\n";

#
#read the gene free region and parse region data into the array
#
print "... ... ...\nRead gene free region into array\n";
print "... ... ... ...\n";
print "... ... ... ... ...\n";
my @chrom=();
my @start=();
my @end=();
my @type=();
my @length=();

open(FR, "<$freeRegion") or die "Can't open $freeRegion\n";   
while ( my $lineFR = <FR> ) {
    my @thisRegionInfo=();
    chomp($lineFR);
    if ( $lineFR =~ /^chr/ ) {
	@thisRegionInfo=split(/\t/, $lineFR);
	push @chrom, $thisRegionInfo[0];
	push @start, $thisRegionInfo[1];
	push @end, $thisRegionInfo[2];
	push @type, $thisRegionInfo[3];
	push @length, $thisRegionInfo[4];
    }
    else{}
}
close FR;

my $totalFrNo=scalar @start;
print "Total gene free regions: $totalFrNo\n";

#
#read the gene free region and parse region data into the array
#
print "... ... ...\nCompare exon coordinates and shift to gene free region\n";
print "... ... ...\n";
print "... ... ... ...\n";
print "... ... ... ... ...\n";

open (OUT, ">$output");
open (EX, "<$exon") or die "Can't open $exon\n";   

my $lastExChrom="";
my $lastExStart="";
my $lastExEnd="";
my $lastExShiftChrom="";
my $lastExShiftStart="";
my $lastExShiftEnd="";
my $lastHitNotes="";
my $lastExAtt="";

while ( my $lineEX = <EX> ) {
    chomp($lineEX);

    my @thisExonInfo=();
    my $thisExChrom="";
    my $thisExStart="";
    my $thisExEnd="";
    my $thisExLength=0;

    my $thisExShiftStart="";
    my $thisExShiftEnd="";
    my $thisExShiftChrom="";
    my $thisExAtt="";

    my $shiftDistance=0;
    my $shiftRegion="";
    my $shiftTag=0;
    my $outTag=0;

    my $hitNotes="";
    
    if ( $lineEX =~ /^chr/ ) {
	@thisExonInfo=split(/\t/, $lineEX);
	$thisExChrom=$thisExonInfo[0];
	$thisExStart=$thisExonInfo[3];
	$thisExEnd=$thisExonInfo[4];
	$thisExLength=$thisExEnd-$thisExStart;
	$thisExAtt=$thisExonInfo[8];

	#print "Gene: $thisExChrom\t$thisExStart\t$thisExEnd\n";
	if($thisExChrom eq $lastExChrom && $thisExStart eq $lastExStart && $thisExEnd eq $lastExEnd){
	    if($lastHitNotes=~/no_Shift/){
		$lineEX=~s/\"\;$/$lastHitNotes/;
		print OUT "$lineEX\n";
	    }
	    else{
		$thisExAtt=~s/\"\;$/$lastHitNotes/;
		print OUT "$lastExShiftChrom\t$thisExonInfo[1]\t$thisExonInfo[2]\t$lastExShiftStart\t$lastExShiftEnd\t$thisExonInfo[5]\t$thisExonInfo[6]\t$thisExonInfo[7]\t$thisExAtt\n";
	    }
	    $shiftTag=0;
	    $lastExAtt=$thisExAtt;
	}
	else{
	    for(my $i=0; $i<$totalFrNo; $i++){
		if($chrom[$i] eq $thisExChrom && $start[$i]>$thisExEnd){
		    for (my $j=0; $j<$totalFrNo; $j++){
			my $down=$i+$j;
			my $up=$i-$j-1;
			if($chrom[$down] eq $thisExChrom && $thisExLength<$length[$down] && $down<$totalFrNo){
			    #print "Down: $chrom[$down]\t$start[$down]\t$end[$down]\n\n";
			    $thisExShiftStart=$start[$down];
			    $shiftDistance=$thisExShiftStart-$thisExStart;
			    $thisExShiftEnd=$thisExEnd+$shiftDistance;
			    $start[$down]=$thisExShiftEnd+1;
			    $length[$down]=$end[$down]-$start[$down];
			    $thisExShiftChrom=$chrom[$down];
			    $shiftRegion=$type[$down];
			    $shiftTag=1;
			    last;
			}
			elsif($chrom[$up] eq $thisExChrom && $thisExLength<$length[$up] && $up>=0){
			    #print "Up:   $chrom[$up]\t$start[$up]\t$end[$up]\n\n";
			    $thisExShiftEnd=$end[$up];
			    $shiftDistance=$thisExShiftEnd-$thisExEnd;
			    $thisExShiftStart=$thisExStart+$shiftDistance;
			    $end[$up]=$thisExShiftStart-1;
			    $length[$up]=$end[$up]-$start[$up];
			    $thisExShiftChrom=$chrom[$up];
			    $shiftRegion=$type[$up];
			    $shiftTag=1;
			    last;
			}
			elsif($chrom[$down] ne $thisExChrom && $chrom[$up] ne $thisExChrom){
			    last;
			}
			else{next;}
		    } 
		    if($shiftTag=1){
			$hitNotes="|gene_coord:" . $thisExChrom . "_" . $thisExStart . "_" . $thisExEnd . "|shiftRegion:" . $shiftRegion . "|shiftDistance:" . $shiftDistance . "\"\;";
			$thisExAtt=~s/\"\;$/$hitNotes/;		
			print OUT "$thisExShiftChrom\t$thisExonInfo[1]\t$thisExonInfo[2]\t$thisExShiftStart\t$thisExShiftEnd\t$thisExonInfo[5]\t$thisExonInfo[6]\t$thisExonInfo[7]\t$thisExAtt\n";
			$shiftTag=0;
		    }
		    else{
			$hitNotes="|no_shift\"\;";
			$thisExShiftChrom=$thisExChrom;
			$thisExShiftStart=$thisExStart;
			$thisExShiftEnd=$thisExEnd;
			$lineEX=~s/\"\;$/$hitNotes/;
			print OUT "$lineEX\n";
			$shiftTag=0;	
		    }
		    $lastExShiftChrom=$thisExShiftChrom;
		    $lastExShiftStart=$thisExShiftStart;
		    $lastExShiftEnd=$thisExShiftEnd;
		    $lastHitNotes=$hitNotes;
		    $lastExAtt= $thisExAtt;
		    last;
		}
		else{next;}
	    }
	}
	$lastExChrom=$thisExChrom;
	$lastExStart=$thisExStart;
	$lastExEnd=$thisExEnd;
    }else{}
}
close EX;
close OUT; 

#Divide the shift file to shift and noShift file
system("grep 'no_shift' $output > $output_noShift");
system("grep -v 'no_shift' $output > $output_shift");
 
print "Search finished\n... ... ...\n"; 
print "... ... ...\n"; 
print "... ... ...\n"; 
print "... ... ...\n"; 
#**************************************************************#

sub usage {
    die qq(
Program: exonFreeRegionShift.pl
Usage: exonFreeRegionShift.pl <-EX exon-GTF file > <-FR gene free region>
\n);
}

sub pr_usage{
    print STDERR "ERROR: $_[0]\n";
    exit;
}
