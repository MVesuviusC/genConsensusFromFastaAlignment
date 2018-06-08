#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;


##############################
# By Matt Cannon
# Date: May 31, 2018
# Last modified: June 8, 2018
# Title: genConsensusFromFastaAlignment.pl
# Purpose: Generate a consensus from fasta alignment
##############################

##############################
# Options
##############################


my $verbose;
my $help;
my $inFile = "";
my $minAf = 0;
my $maxMissing = 50;

# i = integer, s = string
GetOptions ("verbose"           => \$verbose,
            "help"              => \$help,
            "inFile=s"          => \$inFile,
            "minAf=i"           => \$minAf,
            "maxMissing=i"      => \$maxMissing
      )
    or pod2usage(0) && exit;

pod2usage(1) && exit if ($help);


##############################
# Global variables
##############################
my $seqCount = 0;
my %alignmentHash;
my $tooLowCount = 0;

my %ambiguityHash = (
        A     => "A",
        T     => "T",
        C     => "C",
        G     => "G",
        TC    => "Y",
        AG    => "R",
        AT    => "W",
        GC    => "S",
        TG    => "K",
        AC    => "M",
        ATG   => "D",
        AGC   => "V",
        ATC   => "H",
        TGC   => "B",
        ATGC  => "N");


##############################
# Code
##############################

##############################
### check variables input

die "\n\tNo input provided, use --inFile to give input fasta\n\n" if $inFile eq "";
die "\n\t--minAf cannot be greater than 100\n\n" if $minAf > 100;
die "\n\t--maxMissing cannot be greater than 100\n\n" if $maxMissing > 100;


##############################
### Pull in alignment one fasta entry at a time and load the data into a hash

$/  = "\n>";                                                                    # change line input delimiter
open (ALIGN, "<", $inFile) or die "Cannot open input file\n";

while(my $input = <ALIGN>) {
    chomp $input;

    $seqCount++;

    my(undef, @sequences) = split "\n", $input;

    my $sequence = join("", @sequences);
    $sequence = uc($sequence);
    $sequence =~ s/[\t\s]//g;
    $sequence =~ s/[YRWSKMDVHBN]/-/g;                                           # replace ambiguous bases with a dash

    for(my $i = 0; $i < length($sequence); $i++) {                              # load sequence into %alignmentHash
        my $base = substr($sequence, $i, 1);
        $alignmentHash{$i}{$base}++;                                            # hash{base number}{[ATGC-]} = count
    }

    if($verbose) {                                                              # print out progress
        print STDERR $seqCount,"\r";
    }
}

close ALIGN;
$/  = "\n";                                                                     # change input delimiter back

print  ">", $inFile, "\n";

##############################
### Go through each base in order and make consensus

for my $baseNum ( sort {$a <=> $b} keys %alignmentHash) {                       # go through each base in alignment in order
    my @baseArray = ("A", "T", "G", "C");                                       # this has to be in this order to match %ambiguityHash

    my $missing = 0;
    if(exists(($alignmentHash{$baseNum}{"-"}))) {
      $missing = ($alignmentHash{$baseNum}{"-"} / $seqCount) * 100;
    }

    if($missing < $maxMissing) {
      my $topBaseString = "";
      for my $base (@baseArray) {
        if(exists($alignmentHash{$baseNum}{$base})) {
          if(($alignmentHash{$baseNum}{$base} / $seqCount) * 100 >= $minAf) {
            $topBaseString .= $base;
          }
        }
      }
      if($topBaseString eq "") {
        print "N";
        $tooLowCount++;
      } else {
        print $ambiguityHash{$topBaseString};
      }
    }
}
print "\n";

if($verbose) {
    print STDERR "\n", $tooLowCount," bases were below --minAf $minAf\n"
}




##############################
# POD
##############################

#=pod

=head SYNOPSIS

Summary:

  genConsensusFromFastaAlignment.pl - generates a consensus for a specified gene in a specified taxa

Usage:

  perl genConsensusFromFastaAlignment.pl [options] --inFile exampleFasta.fasta > outputConsensus.fasta

=head OPTIONS

=over 4

=item B<--verbose>

        Output status information as the program runs.

=item B<--help>

        Print out program information

=item B<--inFile>

        Input fasta. Must be sequential, not interleaved.

=item B<--minAf> (0)

        Minimum percent of sequences that have a minor variant for ambiguous base call.

=item B<--maxMissing> (50)

        Maximum percent of sequences with a gap for the position to be included in consensus.

=cut
