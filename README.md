# genConsensusFromFastaAlignment.pl

The purpose of this script is to take in a DNA alignment in fasta format and
output a fasta consensus with ambiguity codes. The options determine how many
(what percent) sequences must contain an alternate base for a difference to be
considered "real". The input fasta must be sequential sequences and as of now
ambiguous bases in the input is unsupported (replaced with a gap).

## Usage

perl genConsensusFromFastaAlignment.pl [options] --inFile exampleFasta.fasta > outputConsensus.fasta

## Options

* --verbose

    Output status information as the program runs.

* --help

    Print out program information

* --inFile

    Input fasta. Must be sequential, not interleaved.

* --minAf(0)

    Minimum percent of sequences that have a minor variant for ambiguous base call.

* --maxMissing(50)

    Maximum percent of sequences with a gap for the position to be included in consensus.


## Input

An interleaved fasta.


## Output

A single fasta entry


## To do

* Maybe an option to output a table with base counts at each position.

* Maybe allow ambiguous bases in input fasta.
