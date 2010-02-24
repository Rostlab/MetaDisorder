#!/usr/local/bin/perl -w
$fasta=$ARGV[0];
$prot_name= $ARGV[1];
$output_file=$ARGV[2];
$root="/home/schles/ucon-random-MJ-UR/";
system ("perl $root/calc_energy_1seq.pl 60 $fasta $prot_name");
system ("perl $root/smooth_1seq.pl $prot_name 11 $output_file");
