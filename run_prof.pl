#!/usr/bin/perl
use Cwd;
use File::Copy;

if (@ARGV<2)  {
	die "\nUsage: $0 perl run_psi_prof  [id] [root]\n";
	}
$id=$ARGV[0];$root=$ARGV[1];
$fasta=$root ."/". $id .".f";
$hsspfil=$root ."/".$id."-fil.hssp";
$prof= $root . "/" . "$id-fil.rdbProf";
system ("/usr/pub/molbio/prof/scr/prof.pl $hsspfil fileOut=$prof");

