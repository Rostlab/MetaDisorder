#!/usr/local/bin/perl -w
if (@ARGV<4)  {
	die "\nUsage: $0 [*-fil.hssp] [-fil.rdbprof] [profbval] [results_dir]\n";
	}
$hsspfil=$ARGV[0]; $prof=$ARGV[1];$res=$ARGV[3];$profbval=$ARGV[2];
$fileroot=$hsspfil;
@orla=split('/',$fileroot) ;
$mcf=$orla[$#orla];
if ($mcf=~ /(.*)-fil\.hssp/) {
	$fileroot=$1;
	}
elsif ($mcf=~/(.*)\..*/) {
	$fileroot=$1;
	}
$data= "$fileroot.data";
system ("perl ~schles/norsnet/createDataFileAlign.pl $hsspfil $prof");
print "\nfinished creating data file\n";
system ("perl ~schles/norsnet/norsnet01-randIds.pl 13 $data $profbval $res");
print ("finished running the netwrok. results are in $res\n");
