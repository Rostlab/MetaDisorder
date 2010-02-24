#!/usr/bin/perl -w
##!/usr/local/bin/perl -w
if (@ARGV<3)  {
	die "\nUsage: $0 perl run_ucon_noProfcon  [id] [work_dir] [root]\n";
	}
$id=$ARGV[0];
$work_dir=$ARGV[1];
$root=$ARGV[2];
$temp="$work_dir/$id.eprofcon";
#$servers="$root/servers";
system ("perl $root/ucon-random-MJ-UR/calc_energy_1seq.pl 60 $id $work_dir $root");
system ("perl $root/ucon-random-MJ-UR/smooth_1seq.pl $id 11 $work_dir");
system ("rm $temp");
