#!/usr/bin/perl -w
##!/usr/local/bin/perl -w
if (@ARGV<3)  {
	die "\nUsage: $0 perl run_ucon.pl  [id] [work_dir] [root]\n";
	}
$id=$ARGV[0];
$work_dir=$ARGV[1];
$root=$ARGV[2];
$temp="$work_dir/$id.eprofcon";
#$servers="$root/servers";
system ("perl $root/ucon_only/calc_energy_1seq.pl 200 $id $work_dir $root");
system ("perl $root/ucon_only/smooth_1seq.pl $id 11 $work_dir"); ### note that its 11 for ucon only. May 2008.
system ("rm $temp");
