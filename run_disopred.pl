#!/usr/bin/perl -w
use Cwd;
use File::Copy;

#$tmp="$root/tmp/";
if (@ARGV<4)  {
        die "\nUsage: $0 perl run_disopred.pl  [id] [work_dir] [root] [blast_dir]\n";
        }
$work_dir=$ARGV[1]; ####### do not use it anymore. 
$fileroot=$ARGV[0];
$root=$ARGV[2];
$blast_dir=$ARGV[3];
#$file=$root . "/" . $fileroot . ".f";
$file=$work_dir. "/" . $fileroot . ".f";

#$file2= $root . "/" . $fileroot . ".diso";
#system ("$root/psipred/rundisopred-input-path2 $file $root");
$chk=$work_dir."/"."$fileroot.chk"; ### if check file is provided
#system ("$root/psipred/rundisopred-input-path2-chk $file $root $chk $blast_dir"); ### if check file is provided
#warn ("$root/psipred/rundisopred-input-path2-chk2 $file $root $chk $blast_dir $work_dir"); ### if check file is provided
system ("$root/psipred/rundisopred-input-path2-chk2 $file $root $chk $blast_dir $work_dir"); ### if check file is provided
#$diso="$root/" . $fileroot . ".diso";
#system ("mv $diso $file2");
