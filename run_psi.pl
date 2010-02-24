#!/usr/bin/perl
use Cwd;
use File::Copy;

if (@ARGV<3)  {
	die "\nUsage: $0 perl run_psi_prof  [id] [work] [blast_dir]\n";
	}
system ("export ARCH=LINUX");
$id=$ARGV[0];$root=$ARGV[1];$blast_dir=$ARGV[2];
$fasta=$root ."/". $id .".f";
$hsspfil=$root ."/".$id."-fil.hssp";
$file1= $root ."/"."prof/"."$id-fil.rdbProf";
$saf=$root ."/".$id.".saf";
$hssp=$root ."/".$id.".hssp";
$blastpgp=$root ."/".$id .".blastpgp";;
$chk=$root ."/".$id .".chk";;
warn ("$blast_dir/blastpgp -i $fasta -j 3 -d /data/blast/big -o $blastpgp -C $chk");
system ("$blast_dir/blastpgp -i $fasta -j 3 -d /data/blast/big -o $blastpgp -C $chk");

warn (" $ENV{'HOME'}/server/pub/prof/scr/blast2saf.pl $blastpgp maxAli=3000 eSaf=1 saf=$saf");
system ("perl $ENV{'HOME'}/server/pub/prof/scr/blast2saf.pl $blastpgp maxAli=3000 eSaf=1 saf=$saf");
warn  ("perl $ENV{'HOME'}/server/pub/prof/scr/copf.pl exeConvertSeq=$ENV{'HOME'}/server/pub/prof/bin/convert_seq_big.LINUX $saf hssp fileOut=$hssp");
system ("perl $ENV{'HOME'}/server/pub/prof/scr/copf.pl exeConvertSeq=$ENV{'HOME'}/server/pub/prof/bin/convert_seq_big.LINUX $saf hssp fileOut=$hssp");

system ("perl $ENV{'HOME'}/server/pub/prof/scr/hssp_filter.pl $hssp red=80 fileOut=$hsspfil");



# system ("perl /nfs/home1/pub/molbio/perl/blast2saf.pl $blastpgp maxAli=3000 eSaf=1 saf=$saf");
# system ("/usr/pub/molbio/prof/scr/copf.pl exeConvertSeq=/usr/pub/molbio/prof/bin/convert_seq_big.LINUX $saf hssp fileOut=$hssp");
# system ("/usr/pub/molbio/prof/scr/hssp_filter.pl $hssp red=80 fileOut=$hsspfil");
if ( -e $hsspfil ) {
	system ("rm $blastpgp $saf $hssp");
	}
else {
	#print LOG "hssp_filter.pl couldnt filter so it used the non_filtered file");
	system ("rm $blastpgp $saf");
	system ("mv $hssp $hsspfil");
	}


