#!/usr/bin/perl -w
####### This version of MD does not use PROFcon and instead uses ucon prediction based on statustical potential alone
use Cwd;
use File::Copy;
use Carp qw(cluck);
print_help()                 if ($#ARGV<0 || $ARGV[0]=~/^(-h|help|\?|def)$/i);

foreach $arg (@ARGV){
        #if    ($arg=~/^blast=(.*)$/)             { $fileInBlast =        $1;}
        if    ($arg=~/^fasta=(.*)$/)             { $file =        $1;}
	elsif ($arg=~/^disopred=(.*)$/i)         { $diso    =        $1;} 
	elsif ($arg=~/^hssp=(.*)$/)               { $fileHssp  =        $1;} 
	elsif ($arg=~/^prof=(.*)$/)               { $fileProf  =        $1;}     
	elsif ($arg=~/^profbval_raw=(.*)$/)                    { $fileProfbval  =         $1;}
	elsif ($arg=~/^norsnet=(.*)$/)                  { $fileNors      =         $1;}
	elsif ($arg=~/^ucon=(.*)$/)               { $fileUcon =        $1;}
	elsif ($arg=~/^chk=(.*)$/)               { $fileChk =        $1;}
	elsif ($arg=~/^md_dir=(.*)$/)               { $root_dir =        $1;} ## md directory will be called root from now on. tmp and bin directories will be within root
	elsif ($arg=~/^servers=(.*)$/)               { $print_servers =        $1;}
	elsif ($arg=~/^blast_dir=(.*)$/)               { $dirBlast =        $1;}
	elsif ($arg=~/^out=(.*)$/)               { $fileOut =        $1;}
	elsif ($arg=~/^out_mode=(.*)$/)               { $outMode =        $1;}
	elsif ($arg=~/^debug=(.*)$/)                   { $dbg    =         $1;}
	#elsif ($arg=~/^current_dir=(.*)$/)               { $path =        $1;}
	#elsif (-e $arg)                          { push(@fileIn,$arg);      }
	else {
	    print "*** wrong command line arg '$arg'\n";
	    die;
	}
    }
if (!defined$file) {
	die "cant run without a sequence";
	}
@orla=split(/\//,$file) ;

if (!defined$root_dir) {
    $root=$ENV{PP_MD}||$ENV{PP_PUB}."/md";  
}
else {
    $root=$root_dir; 
} 
if (!defined$dirBlast) {
    $blast_dir=$ENV{PP_BLASTROOT};    

        }
else {
        $blast_dir=$dirBlast;
        }
if (!defined$outMode) {
        $mode_out=0;
        }
else {
        $mode_out=$outMode;
        }

#$fileIn=$fileIn[0];
#die ("missing input $fileIn\n") if (! -e $fileIn);
$mcf=pop(@orla);
if ($mcf=~/(.*)\..*/) {
	$fileroot=$1;
	}
else {
	$fileroot=$mcf;
	}
#$path= "/" . join("/",@orla) . "/";
$win=31;
$win_acc=1;
$win_sec=1;
$ON= 2;
#$tmp=$root . "/tmp";
#$work_dir=".";
#$work_dir=$root;#for server
$rand_num=int(rand(10000000)); $rand_num= "tmpMD".$rand_num;
$id= $rand_num;


if (-d $ENV{PP_WORK}){
	 $work_dir=$ENV{PP_WORK};
}elsif (-d "/dev/shm"){ 
    $work_dir="/dev/shm/";
}else{ 
    $work_dir="/tmp/";
}


$work_dir.= "md/$id";
#die $work_dir;
if (! (-d  $work_dir)){
	warn ("mkdir -p $work_dir") if ($dbg);  
	system ("mkdir", "-p", "$work_dir")==0 or croak ("mkdir -p $work_dir: $?");
} 


$nn=$root . "/nn_files";

$cutoff=0.52;
#$id="3cat"; ############### for testing
$fasta="$work_dir/$id.f";

warn ("cp $file $fasta") if ($dbg);
system ("cp $file $fasta");

#$log_file=$tmp . "/" . $id . ".log";
$log_file=$fileroot . ".log";
###files:
$disopred_file="$work_dir/".$id.".diso";
$profbval_file="$work_dir/".$id."-63-R1-209-237433-ON2";
$ucon_file="$work_dir/".$id.".ucon";
$norsnet_file="$work_dir/$id.norsnet";
$hsspfil_file="$work_dir/$id-fil.hssp";
$prof_file="$work_dir/$id-fil.rdbProf";
$chk_file="$work_dir/$id.chk";

open (LOG, ">$log_file") || die "cant open log file $!";
if ((!defined$fileHssp)||(!defined$fileChk)) {
	print LOG "######### running blast for $id; creating new HSSP and chk files########\n";
	$run_psi= "perl $root/run_psi.pl";
	system ("$run_psi $id $work_dir $blast_dir");
	}
else {
	error_file($fileHssp);
	print LOG "hssp file exists copying $fileHssp $hsspfil_file\n";
	system ("cp $fileHssp $hsspfil_file");
        print LOG "chk file exists copying $fileChk $chk_file\n";
        system ("cp $fileChk $chk_file");
	}
if (!defined$fileProf) {
	print LOG " ######### running ucon-noProfcon for $id ########\n";
	$run_prof= "perl $root/run_prof.pl";
	system ("$run_prof $id $work_dir");
	}
else {
	error_file($fileProf);	
	print LOG "prof exists copying $fileProf $prof_file\n";
	system ("cp $fileProf $prof_file");
	}
if (!defined$diso) {
	print LOG " ######### running disopred for $id ########\n";
	$run_disopred= "perl $root/run_disopred.pl";

	warn ("$run_disopred $id $work_dir $root $blast_dir") if ($dbg);
	system ("$run_disopred $id $work_dir $root $blast_dir");

	}
else {
	error_file($diso);	
	print LOG "disopred file exists copying $diso $disopred_file\n";
	system ("cp $diso $disopred_file");
	}

system ("perl $root/createDataFileAlign.pl $id $hsspfil_file $prof_file $work_dir $root");

if (!defined$fileProfbval) {
	print LOG " ######### running profbval for $id ########\n";
	$run_profbval= "perl $root/run_profbval.pl";
	system ("$run_profbval $id $work_dir $work_dir $root");
	}
else {
	error_file($fileProfbval);
	print LOG "profbval file exists copying $fileProfbval $profbval_file\n";
	system ("cp $fileProfbval $profbval_file");
	}
if (!defined$fileNors) {	
	print LOG " ######### running norsnet for $id ########\n";
	$run_norsnet= "perl $root/run_norsnet.pl";
	system ("$run_norsnet $id $work_dir $work_dir $root");
	}
else {
	error_file($fileNors);
	print LOG "norsnet file exists copying $fileNors $norsnet_file\n";
	system ("cp $fileNors $norsnet_file");
	}
if (!defined$fileUcon) {
	print LOG " ######### running ucon-noProfcon for $id ########\n";
	$run_uconNoProfcon= "perl $root/run_ucon-noProfcon.pl";
	system ("$run_uconNoProfcon $id $work_dir $root");
	}
else {
	error_file($fileUcon);
	print LOG " ucon-noProfcon exists copying $fileUcon $ucon_file\n";
	system ("cp $fileUcon $ucon_file");
	}
close (LOG);
#copy_server_files();
print "### getting all the data from DISOPRED2 #####\n";
get_disopred($disopred_file);

print "### getting all the data from PROFbval #####\n";
get_profbval($profbval_file);
print "### getting all the data from NORSnet #####\n";
get_nors($norsnet_file);
print "### getting all the data from Ucon #####\n";
get_ucon_random($ucon_file);


$IN=$win*21+9+4*$win_acc+4*$win_sec+21+4+3+2+1;@aa=('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y');
$inValid="$work_dir/in_$rand_num";
$inOutValid="$work_dir/in_out_$rand_num";
$parvalid="$work_dir/par_$rand_num";
$validOutFile="$work_dir/val_out_$rand_num";

# we copy the junction file to the work directory to avoid a long path name which casues tha fortrna netutal net to break
$jct_orig="$root/jct110-nu-w31-set5-ac1-sec1-p6f";
$jct = $jct_dest="$work_dir/jct.in";

warn ("cp $jct_orig $jct_dest") if ($dbg);
system ("cp $jct_orig $jct_dest");

undef @pred;undef @result1; undef @result2;undef @md_rel;
$jct_crp=$work_dir . "/jct_crap".$rand_num;
$NNo_tst_err= "NNo_tst_err.dat".$rand_num;
$NNo_yeah= $work_dir . "/NNo-yeah$rand_num.tmp";
$h=1;
$file="$work_dir/$id.data";
;undef @res;undef $end;undef@PREL, undef @otL; undef @otE;undef @otH;undef @RI_A;undef @RI_S;undef @RI_S;undef $expCon1;undef $expCon2;
undef @secC;  $lengthA=$lengthB=$lengthC=0;undef $hydroNet;undef @A;undef @C;undef @D;undef @E;undef @F;undef @G;undef @H;undef @I;undef @K;undef @L;undef @M;undef @N;undef @P;undef @Q;undef @R;undef @S;undef @T;undef @V;undef @W;undef @Y;
foreach $aa (@aa) {
	$c{$aa}=0;
	}
open (FILE, "$file") || die "tinofet $root/DataFiles/$file $!";
<FILE>;
while ($line=<FILE>)  {
	@stuff=split(' ', $line);
	#$resNum=$stuff[0];
	$A=$stuff[1];$C=$stuff[2];$D=$stuff[3];$E=$stuff[4];$F=$stuff[5];$G=$stuff[6];$H=$stuff[7];$I=$stuff[8];
	$K=$stuff[9];$L=$stuff[10];$M=$stuff[11];$N=$stuff[12];$P=$stuff[13];$Q=$stuff[14];$R=$stuff[15];$S=$stuff[16];$expCon1=$stuff[27];$expCon2=$stuff[28];
	$T=$stuff[17];$W=$stuff[18];$Y=$stuff[19];$V=$stuff[20];			
	$otH=$stuff[22];$otE=$stuff[23];$otL=$stuff[24];$PREL=$stuff[25];$RI_A=$stuff[26];
	$lengthA=$stuff[29];$lengthB=$stuff[30];$lengthC=$stuff[31];$res=$stuff[33];
        push (@res,$res); push (@A,$A) ;push (@C,$C);push (@D,$D);push (@E,$E);push (@F,$F);push(@G,$G);push(@H,$H);push(@I,$I);push(@K,$K);push(@L,$L);push (@M,$M);push(@N,$N);push(@P,$P);push(@Q,$Q);push(@R,$R);
	push(@S,$S);push(@T,$T);push(@V,$V);push(@W,$W);push(@Y,$Y);push (@otH,$otH);push (@otE,$otE);push (@secC,$secC);push (@otL,$otL);push (@RI_A,$RI_A);push (@Bnew,$Bnew);push (@PREL,$PREL);;
	$c{$res}++;
	}
close(FILE);
foreach $aa (@aa) {
	$c{$aa}=$c{$aa}/(scalar@res)*100;
	use integer;
	$c{$aa}=$c{$aa}*1;
	no integer;
	}
$sampIn=scalar@A;
open (FOUT, ">$inValid") || die "error0";
printf FOUT "* overall: (A,T25,I8)\nNUMIN                 :      %3d\nNUMSAMFILE            :   %6d\n*",$IN,$sampIn;
print FOUT "\n* samples: count (A8,I8) NEWLINE 1..NUMIN (25I6)\n";
print "#####collecting all positive samples#######\n";
 $end=scalar@res-1;
 undef @compo_all;undef @length;undef @sec_cont;
#@compo_all=compo_all();
#@length=length_prot();
@length=($lengthA,$lengthB,$lengthC);
@sec_cont=sec_cont();
#now, I wish to take from all these sequences only the samples predicted to be in nors 	
loop20:	for ($i=0;$i<scalar@res;$i++) {		
	$k=0;undef @info;
	$lower=$i-($win-1)/2;
	$higher=$i+($win-1)/2;
	$lower_sec=$i-($win_sec-1)/2;
	$higher_sec=$i+($win_sec-1)/2;
	$lower_acc=$i-($win_acc-1)/2;
	$higher_acc=$i+($win_acc-1)/2;		
#secondary structure prediction information
	#push (@info, secondary($lower,$higher,$end));
#loop for solvent accessibility prediction information		
	push (@info, acc($lower_acc,$higher_acc,$end));
	push (@info, secondary($lower_sec,$higher_sec,$end));
	push (@info, compo_win($lower,$higher,$end)); ## 21 because there are the termini residues
	push (@info,sec_cont_win($lower,$higher,$end)); ## 4  because there are the termini residues
#	push (@info,@compo_all); ## 20
	push (@info, @sec_cont); ## 3
	push (@info, $expCon1,$expCon2); ## 2
	push (@info,@length); ## 3
	push (@info, ${$norsP{$id}}[$i], ${$nors1{$id}}[$i], ${$nors2{$id}}[$i]);
	push (@info,${$diso{$id}}[$i]);
	push (@info,${$profbval1{$id}}[$i], ${$profbval2{$id}}[$i]); 
	#push (@info,${$ucon{$id}}[$i]);
	#push (@info,${$ucon_only{$id}}[$i], ${$ucon_only_prob05{$id}}[$i]);
	push (@info,${$ucon_random{$id}}[$i]);
	push (@info,profiles($lower,$higher,$end));
	printf FOUT "ITSAM:%10d\n",$h; 
	presentIt(\@info);
	$h++;
	if ($h==($sampIn + 1))  {
		print FOUT "//";
		last loop20;}
	}
	 
print "\n########creating out-valid $inOutValid ######\n";
open (FOUT, ">$inOutValid") || die "cant open file $!";
printf FOUT "* overall: (A,T25,I8)\nNUMOUT                :        $ON\nNUMSAMFILE            :%9d\n*",$sampIn;
print FOUT "\n* samples: count (I8) SPACE 1..NUMOUT (25I6)\n";
for ($i=1; $i<$h;$i++)  {
	printf FOUT "%8d",$i;
      	printf FOUT "  %5d%6d\n",100,0;
	}
print FOUT "//"; 
print "\t$i\n";
print 1;
close (FOUT);
#next loop0;   ## to avoid recreating sample files!
open (FOUT, ">$parvalid")  || die "can't open file $!";
print FOUT "* I8\n";
printf FOUT "NUMIN                 :      %3d\n",$IN;
print FOUT "NUMHID                :        5\n";
print FOUT "NUMOUT                :        $ON\n";
print FOUT "NUMLAYERS             :        2\n";
printf FOUT "NUMSAM                :%9d\n",$sampIn;
print FOUT "NUMFILEIN_IN          :        1\n";
print FOUT "NUMFILEIN_OUT         :        1\n";
print FOUT "NUMFILEOUT_OUT        :        1\n";
print FOUT "NUMFILEOUT_JCT        :        1\n";
print FOUT "STPSWPMAX             :        0\n";
print FOUT "STPMAX                :        0\n";
print FOUT "STPINF                :        1\n";
print FOUT "ERRBINSTOP            :        0\n";
print FOUT "BITACC                :      100\n";
print FOUT "DICESEED              :   100025\n";
print FOUT "DICESEED_ADDJCT       :        0\n";
print FOUT "LOGI_RDPARWRT         :        1\n";
print FOUT "LOGI_RDINWRT          :        0\n";
print FOUT "LOGI_RDOUTWRT         :        0\n";
print FOUT "LOGI_RDJCTWRT         :        0\n";
print FOUT "* --------------------\n";
print FOUT "* F15.6\n";
print FOUT "EPSILON               :        0.010000\n";
print FOUT "ALPHA                 :        0.100000\n";
print FOUT "TEMPERATURE           :        1.000000\n";
print FOUT "ERRSTOP               :        0.000000\n";
print FOUT "ERRBIAS               :        0.000000\n";
print FOUT "ERRBINACC             :        0.200000\n";
print FOUT "THRESHOUT             :        0.500000\n";
print FOUT "DICEITRVL             :        0.100000\n";
print FOUT "* --------------------\n";
print FOUT "* A132\n";
print FOUT "TRNTYPE               : ONLINE\n";
print FOUT "TRGTYPE               : SIG\n";
print FOUT "ERRTYPE               : DELTASQ\n";
print FOUT "MODEPRED              : sec\n";
print FOUT "MODENET               : 1st,unbal\n";
print FOUT "MODEIN                : win=5,loc=aa\n";
print FOUT "MODEOUT               : KN\n";
print FOUT "MODEJOB               : mode_of_job\n";
print FOUT "FILEIN_IN             : $inValid\n";
print FOUT "FILEIN_OUT            : $inOutValid\n";
#print FOUT "FILEIN_JCT            : ./jct$iter-noucon-win$win-set$set-acc$win_acc-sec$win_sec-prop6f\n";
print FOUT "FILEIN_JCT            : $jct\n";
#print FOUT "FILEIN_JCT            : $dir/jct$i-test\n";
print FOUT "FILEOUT_OUT           : $validOutFile\n";
print FOUT "FILEOUT_JCT           : $jct_crp\n";
print FOUT "FILEOUT_ERR           : $NNo_tst_err\n";
print FOUT "FILEOUT_YEAH          : $NNo_yeah\n";
print FOUT "//\n";

close(FILE);
system ("$ENV{PP_BIN}/NetRun.LINUX $parvalid 1>>$log_file");

if ( -e  $validOutFile ) {
	open (LOG, ">>$log_file") || die "cant open log file $!";
    	print LOG "Meydey ran successfully. removing temporary files: $fasta $prof_file $hsspfil_file $file $inValid $inOutValid $parvalid $jct_crp $NNo_tst_err $NNo_yeah\n";
	close(LOG);
	if (defined$print_servers) {
		copy_server_files();
		}
	system ("rm $fasta $file $prof_file $hsspfil_file $inValid $inOutValid $parvalid $jct_crp $NNo_tst_err $NNo_yeah");
	print "parcing the results\n";
	get_md($validOutFile);
	get_rel();
	create_final();
        if (defined$fileOut) {
                system ("mv $fout $fileOut");
		print "\nfinal output is in $fileOut\n";
                }
	else {
		print "\nfinal output is in $fout\n";
		}
	system ("rm $log_file");
	
	}
else {
	print "$validOutFile doesnt exist!!!\n";
	}
unlink ($jct)||cluck("Cannot remove $jct");				
rmdir ($work_dir)||cluck("Cannot remove $work_dir");				
#=========================================================================================================================================================
sub print_help {
    if ($#ARGV<0 ||			# help
	$ARGV[0] =~/^(-h|help|special)/){
	print  "goal: predict proteon disorder using MD \n";
	print  "use:  run_MD-noPROFcon.pl fasta=(file_fasta) \n \n";
	print  "       will run many methods: psiblast,norsnet, profbval,disopred,ucon\n";
        print  "       can use the following files as input to skip some steps \n";
	print  "opt:  \n";
	printf "%5s %-15s=%-20s %-s\n","","hssp", "x",           "will run psiblast->covert to hssp->filter hssp using 80% reduction, unless specified ";	
	printf "%5s %-15s=%-20s %-s\n","","prof", "x",           "takes profRdb files; will run prof unless specified. "; 
	printf "%5s %-15s=%-20s %-s\n","","profbval_raw", "x",           "raw output from profbval";
	printf "%5s %-15s=%-20s %-s\n","","norsnet", "x",        "output file from norsnet (not raw)";	
	printf "%5s %-15s=%-20s %-s\n","","disopred", "x",           "disopred1 file casp format. runs disopred unless specified";
	#printf "%5s %-15s=%-20s %-s\n","","ucon", "x",           "in this version, the ucon file should not result from contact predictor, only from using statistical potential";
	printf "%5s %-15s=%-20s %-s\n","","out", "x",           "final output file";
	printf "%5s %-15s=%-20s %-s\n","","chk", "x",           "chk file from blast; used as input for disopred. will run blast unless specified";
	printf "%5s %-15s=%-20s %-s\n","","md_dir", "x",         "location of all md files and scripts";
	printf "%5s %-15s=%-20s %-s\n","","blast_dir", "x",           "location of blast";
	exit;
    }
}				

sub shuffle  {
	my($ref) = shift;
	my(@array) = @{$ref};
	my @temp;
	push (@temp,splice(@array,rand(@array),1))
		while @array;
	@array=@temp;
	return @array;
	}
sub generate_rand  {
	my($ref) = shift;
	my(@array) = @{$ref};
	my $x;
	$x= $array[int(rand(scalar(@array)))];
	return $x;
	}
sub presentIt  {
        my($ref) = shift;
        my(@residue) = @{$ref};
	my $l;
	for ($l=0; $l<scalar(@residue); $l++)  {
		$residue[$l]=~s/\s//;
		if (defined($residue[$l])==0)
		 {
		 	print "\n$file\t$i\t$j\t$l";
			}
#		if ($residue[$l]=~/[a-z]|A-Z]/ )  {
#			print "\n$residue[$l]\t$file\t$i\t$j\t$l";
#			}
		printf FOUT "%6d",$residue[$l];
		$k++;
		if ($k==25)  {
			print FOUT "\n";
			$k=0;
			}
		}
	print FOUT "\n"; $k=0;
	
	return;
	}
#functions to obtain properties    
#profiles
sub profiles  {
	my $lower=shift;
	my $higher=shift;
	my $end= shift;
	my @residue; 
	my @array;
	for ($j=$lower; $j<=$higher; $j++)  {
		if (($j<0) ||($j>$end)) {
			undef @residue;
			@residue=(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,100);
			push (@array, @residue);
			}
		else {
			undef @residue;
			@residue=($A[$j],$C[$j],$D[$j],$E[$j],$F[$j],$G[$j],$H[$j],$I[$j],$K[$j],$L[$j],$M[$j],$N[$j],$P[$j],$Q[$j],$R[$j],$S[$j],$T[$j],$W[$j],$Y[$j],$V[$j],0);
			push (@array, @residue);
			}
		}
		return @array;
	}
#secondary structure prediction information
sub secondary  {
	my $lower=shift;
	my $higher=shift;
	my $end= shift;
	my @secon;
	my @array;
	for ($j=$lower; $j<=$higher; $j++)  {
		if (($j<0) ||($j>$end)) {
			undef @secon;
			@secon=(0,0,0,100);
			push (@array, @secon);
			}
		else {
			undef @secon;
			@secon=($otL[$j],$otH[$j],$otE[$j],0);
			push (@array, @secon);
			}
		}
	return @array;
	}
#function for solvent accessibility prediction information		
sub acc  {
	my $lower=shift;
	my $higher=shift;
	my $end= shift;
	my @PRE;
	my @array;
	for ($j=$lower; $j<=$higher; $j++)  {
		if (($j<0) ||($j>$end)) {
			undef @PRE;
			@PRE=(0,0,100,100);
			push (@array, @PRE);
			}
		else {
			undef @PRE;
			@PRE=($PREL[$j],$RI_A[$j],0,0);
			push (@array, @PRE);		
			}
		}
	return @array;
	}

sub max
{
   my $number1=shift;
    my $number2=shift;  
    my $number3=shift; 
    my(@numbers);
    push (@numbers,$number1, $number2,  $number3);
    my ($max);
    $max = $numbers[0];
    foreach my $i (@numbers)
    {
        if ($i >=$max)
        {
            $max = $i;
        }
    }
    return ($max);
}				
sub hydroNet  {
	my $end= shift;
	my $i; my @res_ww;my $count=0; my $sum=0; my $avg; my $hydroNet; my $charge=0;
	my $window; my $length=scalar@res; my $hydro=0;
loop1001:for ($i=0;$i<scalar@res;$i++) {
		$window=5; ## because in loop1010 whenever there is an 'X' i deduct 1 from the $window
		$sum=0;
		if (($res[$i] eq 'E') || ($res[$i] eq 'D') )  {
			$charge--;
			}
		elsif (($res[$i] eq 'K') || ($res[$i] eq 'R') ) {
			$charge++;
			}
		elsif ($res[$i] eq 'X') {
			$length--;
			}			
		my $lower=$i-2;
		my $higher=$i+2;
loop1010:	for ($j=$lower; $j<=$higher; $j++)  { # the hydrophobicity per residue is calculated in window=5
			if (($j<0) ||($j>(scalar@res-1))) {
				next loop1001;
				}
			if ($res[$j] eq "A") {$res_ww[$j] = -1.8;} 
 			elsif ($res[$j] eq "L") {$res_ww[$j] = -3.8;}
 		 	elsif ($res[$j] eq "R") {$res_ww[$j] = 4.5;}
 		 	elsif ($res[$j] eq "K") {$res_ww[$j] = 3.9;}
 		 	elsif ($res[$j] eq "N") {$res_ww[$j] = 3.5;}
 		 	elsif ($res[$j] eq "M") {$res_ww[$j] = -1.9;}
 		 	elsif ($res[$j] eq "D") {$res_ww[$j] = 3.5;}
 			 elsif ($res[$j] eq "F") {$res_ww[$j] = -2.8;}
 			 elsif ($res[$j] eq "C") {$res_ww[$j] = -2.5;}
 			 elsif ($res[$j] eq "P") {$res_ww[$j] = 1.6;}
 			 elsif ($res[$j] eq "Q") {$res_ww[$j] = 3.5;}
			  elsif ($res[$j] eq "S") {$res_ww[$j] = 0.8;}
			  elsif ($res[$j] eq "E") {$res_ww[$j] = 3.5;}
			  elsif ($res[$j] eq "T") {$res_ww[$j] = 0.7;}
			  elsif ($res[$j] eq "G") {$res_ww[$j] = 0.4;}
			  elsif ($res[$j] eq "W") {$res_ww[$j] = 0.9;}
			  elsif ($res[$j] eq "H") {$res_ww[$j] = 3.2;}
			  elsif ($res[$j] eq "Y") {$res_ww[$j] = 1.3;}
			  elsif ($res[$j] eq "I") {$res_ww[$j] = -4.5;}
			  elsif ($res[$j] eq "V") {$res_ww[$j] = -4.2;}
			  else{
			  	$window--; ##not to count X,
				next loop1010;
				}
			  $sum=$sum+($res_ww[$j]+4.5)/9; ## addition 0f 4.5 to make it 0-9 and then dividing by 9 to make it a fraction
			}
		if ($window==0)  {
				next loop1001;
				}
		$hydro=$hydro + $sum/$window; # $hydro is the hydrophobicity of the whole protein
		}
		$netCharge=$charge;
		if ($netCharge<0)  {
			$netCharge=-$netCharge;
			}
		elsif ($netCharge==0)  {
			return 100;
			}
		$hydroNet=$hydro/$netCharge*100/$length;;
		use integer;
		$hydroNet=$hydroNet*1;
		no integer;
		return $hydroNet;
	}
sub normDiff {
	my $totalDiff=shift;
        my($ref) = shift;
        my(@diff) = @{$ref};
	my $l=scalar@diff;
	my $k;my $sum=0;my $sigma;my @diffNew;
	my $avg=$totalDiff/$l;
	for ($k=0;$k<scalar@diff;$k++)  {
		$sum= $sum + ($diff[$k]- $avg)*($diff[$k]- $avg);
		}
	$sigma= sqrt($sum/($l-1));
	for ($k=0;$k<scalar@diff;$k++)  {
		$diffNew[$k]=($diff[$k]-$avg)/$sigma;
		}	
	return (@diffNew);
	}
	
sub getDiff  {
	my $lower=shift;
	my $higher=shift;
	my $end= shift;
	my @diff100;
	my @array;
	for ($j=$lower; $j<=$higher; $j++)  {
		undef @diff100;
		if (($j<0) ||($j>$end)) {
			@diff100=(0,0,0,100);
			}
		else {	
			if ($diffNew[$i]>1.2) {
				@diff100=(100,0,$diff11[$j],0);
				}
			elsif (($diffNew[$i]<1.2) && ($diffNew[$i]>-1.3)) {
				@diff100=(50,50,$diff11[$j],0);
				}
			else {
				@diff100=(0,100,$diff11[$j],0);
				}
				
			}
		push (@array,@diff100);
		}
	return @array;
	}
sub get_nors{
	my @moogla;
	undef @moogla;my $nors_score;$file=shift;my $nors_raw;
	if (!(open (FILE, $file))) {
		print "cant open $file";;
		rm_files($id);
		}
	<FILE>;
	while ($line=<FILE>) {
		chomp$line;
		@moogla=split(/\t/,$line);
		$nors_score=$moogla[4]*100;
		$nors_raw=$moogla[4];
		if ($nors_raw>=0.52) {
			push (@{$nors2states{$id}},'D');
			}
		else {
			push (@{$nors2states{$id}},'-');
			}
		push (@{$norsP{$id}},$nors_score);
		push (@{$nors_raw{$id}},$nors_raw);
		push (@{$nors1{$id}},$moogla[2]);
		push (@{$nors2{$id}},$moogla[3]);
		}
	close (FILE);
	}
	
sub get_disopred {
	my @moogla;
	my $file= shift;my $diso_score;my $diso2states;my $diso_raw; 
       if (!(open (FILE, $file) )){
       		print  "cant open $file . killing the job $!";
       		rm_files($id);
		die
		}
        for ($x=0;$x<6;$x++)  {
                <FILE>;
                }
loop59: while ($line=<FILE>) {
                undef @moogla;$line=~ s/\n//;
                @moogla=split (' ',$line);
                if ($moogla[0] eq 'END')  {
                        last loop59;
			close (FILE)
                        }
		$diso_score=$moogla[2]*100;
		$diso_raw=$moogla[2];
		if ($diso_raw>=0.55) {
			push (@{$diso2states{$id}},'D');
			}
		else {
			push (@{$diso2states{$id}},'o');
			}
                push (@{$diso{$id}}, $diso_score);
                }
	}
sub get_profbval {
	my @moogla;my $file=shift;	
	my $num;my $node1;my $node2;my $muka;my $twostate;
       if (!(open (FILE, "$file"))) {
       		print "cant open $file killing the job$!";
		rm_files($id);
		die;
		}
	for ($f=0; $f<43; $f++)  {
   		 <FILE>;
       		 }	
loop58:	while ($line=<FILE>)  {
        	if ($line=~ /^(.{8})(.{5})(.{4})/)  {
                	$num=$1;$node1=$2; $node2=$3;
			$node1=~ s/\s//g;$node2=~ s/\s//g;$twostate='-';
			$muka=$node1/($node1+$node2);
			push (@{$profbval1{$id}},$node1);
			push (@{$profbval2{$id}},$node2);
			$muka=($node1-$node2+100)/200;
			 push (@{$profbval_prob{$id}},$muka);
			if ($muka>=0.465) {
				$twostate='D';
				}
			push (@{$profbval2_states{$id}},$twostate);		
			}
		}
	close (FILE);
	}
sub get_ucon {
	my @moogla;my $file; my $num;my $prob;my $val; my $line;
	$file=$id.".ucon";
       	open (FILE, "$root/ucon/$file") || die "$root/ucon/$file no siyut ahhhhhhhhhh $!";
	<FILE>;
loop57:	while ($line=<FILE>)  {
		chomp $line;
		if ($line=~/END/) {
			close(FILE);
			last;
			}
		@moogla=split (/\t/,$line);
		$val=100*$moogla[3];
		#push (@{$ucon{$id}}, $val);###### just to avoid warning messages for now
		}
	close(FILE);
	}
sub get_ucon_only {
	my @moogla;my $file; my $num;my $prob;my $val; my $line;
	$file=$id.".ucon";
       	open (FILE, "$root/ucon_only/$file") || die "$root/ucon_only/$file no siyut ahhhhhhhhhh $!";
	<FILE>;
loop67:	while ($line=<FILE>)  {
		if ($line=~/END/) {
			close(FILE);
			last;
			}
		chomp $line;
		@moogla=split (/\t/,$line);
		$val=$moogla[4];
		#push (@{$ucon_only{$id}}, $val);###### just to avoid warning messages for now
		}
	close(FILE);
	}
sub get_ucon_only_prob05 {
	my @moogla;my $file; my $num;my $prob;my $val; my $line;
	$file=$id.".ucon";
       	open (FILE, "$root/ucon_only_prob05/$file") || die "$root/ucon_only_prob05/$file no siyut ahhhhhhhhhh $!";
	<FILE>;
loop77:	while ($line=<FILE>)  {
		if ($line=~/END/) {
			close(FILE);
			last;
			}
		chomp $line;
		@moogla=split (/\t/,$line);
		$val=$moogla[4];
		#push (@{$ucon_only_prob05{$id}}, $val);###### just to avoid warning messages for now
		}
	close(FILE);
	}		
sub compo_all {
        my $res;
        my @compo;my $aa;my %count;my %compo;
	my $length=scalar@res;my @array;
        my @aa=('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y');
        foreach $aa (@aa) {
                $count{$aa}=0;
                $compo{$aa}=0;
                }
        for ($j=0; $j<scalar@res; $j++)  {
                 $res=$res[$j];
	#	if (defined$count{$res}==0) {
#			print $res;
#			}
              	 $count{$res}++;
                }
        foreach $aa (@aa) {
                $compo{$aa}=int(($count{$aa}/$length*100));
                push (@array, $compo{$aa});
                }
        return @array;
        }
sub sec_cont {
        my $sec; my %sum;my %compo;my $j;
        my $length=scalar@otL;my @array;
        my @sec=('H','E','L');
        foreach $sec (@sec) {
                $sum{$sec}=0;
                $compo{$sec}=0;
                }
        for ($j=0; $j<scalar@otL; $j++)  {
                 $sum{'L'}= $sum{'L'}+$otL[$j];
		$sum{'E'}= $sum{'E'}+$otE[$j];
		$sum{'H'}= $sum{'H'}+$otH[$j];
                }
        foreach $sec (@sec) {
                $compo{$sec}=int(($sum{$sec}/$length));
                push (@array, $compo{$sec});
                }
        return @array;
        }
sub compo_win  {
        my $lower=shift;
        my $higher=shift;
        my $end= shift;my $res;
        my $aa;my %count;my %compo; my @array;
	my @aa=('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','Z');
	foreach $aa (@aa) {
		$count{$aa}=0;
		$compo{$aa}=0;
		}
        for ($j=$lower; $j<=$higher; $j++)  {
                if (($j<0) ||($j>$end)) {
			$res='Z';
                        }
                else {

			$res=$res[$j];
			}
               $count{$res}++;
		}
	foreach $aa (@aa) {
		$compo{$aa}=int(($count{$aa}/$win*100));
		push (@array, $compo{$aa});		
                }
        return @array;
        }
sub sec_cont_win  {
        my $lower=shift;
        my $higher=shift;
        my $end= shift;my $res;
        my $sec;my %count;my %compo; my @array;my $count=0;;my $j;my %sum;
        my @sec=('H','E','L');
        foreach $sec (@sec) {
                $sum{$sec}=0;
                $compo{$sec}=0;
                }
        for ($j=$lower; $j<=$higher; $j++)  {
                if (($j<0) ||($j>$end)) {
			{
				@array=(0,0,0,100);
				return @array;
				}
			}
                else {
               		 $sum{'L'}= $sum{'L'}+$otL[$j];
			$sum{'E'}= $sum{'E'}+$otE[$j];
			$sum{'H'}= $sum{'H'}+$otH[$j];
	      		$count++;
			}
		}
        foreach $sec (@sec) {
                $compo{$sec}=int(($sum{$sec}/$count));
                push (@array, $compo{$sec});
                }
	push (@array,0);
        return @array;
        }
sub get_md {
my $result1;my $result2;my $w;
my $file=shift;
open (F, "$file") ||  die "cant open $file, wasnt supposed to be here then...\n";
for ($r=0; $r<43; $r++) {
        <F>;
        }
while ($line=<F>) {
        if ($line=~ /^(.{8})(.{5})(.{4})/)  {
                $result2=$3;
                $result1=$2;
                $result2=~ s/\s//g;$result1=~ s/\s//g;
                push (@result1,$result1);push (@result2,$result2);
                $w=$result1/($result1+$result2);
		if ($w>=$cutoff) {
			push (@{$md2states{$id}},'D');
			}
		else {
			push (@{$md2states{$id}},'-');
			}
                push (@pred, $w);
		}
	}
close(F);
#system ("rm $validOutFile_orig");
}
sub get_ucon_random {
	my @moogla;my $file=shift;my $muka; 
	my $num;my $prob;my $val; my $line;
       	if (!(open (FILE, "$file"))){
		 print "$file no siyut ahhhhhhhhhh $!";
		 rm_files($id);
		 }
	<FILE>;
loop57:	while ($line=<FILE>)  {
		chomp $line;
		if ($line=~/END/) {
			close(FILE);
			last;
			}
		@moogla=split (/\t/,$line);
		$val=100*$moogla[3];
		push (@{$ucon_random{$id}}, $val);
		push (@{$ucon_prob{$id}}, $moogla[3]);
		if ($moogla[3]>=0.58) {
			$muka='D';
			}
		else {
			$muka='-';
			}
		push (@{$ucon2_states{$id}}, $muka);
		}
	close(FILE);
	}
sub rm_files {
	my $id=shift;
	my $fold=$root;
	my $foldTmp=$work_dir;
	open (LOG, ">>$log_file") || die "cant open log file $!";
#	print LOG "rm $fold/*$id* \n";
	print LOG "rm $foldTmp/*$id*\n";system ("rm $foldTmp/*$id*"); ###the rm line
#	print LOG "$foldSer/*$id*\n";
	close (LOG);
#	system ("rm $fold/$id*");system ("rm $foldSer/$id*");
	die;
	}
sub error_file {
	my $file=shift;
	if (!( -e $file  )) {
		print LOG "$file does not exist. possible typo? \n";
		print LOG "removing files of $id \n";
		rm_files($id);
		die;
		}
	}	
sub create_final {
	$fout= "$root/" . "$fileroot.md";
	#$fout= "./" . "$fileroot.md";
	open (FOUT, ">$fout") || die "cant open  $fout";
	if ($mode_out==0) {
#		print FOUT "num\taa\tnorsnet\tnorsnet\tdiso\tdiso\tbval1\tbval2\tucon\tnode1\tnode2\tmd\tmd\n";
       		print FOUT "num\taa\tnode1\tnode2\tscore\tmd\n";
		for ($i=0;$i<scalar@res;$i++) {
			$num=$i+1;
		#	printf FOUT "%4d\t%s\t%1.2f\t%s\t%1.2f\t%s\t%d\t%3d\t%3d\t%3d\t%3d\t%1.3f\t%s\n",$num,$res[$i],${$norsP{$id}}[$i],${$nors2states{$id}}[$i],${$diso{$id}}[$i],${$diso2states{$id}}[$i],${$profbval1{$id}}[$i], ${$profbval2{$id}}[$i],0,$result1[$i],$result2[$i],$pred[$i],${$md2states{$id}}[$i];
			printf FOUT "%4d\t%s\t%3d\t%3d\t%1.3f\t%s\n",$num,$res[$i],$result1[$i],$result2[$i],$pred[$i],${$md2states{$id}}[$i];
			}
		close(FOUT);
		print "raw results are in $validOutFile\n";
		open (LOG, ">>$log_file") || die "cant open log file $!";
		print LOG "$work_dir/*$id*\n";
    		print LOG "removing temporary result file:$validOutFile \n";
	
		}
	elsif($mode_out==1) {
		#print FOUT "num\taa\tNORSnet\tNORSnet2st\tPROFbval\tPROFbval2st\tUcon\tUcon2st\tMD\tMD2st\n";
                printf FOUT '%s %s %s %s %s %s %s %s %s %s %s',"Number","Residue","NORSnet","NORS2st","PROFbval","bval2st","Ucon","Ucon2st","MD_raw "," MD_rel"," MD2st ";
		print FOUT "\n";
		for ($i=0;$i<scalar@res;$i++) {
                        $num=$i+1;
                       printf FOUT "%5d\t%s\t%1.2f\t%s\t%1.2f\t%s\t%1.2f\t%s\t%1.3f\t%d\t%s\n",$num,$res[$i],${$nors_raw{$id}}[$i],${$nors2states{$id}}[$i],${$profbval_prob{$id}}[$i], ${$profbval2_states{$id}}[$i],${$ucon_prob{$id}}[$i],${$ucon2_states{$id}}[$i],$pred[$i],$md_rel[$i],${$md2states{$id}}[$i];
                        }
		print FOUT "\n\nKey for output\n";
		print FOUT "----------------\n";
		print FOUT "Number - residue number\n";
		print FOUT "Residue - amino-acid type\n";
		print FOUT "NORSnet - raw score by NORSnet (prediction of unstructured loops)\n";
		print FOUT "NORS2st - two-state prediction by NORSnet; D=disordered\n";
		print FOUT "PROFbval - raw score by PROFbval (prediction of residue flexibility from sequence)\n";
		print FOUT "Bval2st - two-state prediction by PROFbval\n";
		print FOUT "Ucon - raw score by Ucon (prediction of protein disorder using predicted internal contacts)\n";
		print FOUT "Ucon2st - two-state prediction by Ucon\n";
		print FOUT "MD - raw score by MD (prediction of protein disorder using orthogonal sources)\n";
		print FOUT "MD2st - two-state prediction by MD\n";
		print FOUT "MD_rel - reliability of the prediction by MD; values range from 0-9. 9=strong prediction\n";
                close(FOUT);
                print "raw results are in $validOutFile\n";
                open (LOG, ">>$log_file") || die "cant open log file $!";
                print LOG "$work_dir/*$id*\n";
                print LOG "removing temporary result file:$validOutFile \n";
		}

	close (LOG);
#	rm_files($id);
	system ("rm $work_dir/$id*");print "rm $work_dir/$id*\n";####another rm line
	system ("rm $validOutFile");
	}	
sub copy_server_files {
	$disopred_final="$work_dir/".$fileroot.".diso";
	$ucon_final="$work_dir/".$fileroot.".ucon";
	$norsnet_final="$work_dir/$fileroot.norsnet";
	#$hsspfil_final="$root/$id-fil.hssp";
	$prof_final="$work_dir/$fileroot-fil.rdbProf";
	system ("cp $disopred_file $disopred_final");
	system ("cp $ucon_file $ucon_final");
	system ("cp $norsnet_file $norsnet_final");
	system ("cp $prof_file $prof_final");
	}
#sub rm_server_files {
#	$disopred_final="$work_dir/".$fileroot.".diso";
#	$ucon_final="$work_dir/".$fileroot.".ucon";
#	$norsnet_final="$work_dir/$fileroot.norsnet";
#	#$hsspfil_final="$root/$id-fil.hssp";
#	$prof_final="$work_dir/$fileroot-fil.rdbProf";
#	system ("rm $disopred_final");
#	system ("rm $ucon_final");
#	system ("rm $norsnet_final");
#	system ("rm $prof_final");
#	}
sub get_rel {
	my $i;my $tmp; my $rel;my %iter;my %min;my %max; 
	my $set=5;my $result1;my $result2;
	$iter{1}=696; $iter{2}=68; $iter{3}=106; $iter{4}=84; $iter{5}=110; $iter{6}=37;
	$min{1}=0.14;$max{1}=0.92; $min{2}=0.24;$max{2}=0.72; $min{3}=0.22;$max{3}=0.76;
	 $min{4}=0.26;$max{4}=0.76; $min{5}=0.22;$max{5}=0.76; $min{6}=0.30;$max{6}=0.72;
	for ($i=0;$i<scalar@result1;$i++) {
		$result1=$result1[$i];$result2=$result2[$i];
                #$tmp=int(100*($result1/($result1+$result2)));  ## In this version i calculate the reliability relative to the actual score. luckily the min/max values are the same so i use the same liniar function
               $tmp=$result1/($result1+$result2);		#the same as the line above but i had to divide by 100 
		 $score=$tmp;
                if ($tmp>$max{$set}) {$rel=9;}
                elsif ($tmp<$min{$set}) {$rel=0;}
                else {
                       if ($tmp>=$cutoff) {
                              $rel=int(10*($tmp-$cutoff)/($max{$set}-$cutoff));
                              }
                        else {
                              $rel=int(10*($cutoff-$tmp)/($cutoff-$min{$set}));
                             }
                       }
                  push (@md_rel, $rel);
		}
	}
