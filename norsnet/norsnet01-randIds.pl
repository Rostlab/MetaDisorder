#!/usr/local/bin/perl -w
#create the first in_test. the nodes' correspond to the amino acids in the following order:
#1-A,2-C,3-D,4-E,5-f,6-g,7-h,8-i,9-k,10-l,11-M,12-n,13-p,14-q,15-r,16-s,17-t,18-w,19-y,20-v,21-z
#the nodes of secondary structure are:
#1-G,2-H,3-I,4-E,5-B,6-T,7-S,8-L
use Cwd;
use File::Copy;
if (@ARGV<4)  {
	die "\nUsage: $0 perl profPred.pl [window] [data] [profbval] [results dir]\n";
	}
#$stretch=$ARGV[0];;
$win=$ARGV[0];
$file=$ARGV[1];
$jct="/home/schles/norsnet/jct50-5HN-win13-ComLenAcHNexNBba2rel13";
$dir="/home/schles/norsnet/work";;;
$resultsdir=$ARGV[3];
$profbval= $ARGV[2];
#$profbval= "$id-63-R1-209-237433-ON2";
#$thre=$ARGV[4];
$ON= 2;$sampIn=0;
$IN=416;$rand_num=rand(10000);
$tempo="/home/schles/norsnet/work/tempo$rand_num";
$tempo50="/home/schles/norsnet/work/$file";
system ("cat $tempo50 | wc > $tempo");
open (F, "$tempo") || die "siyut ahhh $tempo";
$line=<F>;
 if ($line=~/^\s*(\d+)\s*/) {
	$num2=$1;
	}
close (F);
system ("rm $tempo");
$sampIn=$num2-1;
$id=$file;
print "\n####working on $id\n";
$id=~ s/\.data//;
$oldid=$id;
$outFinal= "$resultsdir/$id.norsnet";
$rand2=$rand_num=int(rand(100000));
$id="tmp$rand2";
print "$sampIn samples\n";if ($sampIn==0)  {next loop2}	 
$inTest="$dir/$id-in_test";
$inOutTest="$dir/$id-in-out-test";
#$jctFile="jct$u-samp10states-$sampIn-Balanced-win9";
if ($id=~/help/)  {
	$idtemp=$id;
	$idtemp=~s/help//;
	$partest="$dir/$idtemp-partest";
	}
else {
	$partest="$dir/$id-partest";
	}
$testOutFile="$resultsdir/$id-50-5HN-416-NBba2rel13";
$testOutFile2="$resultsdir/$oldid-50-5HN-416-NBba2rel13";	
print "1\n";
open (FOUT, ">$inTest") || die "error0";
printf FOUT "* overall: (A,T25,I8)\nNUMIN                 :      %3d\nNUMSAMFILE            :   %6d\n*",$IN,$sampIn;
print FOUT "\n* samples: count (A8,I8) NEWLINE 1..NUMIN (25I6)\n";
print "#####collecting all samples#######\n";
$h=1;
undef @res;undef $end;undef@PREL;  undef @otL; undef @otE;
undef @otH;undef @RI_A;$expCon1=$expCon2=0;undef @RI_S;undef @outPut;
undef @secC;  $lengthA=$lengthB=$lengthC=0;
undef @A;undef @C;undef @D;undef @E;undef @F;undef @G;undef @H;undef @I;undef @K;undef @L;
undef @M;undef @N;undef @P;undef @Q;undef @R;undef @S;undef @T;undef @V;undef @W;undef @Y;undef @GS;
chomp($id);undef @RI_S;undef $hydroNet;$totalDiff=0;undef @diff1;undef @diffNew;undef @diff11;
;undef @net2;undef @net3;
#$net1=$id . ".profbval";	
if (!(open (FINO, "$profbval"))) {
	print "no profbval file? cant open $profbval $!";
	next loop2;
	}
for ($f=0; $f<43; $f++)  {
  		 <FINO>;
       		 }	
loop1221:while ($line=<FINO>)  {
       	if ($line=~ /^(.{8})(.{5})(.{4})/)  {
               	$uno=$1;$duwe=$2; $tre=$3;
		$uno=~ s/\s//g;$duwe=~ s/\s//g;$tre=~ s/\s//g;
		push (@net2,$duwe); push (@net3,$tre); 
		 $diff1=$duwe-$tre;push (@diff1,$diff1);
		 $diff11=$diff1+100;
		 use integer;
		 $diff11=$diff11/2;
		 no integer;  
		 push (@diff11,$diff11);
		$totalDiff=$totalDiff+$diff1;
		}
	}
close (FINO);
@diffNew=normDiff($totalDiff,\@diff1);
$cA=$cC=$cD=$cE=$cF=$cG=$cH=$cI=$cK=$cL=$cM=$cN=$cP=$cQ=$cR=$cS=$cT=$cV=$cW=$cY=$compo=0;
open (FILE, "$dir/$file") || die "$dir/$file $!";
<FILE>; ## the format of these Datafiles is a bit different
while ($line=<FILE>)  {
	@stuff=split(' ', $line);
	
	$A=$stuff[1];$C=$stuff[2];$D=$stuff[3];$E=$stuff[4];$F=$stuff[5];$G=$stuff[6];$H=$stuff[7];$I=$stuff[8];
	$K=$stuff[9];$L=$stuff[10];$M=$stuff[11];$N=$stuff[12];$P=$stuff[13];$Q=$stuff[14];$R=$stuff[15];$S=$stuff[16];
	$T=$stuff[17];$W=$stuff[18];$Y=$stuff[19];$V=$stuff[20];			
	$otH=$stuff[22];$otE=$stuff[23];$otL=$stuff[24];$PREL=$stuff[25];$RI_A=$stuff[26];
	$expCon1=$stuff[27];$expCon2=$stuff[28];$lengthA=$stuff[29];$lengthB=$stuff[30];$lengthC=$stuff[31];$outPut=$stuff[32];$res=$stuff[33];
       	push (@res,$res);
	push (@A,$A) ;push (@C,$C);push (@D,$D);push (@E,$E);push (@F,$F);push(@G,$G);push(@H,$H);push(@I,$I);push(@K,$K);push(@L,$L);
	push (@M,$M);push(@N,$N);push(@P,$P);push(@Q,$Q);push(@R,$R);push(@S,$S);push(@T,$T);push(@V,$V);push(@W,$W);push(@Y,$Y);			
      	push (@otH,$otH);push (@outPut,$outPut);
       	push (@otE,$otE);
      	push (@secC,$secC);
       push (@otL,$otL);
        push (@RI_A,$RI_A);
      	push (@Bnew,$Bnew);
       	push (@PREL,$PREL);;$end=scalar@PREL-1;
	if ($res eq 'A')  {
		$cA++;
		}
	elsif($res eq 'C')  {
		$cC++;
		}
	elsif($res eq 'D')  {
		$cD++;
		}			
	elsif($res eq 'E')  {
		$cE++;
		}			
	elsif($res eq 'F')  {
		$cF++;
			}			
	elsif($res eq 'G')  {
		$cG++;
		}			
	elsif($res eq 'H')  {
		$cH++;
		}
	elsif($res eq 'D')  {
		$cD++;
		}			
	elsif($res eq 'E')  {
		$cE++;
		}			
	elsif($res eq 'F')  {
		$cF++;
		}			
	elsif($res eq 'G')  {
		$cG++;
		}			
	elsif($res eq 'H')  {
		$cH++;
		}			
	elsif($res eq 'I')  {
		$cI++;
		}			
	elsif($res eq 'K')  {
		$cK++;
		}			
	elsif($res eq 'L')  {
		$cL++;
		}
	elsif($res eq 'M')  {
		$cM++;
		}			
	elsif($res eq 'N')  {
		$cN++;
		}			
	elsif($res eq 'P')  {
		$cP++;
		}			
	elsif($res eq 'Q')  {
		$cQ++;
		}						
	elsif($res eq 'R')  {
		$cR++;
		}
	elsif($res eq 'S')  {
		$cS++;
		}			
	elsif($res eq 'T')  {
		$cT++;
		}			
	elsif($res eq 'V')  {
		$cV++;
		}			
	elsif($res eq 'W')  {
		$cW++;
		}		
	elsif($res eq 'Y')  {
		$cY++;
		}		
	}	
$compo=scalar@res;
if ($compo==0)  {
	print "\n$file is too short, has only $compo residues \n";
	die;
	}
$compo=scalar@res;				
$cA=$cA/$compo*100;$cC=$cC/$compo*100;$cD=$cD/$compo*100;$cE=$cE/$compo*100;$cF=$cF/$compo*100;$cG=$cG/$compo*100;$cH=$cH/$compo*100;
$cI=$cI/$compo*100;$cK=$cK/$compo*100;$cL=$cL/$compo*100;$cM=$cM/$compo*100;$cN=$cN/$compo*100;$cP=$cP/$compo*100;$cQ=$cQ/$compo*100;
$cR=$cR/$compo*100;$cS=$cS/$compo*100;$cT=$cT/$compo*100;$cV=$cV/$compo*100;$cW=$cW/$compo*100;$cY=$cY/$compo*100;
use integer;
$cA=$cA*1;$cC=$cC*1;$cD=$cD*1;$cE=$cE*1;$cF=$cF*1;$cG=$cG*1;$cH=$cH*1;$cI=$cI*1;$cK=$cK*1;$cL=$cL*1;$cM=$cM*1;$cN=$cN*1;
$cP=$cP*1;$cQ=$cQ*1;$cR=$cR*1;$cS=$cS*1;$cT=$cT*1;$cV=$cV*1;$cY=$cY*1;$cY=$cY*1;
no integer;
$hydroNet=hydroNet($end);
if ($hydroNet>=8) {$hydroNet=0}
else {$hydroNet=100}		
loop4:for ($i=0;$i<scalar@PREL; $i++)  {		
	$k=0;undef @info;
	if ($h==($sampIn + 1))  {
		last loop4;
		}
	$lower=$i-($win-1)/2;
	$higher=$i+($win-1)/2;
#profiles information
#	push (@info,profiles($lower,$higher,$end));
#secondary structure prediction information
	push (@info, secondary($lower,$higher,$end));
#loop for solvent accessibility prediction information		
	push (@info, acc($lower,$higher,$end));
# global information
#	push (@info,$RI_A[$i] );
	push (@info, $expCon1,$expCon2);	
#	push (@info, $Helix, $Beta, $Loop);
	push (@info, $lengthA,$lengthB,$lengthC);
	push (@info,$cA,$cC,$cD,$cE,$cF,$cG,$cH,$cI,$cK,$cL,$cM,$cN,$cP,$cQ,$cR,$cS,$cT,$cV,$cW,$cY);
	push (@info,profiles($lower,$higher,$end));
	push (@info,$hydroNet);
	push (@info, getDiff($lower,$higher,$end));
	printf FOUT "ITSAM:%10d\n",$h;
	presentIt(\@info);
	$h++;
#	presentIt(\@info);		
#	print FOUT "\n" unless $k==25;
	}
print FOUT "//"; 
close(FOUT);
print "\nnumber of samples saved in memory:$h\n"; 
#Mark2:	
print "\n########creating out-test file ######\n";
open (FOUT, ">$inOutTest") || die "cant open file $!";
printf FOUT "* overall: (A,T25,I8)\nNUMOUT                :        $ON\nNUMSAMFILE            :%9d\n*",$sampIn;
print FOUT "\n* samples: count (I8) SPACE 1..NUMOUT (25I6)\n";
for ($i=0; $i<scalar@outPut;$i++)  {
	$ii=$i+1;
	if ($outPut[$i]==100)  {
		$GS[$i]='G';
		}
	else {
		$GS[$i]='-';
		}
	printf FOUT "%8d",$ii; $m=100- $outPut[$i];
      		printf FOUT "  %5d%6d\n",$outPut[$i],$m;
	}
print FOUT "//"; 
print "\t$i\n";
print 1;
close (FOUT);
open (FOUT, ">$partest")  || die "can't open file $!";
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
print FOUT "EPSILON               :        0.001000\n";
print FOUT "ALPHA                 :        0.010000\n";
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
print FOUT "FILEIN_IN             : $inTest\n";
print FOUT "FILEIN_OUT            : $inOutTest\n";
print FOUT "FILEIN_JCT            : $jct\n";
print FOUT "FILEOUT_OUT           : $testOutFile\n";
print FOUT "FILEOUT_JCT           : jct_crap\n";
print FOUT "FILEOUT_ERR           : NNo_tst_err.dat\n";
print FOUT "FILEOUT_YEAH          : NNo-yeah1637.tmp\n";
print FOUT "//\n";
system ("/home2/schles/palm-share/antigenicity/NET/for/NetRun2.LINUX $partest");
system ("rm $inTest");system ("rm $inOutTest");
system ("rm $partest");undef @nors40; undef @nors40f;undef @result1;undef @result2;undef @nors59; undef @nors59f;undef @pred;
print "cleaning up. removing : $inTest $inOutTest $partest\n";
open (F, "$testOutFile") ||  next;
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
		push (@pred, $w);
		if ($w<0.40) {
			push (@nors40,'-');
			}
		else {
			push (@nors40,'N');
			}
		if ($w<0.59) {
			push (@nors59,'-');
			}
		else {
			push (@nors59,'N');
			}
		}
	}
		
close (F);
print "mv $testOutFile $testOutFile2\n";
system ("mv $testOutFile $testOutFile2");
$v=$sofer=0;
for ($r=0;$r<scalar@pred;$r++) {
$nors40f[$r]='-';
	if ($nors40[$r] eq '-') {
		$sofer=0;$v=0;
		next;
 		}
	else {
 		$sofer++;
		if ($sofer>30){ 
			if ($v==0)  {
				$v++;
				#$nors=1;
				$start=$r-30;
				$end=$r;
				for ($u=$start;$u<=$end;$u++) {
					$nors40f[$u]='N';
					}
				}
			else {
				$nors40f[$r]='N';
				}
			}			
                 }
	}
$v=$sofer=0;
for ($r=0;$r<scalar@pred;$r++) {
	$nors59f[$r]='-';
	if ($nors59[$r] eq '-') {
		$sofer=0;$v=0;
		next;
 		}
	else {
 		$sofer++;
		if ($sofer>30){ 
			if ($v==0)  {
				$v++;
				#$nors=1;
				$start=$r-30;
				$end=$r;
				for ($u=$start;$u<=$end;$u++) {
					$nors59f[$u]='N';
					}
				}
			else {
				$nors59f[$r]='N';
				}
			}			
                 }
	}
open (FOUT2, ">$outFinal") || die "cant open $outFinal $!";
print FOUT2 "pos\tres\tnode1\tnode2\tpred\tn40\tn40fil\tn59\tn59fil\n";
for ($r=0;$r<scalar@pred;$r++) {
	$r1=$r+1;
	printf FOUT2 "$r1\t$res[$r]\t$result1[$r]\t$result2[$r]\t%1.2f\t$nors40[$r]\t$nors40f[$r]\t$nors59[$r]\t$nors59f[$r]\n",$pred[$r];
	}
system ("rm /home/schles/norsnet/work/$file");
	
#function  tha shuffles the array randomly
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
#functions to obtain properties    
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
			@secon=(100,100,100);
			push (@array, @secon);
			}
		else {
			undef @secon;
			@secon=($otL[$j],$otH[$j],$otE[$j]);
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
			@PRE=(100,100);
			push (@array, @PRE);
			}
		else {
			undef @PRE;
			@PRE=($PREL[$j],$RI_A[$j]);
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
sub gaa {
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
			if ($res[$j] eq 'C') {
				@residue=(100,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
				}
			elsif ($res[$j] eq 'F') {
				@residue=(0,100,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
				}
			elsif ($res[$j] eq 'M') {
				@residue=(0,0,100,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
				}
			elsif ($res[$j] eq 'I') {
				@residue=(0,0,0,100,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
				}
			elsif ($res[$j] eq 'G') {
				@residue=(0,0,0,0,100,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
				}
			elsif ($res[$j] eq 'P') {
				@residue=(0,0,0,0,0,100,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
				}				
			elsif ($res[$j] eq 'A') {
				@residue=(0,0,0,0,0,0,100,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
				}				
			elsif  ($res[$j] eq 'H'){
				@residue=(0,0,0,0,0,0,0,100,0,0,0,0,0,0,0,0,0,0,0,0,0);
				}
			elsif ($res[$j] eq 'Q') {
				@residue=(0,0,0,0,0,0,0,0,100,0,0,0,0,0,0,0,0,0,0,0,0);
				}
			elsif ($res[$j] eq 'R') {
				@residue=(0,0,0,0,0,0,0,0,0,100,0,0,0,0,0,0,0,0,0,0,0);
				}
			elsif ($res[$j] eq 'W') {
				@residue=(0,0,0,0,0,0,0,0,0,0,100,0,0,0,0,0,0,0,0,0,0);
				}
			elsif ($res[$j] eq 'L') {
				@residue=(0,0,0,0,0,0,0,0,0,0,0,100,0,0,0,0,0,0,0,0,0);
				}
			elsif ($res[$j] eq 'Y') {
				@residue=(0,0,0,0,0,0,0,0,0,0,0,0,100,0,0,0,0,0,0,0,0);
				}
			elsif ($res[$j] eq 'T') {
				@residue=(0,0,0,0,0,0,0,0,0,0,0,0,0,100,0,0,0,0,0,0,0);
				}
			elsif ($res[$j] eq 'K') {
				@residue=(0,0,0,0,0,0,0,0,0,0,0,0,0,0,100,0,0,0,0,0,0);
				}
			elsif ($res[$j] eq 'S') {
				@residue=(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,100,0,0,0,0,0);
				}
			elsif ($res[$j] eq 'N') {
				@residue=(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,100,0,0,0,0);
				}
			elsif ($res[$j] eq 'V') {
				@residue=(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,100,0,0,0);
				}
				
			elsif ($res[$j] eq 'D') {
				@residue=(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,100,0,0);
				}
				
			elsif ($res[$j] eq 'E') {
				@residue=(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,100,0);
				}
			else {
				@residue=(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
				}				
			push (@array, @residue);
			}
		}
		return @array;
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
