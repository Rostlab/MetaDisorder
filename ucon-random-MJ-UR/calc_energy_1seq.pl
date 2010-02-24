#!/usr/local/bin/perl -w
$separation=$ARGV[0]; $id=$ARGV[1]; $work_dir=$ARGV[2]; $root=$ARGV[3];
$file="$work_dir/$id.f";
$fout="$work_dir/$id.eprofcon";
$dir="$root/ucon-random-MJ-UR";
get_values();
$seq='';undef %E;undef @seq;
open (FILE, "$file") || die "cant open file $!";
<FILE>;
while ($line=<FILE>) {
	$line=~ s/\n//;
	$seq= $seq . $line;
	$seq=~s/-//g; $seq=~tr/a-z/A-Z/;$seq=~s/ //g;
	@seq= split ('',$seq);
	}
close (FILE);
@aa=('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y');
loop5:for ($t=0;$t<scalar@seq;$t++) {
	foreach $aa (@aa) {
		if ($seq[$t] eq $aa) {
			next loop5;
			}
		}
	$seq[$t]='X';
	}
for ($n=0;$n<scalar@seq;$n++) {
	$aa3=$seq[$n];$aa1=$n;
	$start=$aa1-$separation;$end= $aa1+$separation;
	for ($m=$start;$m<=$end;$m++) {
		if ($m==$n) {
			next;
			}
		if (($m<0) || ($m>$#seq)) {
			next;
			}
		$aa4=$seq[$m];#$aa2=$m;
		if (($aa3 eq 'X') || ($aa4 eq 'X')) {
			$tempE=0;
			}
		else {
			$tempE= $value{$aa3}{$aa4};					
			push (@{$E{$aa1}},$tempE);
			#push (@{$E{$aa2}},$tempE); i use 2 loops so i dont need it here
			}
		}
	}
######################################## end of the prof section
open (FOUT,">$fout") || die "$fout";
for ($i=0;$i<scalar@seq;$i++) {
	 $aa=$i+1;$total=0;
         foreach  $e (@{$E{$aa}}) {
	 	if ((defined$total==0)|| (defined$e==0)) {
			print "$id total=$total e=$e\n";
			}
         	$total=$total+$e;
         	}
        if (scalar@seq<=$separation) {
                $totRes= scalar@seq;
                }
        elsif ((scalar@seq<=$separation*2)&&(scalar@seq>$separation)){
             if (($aa+$separation)>$#seq) {
              		 if (($aa-$separation)<=0)  {
                                  $totRes=scalar@seq
                                  }

                         elsif (($aa-$separation)>0) {
                                 $totRes=$separation+scalar@seq-$aa+1;
                                 }
                        }
             else {
                        $totRes=$aa+ $separation;
                       }
             }
       else {
             if (($aa-$separation)<=0) {
                        $totRes=$separation+$aa;
                        }
             elsif (($aa+$separation)>$#seq) {
                        $totRes=$separation+scalar@seq-$aa+1;
                        }
             else {
                       $totRes=$separation*2+1;
                       }
        }
	$energy=$total/$totRes;
	if ($totRes==($separation*2+2)) {
		print 1;
		}
	if ((defined$seq[$i]==0) ||(defined$energy==0) ||(defined$totRes==0)) {
		print "file=$file aa=$aa\tseq=$seq[$i]\tenergy=$energy\ttotRes=$totRes\n";
		}
	print FOUT "$aa\t$seq[$i]\t$energy\t$totRes\n";
		#print "$aa\t$seq[$i]\t$energy\t$totRes\t$secC[$i]\t$otH[$i]\t$otE[$i]\t$otL[$i]\n";
	}
close (FOUT);





sub get_values {
        open (F, "$dir/mj-upper-right-matrix.txt") || die "can't open file $!";
        $line=<F>;chomp$line;#$line=~ s/\s//g;
        @param=split (/\t/,$line);
        while ($line=<F>) {
                chomp$line;#$line=~ s/\s//g;
                undef @number;
                @number=split (/\t/,$line);
                $l=$number[0];push (@l,$l);
                for ($z=1;$z<scalar@number;$z++) {
			if ($number[$z]!~/\d/) {
				next;
				}
                        $value{$l}{$param[$z]}=$value{$param[$z]}{$l}=$number[$z];
                        }
                }
	$value{"X"}{"X"}=0;
        close (F);
        }
