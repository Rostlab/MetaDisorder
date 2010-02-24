#!/usr/bin/perl
$dir=".";
opendir  (DIR, "$dir/hssp-fil") || die "cant open file $!";
@hssp_un= (grep /\-fil.hssp/, readdir(DIR));
@hssp=shuffle(\@hssp_un);
closedir(DIR);
foreach $hssp (@hssp) {
	$prof=$hssp;
	$prof=~ s/\.hssp/\.rdbProf/;

	$profbval=$hssp;$profbval=~s/-fil\.hssp/-63-R1-209-237433-ON2/;
	$norsnet=$profbval;
#	$norsnet=~s/-63-R1-209-237433-ON2/-50-5HN-416-NBba2rel13/;
	$norsnet=~s/-63-R1-209-237433-ON2/\.norsnet/;
	$eef="$dir/norsnet/".$norsnet;
    	if ( -e $eef ) {
        print "### Skipping file $eef since it exists\n";
        next;   }
	system (" perl /home/schles/norsnet/runNORSnetAlign.pl $dir/hssp-fil/$hssp $dir/prof/$prof $dir/profbval/$profbval $dir/norsnet");	
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
