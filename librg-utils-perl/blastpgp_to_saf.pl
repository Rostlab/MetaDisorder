#!/usr/bin/perl -w
use Carp qw| cluck :DEFAULT |;
use Data::Dumper qw||;
use List::MoreUtils;
our $debug;

if ($#ARGV < 0) { print  "*** ERROR blastpgp_to_saf.pl : no arguments given\n"; print  FHERROR "*** ERROR blastpgp_to_saf.pl : no arguments given\n"; die; }

foreach $arg (@ARGV){
        if    ($arg=~/^fileInBlast=(.*)$/)       { $fileInBlast =        $1;}
        elsif ($arg=~/^fileInQuery=(.*)$/)       { $fileInQuery =        $1;}
	elsif ($arg=~/^fileOutSaf=(.*)$/)        { $fileOutSaf  =        $1;}
	elsif ($arg=~/^fileOutRdb=(.*)$/)        { $fileOutRdb  =        $1;}
	elsif ($arg=~/^fileOutErr=(.*)$/)        { $fileOutErr  =        $1;}
	elsif ($arg=~/^red=(.*)$/)               { $filterThre  =        $1;}
	elsif ($arg=~/^maxAli=(.*)$/)            { $maxAli      =        $1;}
        elsif ($arg=~/^tile=(.*)$/)              { $alignTiling =        $1;}
        elsif ($arg=~/^debug=(.*)$/)             { $debug =              $1;}
	else {
	    print FHERROR "*** wrong command line arg '$arg'\n";
	    die;
	}
}

if ($#ARGV < 0) { print  "*** ERROR blastpgp_to_saf.pl : no arguments given\n"; 
		  print  FHERROR "*** ERROR blastpgp_to_saf.pl : no arguments given\n"; die; }

#if (! defined $fileOutErr)  { 
#    print "*** ERROR blastpgp_to_saf : fileOutErr output filename is not definded\n";
#    die;}

if( defined( $fileOutErr ) )
{
	open(FHERROR,">", $fileOutErr) || die "*** ERROR could not open error log file=$fileOutErr for writing";
}
else
{
	open( FHERROR, '>&', \*STDERR ) || die( "$!" );
}
$inicheck=0;
if (! defined $fileInBlast)  { 
    print FHERROR "*** ERROR blastpgp_to_saf : blast input file name is not defined\n"; 
    $inicheck++;}
else { if (! -e $fileInBlast){ 
    print FHERROR "*** ERROR blastpgp_to_saf : no input blast file $fileInBlast found\n";
    $inicheck++;} }

if (! defined $fileInQuery)  { 
    print FHERROR "*** ERROR blastpgp_to_saf : query input file name is not defined\n";
    $inicheck++;}
else{ if(! -e $fileInQuery)  { 
    print FHERROR "*** ERROR blastpgp_to_saf : no input query file $fileInQuery found\n";
    $inicheck++;} }

if (! defined $fileOutSaf)  { 
    print FHERROR "*** ERROR blastpgp_to_saf : SAF output filename is not definded\n";
    $inicheck++;}
if (! defined $fileOutRdb)  { 
    print FHERROR "*** ERROR blastpgp_to_saf : blastRdb output filename is not definded\n";
    $inicheck++;}

if (! defined $filterThre)  { 
    print FHERROR "*** ERROR blastpgp_to_saf : filter threshold is not definded\n";
    $inicheck++;}
if (! defined $maxAli)      { 
    print FHERROR "*** ERROR blastpgp_to_saf : maximum number of aligned sequences to be included in saf output file is not definded\n";
    $inicheck++;}
if (! defined $alignTiling) { 
    print FHERROR "*** ERROR blastpgp_to_saf : tiling method of Blast alignment  is not definded\n";
    $inicheck++;}

die           if ($inicheck != 0);

($Lok,$msg)=   &blastp_to_saf($fileInBlast,$fileInQuery,$fileOutSaf,$fileOutRdb,$filterThre,$maxAli,$alignTiling);

if (! $Lok )    { print FHERROR "*** ERROR blastpgp_to_saf.pl : $msg\n"; die;}
if ($Lok == 2 ) { if( $debug ){ print FHERROR "*** WARNING blastpgp_to_saf.pl : $msg\n"; } exit(0); }

exit(0);


#=============================================================================================
sub blastp_to_saf {
    my ($Lok,$msg);
    my $sbr="blastp_to_saf";	# 
    local ($blastfile,$queryfile,$fileout,$rdb,$filter,$maxAli,$tile)=@_; # 
    local (@query, %sequences,@alignedids,@namesSort,%rdb_lines); # 
    local (@tmp_seq,@inserted_query,@seq,@alignedNames ); # 
    local ($key,$first,$last,$fhin,$local_counter,$beg,$index,$iter,$Score_count); 

    $fhin= "FHIN"; 		# 
                                       #------------------- gets the query sequence and its length
    undef @query;		# 
    open($fhin,$queryfile) || return(0,"*** ERROR sbr: could not open $queryfile  - no such file");
    $queryName='query';		# 
    while(<$fhin>){
	next   if( $_=~/^\n/ || $_=~/\// || $_=~/>/ || $_=~/#/ );
	$_=~s/\s+//g;
	@tmp=split(//,$_);
	push @query, @tmp;
    }
    close $fhin;
    $queryLength=$#query+1;
                                #..........................................................
    
    if ( $rdb ne '0' && $rdb ne '' ){ 
	$rdbFile=$rdb;
	($Lok,$msg)= &printRdbHeader();
	if (! $Lok){ return(0,"*** ERROR $sbr : $msg"); }
    }

                                #------------------- finds number of iterations in blast file
    my $blastplus_flag = 0;

    open($fhin,$blastfile) || return(0,"*** ERROR $sbr: failed to open blast file $blastfile\n");
    $iter=0; $nohits=0;
    while(<$fhin>){
	# BLASTP 2.2.18 [Mar-02-2008]
	# PSIBLAST 2.2.25+
	if($_=~/^\w+\s+\d+\.\d+\.\d+\+/o){ $blastplus_flag = 1; }
	if($_=~/No hits found/){$nohits=1; last;}
	if( !$blastplus_flag && $_=~/^Searching../o ){       $iter++; }
	if(  $blastplus_flag && $_=~/Results from round/o ){ $iter++; }
    }
    close $fhin;

	

    if ($nohits eq '1'){
	$sequences{''}=''; $alignedNames='';
	($Lok,$msg)=&print_saf_file($queryName,$queryLength,\%sequences,\@alignedNames,@query);
	return(0,"*** ERROR $sbr : $msg")            if (! $Lok );
	return(2,"*** WARNING $sbr : no hits found in Blastpgp search ");  
    }
    return(0,"*** ERROR $sbr: blast file format not recognized")              if ($iter == 0);
                                #............................................................

                                #--------------------------------skips to the last iteration
    open($fhin,$blastfile) || return(0,"*** ERROR $sbr: failed to open blast file $blastfile\n");
    $local_counter=0;
    while(<$fhin>){
	if( !$blastplus_flag && $_=~/^Searching../o ){       $local_counter++; }
	if(  $blastplus_flag && $_=~/Results from round/o ){ $local_counter++; }
	last if($local_counter == $iter);
    }
   
    #............................................................

	my $bug_query = 0;
    undef @alignedNames; undef @alignedids; $Score_count=0; undef %multi_aligned; $global_count=0; undef %rdb_lines;
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    # lkajan: processing of a record is done when the beginning of the /following/ record is reached.
    # lkajan: the last record therefore is processed when we reach the end of the input file - we do not rely on the
    # lkajan: Number of Hits line any more - this allows truncated blast input (some user in Yanay's lab has such)
    while(1)
    {
	my $bline = <$fhin>;
	if( !$bline || $bline=~/^>/o )
	{
	    if($global_count > 0){
		undef %u_endings;
		@seq=@{ $multi_aligned{'1'} };
		if($tile ne "0"  && $Score_count > 1){
		    $u_endings{'1'}[0]=$endings{'1'}[0]; $u_endings{'1'}[1]=$endings{'1'}[1];
		    foreach $itnum (2 .. $Score_count){
			@temp=@{ $multi_aligned{$itnum} };
			$iffy=1;
			foreach $it ( 1 .. ($itnum-1)){	
			    if (defined $u_endings{$it}){
				if ($endings{$itnum}[0] >= $u_endings{$it}[0]  && $endings{$itnum}[0] <= $u_endings{$it}[1] ){ $iffy=0;last;}
				if ($endings{$itnum}[1] >= $u_endings{$it}[0]  && $endings{$itnum}[1] <= $u_endings{$it}[1] ){ $iffy=0;last;}
			    }
			}
			if ($iffy==1){
			    $u_endings{$itnum}[0]=$endings{$itnum}[0]; $u_endings{$itnum}[1]=$endings{$itnum}[1];
			    $one=$u_endings{$itnum}[0]; $two=$u_endings{$itnum}[1];
			    @seq[$one .. $two]=@temp[$one .. $two];    
			}	
			else{delete $rdb_lines{$id}{$itnum}; }
		    }
		}
		foreach $elem (@seq){
		    if(($elem ne ".") && ($elem !~ /[a-z_A-Z]/)){
			$elem=".";
		    }
		}
		$sequences{$id}=[ @seq ];
	    }
	    $global_count++;

	    undef @seq; $Score_count=0; undef %multi_aligned; undef %endings; 
	    undef %u_endings;
							                   #--- getting name of aligned sequence
	    for($it=0;$it<=$queryLength-1;$it++){   $seq[$it]="."; }       #initialising array seq
	
	    if( $bline )
	    {
# lkajan: PE, SV may be on the > line when it is shorter
# blast:
#>sp|Q9TUI5|MT4_CANFA Metallothionein-4 OS=Canis familiaris GN=MT4
#          PE=3 SV=1
#          Length = 62
# blast+:
#> sp|Q9TUI5|MT4_CANFA Metallothionein-4 OS=Canis familiaris GN=MT4 
#PE=3 SV=1
#Length=62
	    	if( not $bline=~s/^> *//o ){ confess("unrecognized line '$bline'"); }
	    	$id=$bline; if( not $id=~s/^(\S*)\s+(.*)\s*$/$1/o ){ confess("failed to extract id from '$id'"); }
	    	$protDspt=$2;  chomp $protDspt; 
	    	push @alignedids, $id;
	    }
	    else { last; }
	    
	}                                                   #------------------------------------
	if ($bline=~/^ Score/o){ 
	       $Score_count++; 
	       for($it=0;$it<=$queryLength-1;$it++){ $block_seq[$it]="."; }
	       undef @ali_para;
	}

	next                   if ( $Score_count > 1 && $tile==0 );
#-----------------------------------------------------------------------------------------------------
	if ( $rdb ne '0' && $rdb ne ''){
# 
	    if ( $bline=~ /\s*Length/){ $len2=$bline;if( $len2 !~ /Length\s*=\s*([0-9]+)/o ){ confess( "unrecognized length '$len2'"); } $len2 = $1; }
# lkajan: whitespaces vary between blast and blast+
# Score = 84.3 bits (207),  Expect = 2e-15, Method: Compositional matrix adjust.
	    if ( $bline=~ /Score/ ){
		$lali=$pid=$sim=$bitScore=$expect=''; $gap=0;
		$line=$bline;
		chomp $line; @tmp=split(/,/o,$line);
		push @ali_para,@tmp;
	    }
# 2nd: blast+
# Identities = 57/62 (91%), Positives = 59/62 (95%)
# Identities = 57/62 (91%), Positives = 59/62 (95%), Gaps = 0/62 (0%)
	    if ( $bline=~ /Identities/){
		$line=$bline; chomp $line;
		@tmp=split(/,/o,$line); push @ali_para,@tmp;
		foreach $param (@ali_para){
		    $param=~s/\s+//g;
		    if ($param=~/Score/){ $bitScore=$param; if( not $bitScore=~s/^.+=(.+)bits.*$/$1/o ){ confess($bitScore); } }
		    elsif ($param=~/Expect/){@temp=split(/=/,$param); $expect=$temp[1];}
		    elsif ($param=~/Identities/){$lali=$param; if( not $lali=~s/.+\/([0-9]+)\(([0-9]+)%\).*$/$1/o ){ confess( $lali ); }
						 $pid=$2;}
		    elsif ($param=~/Positives/){$sim=$param; if( not $sim=~s/^.+\(([0-9]+)%\).*$/$1/o ){ confess( $sim ); }}
		    elsif ($param=~/Gaps/){$gap=$param; if( not $gap=~s/^.+=([0-9]+)\/.*$/$1/o ){ confess($gap); }}
		}
		$lali=$lali-$gap;
		$qLength=$queryLength;
		$rdb_lines{$id}{$Score_count}=[$qLength,$len2,$lali,$pid,$sim,$gap,$bitScore,$expect,$protDspt];
	    } 	
	}
#-----------------------------------------------------------------------------------------------------
	# lkajan: below is probably a blast bug, with the Query: 0 line...
	# 0      1   2
	# Query: 0   ----                                                        
	#
	# Sbjct: 64  PAQG                                                         67
	# Query: 0                                                               
	#
	# Sbjct: 67                                                               67
	# Query: 1713 HVETRWHCTVCEDYDLCINCYNTKSHAHKMVKWGLGLDDEGSSQGEPQSKSPQESRRVSI 1772
	# Query: 1890 PGTPTQQPSTPQT 1902
	# blast
	#Query: 1  MDPRECVCMSGGICMCGDNCKCTTCNCKTCRKSCCPCCPPGCAKCARGCICKGGSDKCSC 60
	#          MDP EC CMSGGIC+CGDNCKCTTCNCKTCRKSCCPCCPPGCAKCA+GCICKGGSDKCSC
	#Sbjct: 1  MDPGECTCMSGGICICGDNCKCTTCNCKTCRKSCCPCCPPGCAKCAQGCICKGGSDKCSC 60
	# blast+
	#Query  1   MDPRECVCMSGGICMCGDNCKCTTCNCKTCRKSCCPCCPPGCAKCARGCICKGGSDKCSC  60
	#           MDP EC CMSGGIC+CGDNCKCTTCNCKTCRKSCCPCCPPGCAKCA+GCICKGGSDKCSC
	#Sbjct  1   MDPGECTCMSGGICICGDNCKCTTCNCKTCRKSCCPCCPPGCAKCAQGCICKGGSDKCSC  60

	# lkajan: do not confuse with 'Query=  MT4_HUMAN' - this may appear in blast+ after each round header
	if($bline=~/^Query:?\s+\d+/o)
	{
		@tmp=split(/\s+/o,$bline);
		if( $tmp[1] != 0 )
		{
			$bug_query = 0;
			undef @aligned; undef @inserted_query;
			$beg=$tmp[1]-1; $end=$tmp[3]-1;
			if (! defined $endings{$Score_count}[0] || $endings{$Score_count}[0] < 0 ){ $endings{$Score_count}[0]=$beg;}
			$endings{$Score_count}[1]=$end;
			@inserted_query=split(//o,$tmp[2]);
		}
		else
		{
			$bug_query = 1;
		}
	}
	if( !$bug_query && $bline =~ /^Sbjct/o ){
	    @tmp=split(/\s+/o,$bline); 
	    @aligned=split(//o,$tmp[2]);
						     #getting rid of insertions at query sequence
	    cluck( " *** ERROR sbr: blastp_to_saf in lenghts for $id" )  if ($#inserted_query != $#aligned);
	    $local_counter=0;
	    undef @tmp_seq;
	    for($it=0;$it <= $#inserted_query; $it++){
		if ($inserted_query[$it] =~ /[a-z_A-Z]/){
		    $tmp_seq[$local_counter]=$aligned[$it];
		    $local_counter++;
		}
	    }
	    #@aligned=@tmp_seq;                     #-------------------------------------------
#+++++++++++=
	    
	    @block_seq[$beg .. ($beg+$#tmp_seq)]=@tmp_seq;    #----alignig part of the subject seguence
	    
	    $multi_aligned{$Score_count}=[ @block_seq ]; 
	}
    }
    close $fhin;

    @alignedids = sort{ $a cmp $b }@alignedids;

				                      #getting rid of repeats in the list    
    undef @namesSort;
    push @namesSort, $queryName;
#    $Lname=$queryName; $Lname=~tr/A-Z/a-z/;
#    $Cname=$Lname; $Cname=~tr/a-z/A-Z/;
#    foreach $it (@alignedids){
#	$rflag=0; @tmp=split(/\|/,$it);
#	if( $it eq $queryName || $it eq $Lname || $it eq $Cname){ $rflag=1;}
#	else { 
#	    foreach $elem (@tmp){
#		if($elem eq $queryName || $elem eq $Lname || $elem eq $Cname){$rflag=1;}
#	    }
#	}
#	if($rflag == 0){push @namesSort, $it;}
#    }

    push @namesSort, @alignedids;
    @namesSort = List::MoreUtils::uniq( @namesSort );
    @alignedNames = splice( @namesSort, 1 );

    if( $debug ){ warn( Data::Dumper::Dumper( \@alignedNames ) ); }
                                                       #-------------------------------------
 
                                                       #--filtering alignment
    
    if($filter != 100){
	($Lok,$msg)=&saf_filter(\@alignedNames,\%sequences,$filter,$maxAli);
	return(0,"*** ERROR $sbr : $msg")               if(! $Lok );
    }
  
                          #----------------------------------- printing out the resulting files
    
    if ($rdb ne '0' ){
	foreach $id (@alignedNames){
	    $identifier=$id;
	    foreach $score (sort keys %{ $rdb_lines{$id} }){		
		#if($score > 1){ foreach $it ( @{ $rdb_lines{$id}{$score} }){$it='!'.$it;} };
		($qLength,$len2,$lali,$pid,$sim,$gap,$bitScore,$expect,$protDspt)= @{ $rdb_lines{$id}{$score} };
		if($score > 1){ $expect='!'.$expect;
				$protDspt='!'.$protDspt;
				$identifier='!'.$identifier;
		}
		print FHRDB "$identifier\t$len2\t$pid\t$sim\t$lali\t$gap\t$bitScore\t$expect\t$protDspt\n";  
	    }
	}
    } 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      short naming format
#    undef %short_alignedNames;
#    foreach $key (@alignedNames){
#	if($key =~ /trembl/ || $key =~/swiss/){
#	    @tmp=split(/\|/,$key); $short=$tmp[2];
#	}
#	else { $short=$key;}
#	$short_alignedNames{$key}=$short;
#    }
    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ($Lok,$msg)=&print_saf_file($queryName,$queryLength,\%sequences,\@alignedNames,@query);
    return (0,"*** ERROR $sbr : $msg")                        if ( ! $Lok );
    return (2,"*** WARNIG $sbr : no hits found after filtering")             if ( $#alignedNames < 0);
    if (-e $fileout) { return (1,"ok"); }
    else { return (0,"*** ERROR $sbr : failed to produce the saf output file"); }
}				# end of blastp_to_saf

#==============================================================================

sub saf_filter{
    my $sbr='saf_filter';
    my ($Lok,$msg);
    local($array_name,$ali_hash,$red,$bound)=@_;
    my @ord_ali_list=@$array_name;
    my %alignment=%$ali_hash;
    my @new_aligned_names=();
    my ($count,$ct,$same);

    @last= @query;
    $len=$#last;
    $count=0;
    Floop:for($index=0; $index <= $#ord_ali_list; $index++){
	if ( ! defined $alignment{$ord_ali_list[$index]}) { return(0,"*** ERROR $sbr : alignment to be filtered is not defined\n")}
	@maybe=@{ $alignment{$ord_ali_list[$index]} };
	$ct=$same=0;
	foreach $itres (0..$len){
	    next if ($maybe[$itres] !~ /[a-zA-Z]/);
	    next if ( $last[$itres] !~ /[a-zA-Z]/);
	    ++$same if ($maybe[$itres] eq $last[$itres]);
	    ++ $ct;
	}
	if ( ( $ct > 0 && (100*$same/$ct)<$red ) || $ct==0 ){
	    @last=@maybe; $count++;
	    push @new_aligned_names, $ord_ali_list[$index];
	    last Floop    if ($count >= $bound); 
	}
    }
    @$array_name=@new_aligned_names;
    return(1,"filtering is done");
}
	
#==============================================================================================
#==============================================================================================

sub printRdbHeader{
    my $sbr='printRdbHeader';
    my ($Lok,$msg);

    $header="
#PERL-RDB
#SEQLENGTH\t $queryLength
#ID\t:\tidentifier of the aligned (homologous) protein
#LSEQ2\t:\tlength of the entire sequence of the aligned protein
#LALI\t:\tlength of the alignment excluding insertions and deletions
#%IDE\t:\tpercent indentity
#%SIM\t:\tpercent similarity
#LGAP\t:\ttotal gap length
#BSCORE\t:\tblast score (bits)
#BEXPECT\t:\tblast expectation value
#PROTEIN\t:\tone-line description of aligned protein
#'!'\t:\tindicates lower scoring alignment that is combined with 
#the higher scoring adjacent one  
##ID\tLSEQ2\t%IDE\t%SIM\tLALI\tLGAP\tBSCORE\tBEXPECT\tPROTEIN\n";

    open(FHRDB,">$rdbFile") || return(0, "*** ERROR $sbr :  could not open rdbfile=$rdbFile for writing\n");
    print FHRDB $header;
    return( 1,'ok');
}
#==============================================================================================
#==============================================================================================
sub print_saf_file{
    my $sbr='print_saf_file';
    my ($Lok,$msg);
    local($queryName,$queryLength,$alignment,$aliNames,@query)=@_;
    my %sequences=%$alignment;
    my @alignedNames=@$aliNames;
    my ($fhout,$pages,$nameField);

    $fhout="FHOUT"; 
    $pages=int $queryLength/50;
    if ($queryLength%50 != 0){$pages++;}

    open($fhout,">", $fileout ) || return(0,"*** ERROR $sbr: failed to open fileout=$fileout for writing");
    print $fhout "# SAF (Simple Alignment Format)\n";
    print $fhout "#\n";
    $nameField=0;
    @tmp=split(//,$queryName);
    if ($#tmp+2>$nameField){$nameField=$#tmp+2;}
    foreach $key (@alignedNames){
	@tmp=split(//,$key);
	if ($#tmp+2>$nameField){$nameField=$#tmp+2;}
    }


    for($it=1;$it<=$pages;$it++){
	$beg=($it-1)*50; $end=$it*50-1;
	printf $fhout "%-${nameField}.${nameField}s", $queryName;
	for($index=0;$index<50;$index=$index+10){
	    $first=$beg+$index; $last=$first+9;
	    if ($last <= $queryLength -1 ){print $fhout  @query[$first .. $last]," " ;}
	    else { print $fhout @query[$first .. ($queryLength-1)]; }
	}
	print $fhout  "\n";
	foreach $key (@alignedNames){
	    printf $fhout "%-${nameField}.${nameField}s", $key;
	    for($index=0;$index<50;$index=$index+10){
	       $first=$beg+$index; $last=$first+9; 
	       if ($last <= $queryLength -1 ){
		   print $fhout @{ $sequences{$key}}[$first .. $last]," " ;
	       }
	       else { print $fhout @{ $sequences{$key}}[$first .. ($queryLength-1)]," " ; }
	    }
	    print $fhout  "\n";
	}
	print $fhout "\n";
    }
    close $fhout;
    return (1,'ok');
}
#======================================================================





# vim:ai:ts=4:
