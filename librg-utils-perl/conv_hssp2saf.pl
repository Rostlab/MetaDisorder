#!/usr/bin/perl
no warnings 'deprecated';
#------------------------------------------------------------------------------#
#	Copyright				        	1998	       #
#	Burkhard Rost		rost@EMBL-Heidelberg.DE			       #
#	Wilckensstr. 15		http://www.embl-heidelberg.de/~rost/	       #
#	D-69120 Heidelberg						       #
#				version 0.1   	May,    	1998	       #
#------------------------------------------------------------------------------#
$[ =1 ;				# sets array count to start at 1, not at 0

use RG::Utils::Conv_hssp2saf;

($Lok,$msg)=
    &RG::Utils::Conv_hssp2saf::conv_hssp2saf(@ARGV);

die( "*** package RG::Utils::Conv_hssp2saf returned ERROR:\n".$msg." " ) if (! $Lok);

exit(0);

