#!/usr/bin/perl -w
no warnings 'deprecated';

#------------------------------------------------------------------------------#
#	Copyright				        	1998	       #
#	Burkhard Rost		rost@EMBL-Heidelberg.DE			       #
#	Wilckensstr. 15		http://www.embl-heidelberg.de/~rost/	       #
#	D-69120 Heidelberg						       #
#				version 0.1   	May,    	1998	       #
#------------------------------------------------------------------------------#
$[ =1 ;				# sets array count to start at 1, not at 0

use RG::Utils::Copf;

#$pack=   "pack/copf.pm";
# lkajan: packName does not seem to be used in copf.pm
($Lok,$msg)=
    &RG::Utils::Copf::copf( 'dirBin=/usr/bin', 'dirHome=/usr/share/librg-utils-perl', @ARGV);

die( "*** RG::Utils::Copf::copf returned ERROR:\n".$msg." " ) if (! $Lok);

exit(0);
