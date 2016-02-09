#!/usr/bin/perl -w

#------------------------------------------------------------------------------#
#	Copyright				        	1998	       #
#	Burkhard Rost		rost@EMBL-Heidelberg.DE			       #
#	Wilckensstr. 15		http://www.embl-heidelberg.de/~rost/	       #
#	D-69120 Heidelberg						       #
#				version 0.1   	May,    	1998	       #
#------------------------------------------------------------------------------#
				# sets array count to start at 1, not at 0

use RG::Utils::Hssp_filter;

# lkajan: trying to move away from $pack
#$pack="hssp_filterPack";
($Lok,$msg)=
    &RG::Utils::Hssp_filter::hssp_filterSbr( 'dirBin=/usr/bin', 'dirHome=/usr/share/librg-utils-perl', @ARGV);

die( "*** RG::Utils::Hssp_filter::hssp_filterSbr returned ERROR:\n".$msg." " ) if (! $Lok);

exit(0);

__END__




