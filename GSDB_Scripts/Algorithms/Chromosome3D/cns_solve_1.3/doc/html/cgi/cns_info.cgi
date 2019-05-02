#!/usr/bin/perl
# <!--
$this_cgi = 'cns_info.cgi';
$modified = 'Last modified 5/11/98 J. Reklaw';
#
# CGI file for displaying on-line help information
#
# Written by: Paul Adams and Jesse Reklaw
#
# copyright Yale University
#

require 'cns.lib';

##### PROGRAM FLOW OF CONTROL #####

$argstr=&init_cgi();
%input = &parse_args ($argstr);
$variable = $input{var};
$type = $input{type};

open( LS, "ls $info_dir |" );
@files = <LS>;
close LS;

@infolist=('common.idb');
foreach $file (@files) {
    $file =~ s/\s//g;
    if ( $file =~ /\w+\.idb$/ && $file !~ /common\.idb$/ ) {
	unshift(@infolist,$file);
    }
}

&init_html('CNSsolve Information Window');

&info_begin;

$found=0;
foreach $file ( @infolist ) {
    if ( $file =~ /common\.idb$/ || $file =~ /$type\.idb$/ ) {
	open (INFO, "< $info_dir/$file");
	while ( $line = <INFO> ) {
	    if ( $line =~ /vars:(.*)/ ) {
		@vars = split(' ',$1);
		foreach $var (@vars) {
		    if ( $variable =~ /^$var$/ ) {
			$found=1;
		    }
		}
		if ( $found ) {
		    while ( ($line = <INFO>) !~ /^#\s*$/ ) {
			if ( $line =~ /info:\s*(.*)/ ) {
			    print "$1\n";
			}
			else {
			    print "$line";
			}
		    }
		    close (INFO);
		    &info_end;
		    exit;
		}
	    }
	}
	close (INFO);
    }
}

if ( ! $found ) {
    print "No information was found for this input field\n";
}

exit;

sub info_begin {

    print '
  <BODY BGCOLOR=white>
  <CENTER><FONT SIZE=+1>CNSsolve Information</FONT></CENTER>
  <BR>
  ';

}

sub info_end {

    print '
  <BR> <FORM> <CENTER>
  <INPUT TYPE=BUTTON VALUE="Dismiss" onClick="self.close()">
  </CENTER> </FORM> </BODY> </HTML>
  ';

}

# EOF --> <H3>An error has occurred.</H3>
