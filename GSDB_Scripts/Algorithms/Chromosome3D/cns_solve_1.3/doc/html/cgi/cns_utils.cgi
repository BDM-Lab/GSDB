#!/usr/bin/perl
# <!--
$this_cgi = "cns_utils.cgi";
$modified = "Last modified 12/08/98 P. Adams";
#
# CGI file for generating utility file listing
#
# Written by: Paul Adams
#
# copyright Yale University
#

require 'cns.lib';

### INITIALIZE SCRIPT AND HTML PAGE...

$query = &init_cgi();
$cns_view .= "?$query";
$cns_save .= "?$query";
$jscode  = &jscode_cnsView();
$jscode .= &jscode_cnsSave();
&init_html( 'CNSsolve Utility Files' , $jscode );
print "<BODY BGCOLOR=$CNS_BGCOLOR>";

### Get utility file names...

# PROBLEMS WITH PERL'S readdir(), USE SYSTEM CALL ls INSTEAD...
# KLUDGE! THIS WON'T WORK WITH SYMBOLIC LINKS!

print "<HR WIDTH=\"100%\">\n";

open( LS, "ls $util_dir|" );
@files = <LS>;
close LS;
	
print "<TABLE CELLSPACING=0 CELLPADDING=2 BORDER=2>\n";

print "<H2><U>Utility programs</U></H2>\n";

$count=0;
	
foreach $file (@files) {
	    
    if ( $count == 0 ) {
	print "<TR>\n";
    }
    
    $file =~ s/\s//g;	# Eliminate extraneous whitespace
    $size = -s "$util_dir/$subdir/$file";
    
    $size=int($size/1024)+1;
    $units="Kb";
    
    if ( $size > 999 ) {
	$size=int($size/1024)+1;
	$units="Mb";
    }

    print "<TD>\n";
    print "<TABLE CELLSPACING=0 CELLPADDING=2 BORDER=2>\n";
    print "<TR><TD><A HREF=\"javascript:cnsView(\'utils/$subdir/$file\')\">\n";
    print "<IMG SRC=\"$icon_dir/view_small.gif\" ALT=\"VIEW\"\n";
    print "     BORDER=0 HEIGHT=20 WIDTH=43 ALIGN=CENTER></A></TD>\n";
    print "<TD><A HREF=\"javascript:cnsSave(\'utils/$subdir/$file\')\">\n";
    print "<IMG SRC=\"$icon_dir/save_small.gif\" ALT=\"SAVE\"\n";
    print "     BORDER=0 HEIGHT=20 WIDTH=43 ALIGN=CENTER></A></TD>\n";
    print "<TD NOWRAP><B>$file</B>&nbsp;(<I>$size $units</I>)</TD></TR>\n";
    print "</TABLE>\n";
	    
    $count++;
	    
    if ( $count == 1 ) {
	print "</TR>\n";
	$count=0;
    }
	    
}

print "</TABLE>\n";

print "<HR WIDTH=\"100%\">\n";

### Print blank table to leave some space at end of page

print "<TABLE>\n";
foreach (1..12) { print "  <TR><TD>&nbsp;</TD></TR>\n"; }
print "</TABLE></BODY></HTML>\n";

# EOF --> <H3>An error has occurred.</H3>
