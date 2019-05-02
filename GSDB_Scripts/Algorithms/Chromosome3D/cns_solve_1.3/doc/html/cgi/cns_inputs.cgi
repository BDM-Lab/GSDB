#!/usr/bin/perl
# <!--
$this_cgi = "cns_inputs.cgi";
$modified = "Last modified 5/29/98 P. Adams";
#
# CGI file for generating input file listing
#
# Written by: Jesse Reklaw and Paul Adams
#
# copyright Yale University
#

require 'cns.lib';

### INITIALIZE SCRIPT AND HTML PAGE...

$query = &init_cgi();
$cns_form .= "?$query&action=inp2form";
$cns_view .= "?$query";
$jscode  = &jscode_cnsView();
$jscode .= &jscode_cnsForm();
$jscode .= &jscode_cns_userForm();
&init_html( 'CNSsolve Script Files' , $jscode );
print "<BODY BGCOLOR=$CNS_BGCOLOR>";


### USER-FILE SUBMISSION FORM 

$form_name = 'user_form';

print "
  <P><A NAME=\"own-file\"></A>
     <HR WIDTH=\"100%\"></P>
  <CENTER><P><FONT SIZE=+3>Edit file</FONT></P></CENTER>
  <FORM NAME=\"$form_name\"
        METHOD=POST
        ACTION=\"$cns_form\"
        TARGET=\"$form_name\"
        ENCTYPE=\"multipart/form-data\">
    <P><FONT SIZE=+1>Enter the name of your CNSsolve input file:</FONT><BR>
       <INPUT NAME=\"cns_file\" TYPE=\"file\" SIZE=60></P>
    <TABLE WIDTH=\"100%\" CELLSPACING=5 CELLPADDING=0><TR><TD ALIGN=RIGHT>
       <TABLE CELLSPACING=0 CELLPADDING=2 BORDER=2><TR>
       <TD><A HREF=\"javascript:cns_userForm(\'$form_name\',\'$cns_form\',\'off\')\">
       <IMG SRC=\"$icon_dir/edit_small.gif\" ALT=\"EDIT\"
        BORDER=0 HEIGHT=20 WIDTH=43 ALIGN=CENTER></A></TD>
       <TD><A HREF=\"javascript:cns_userForm(\'$form_name\',\'$cns_form\',\'on\')\">
       <IMG SRC=\"$icon_dir/edit_help_small.gif\" ALT=\"EDIT+HELP\"
        BORDER=0 HEIGHT=20 WIDTH=75 ALIGN=CENTER></A></TD>
       </TR></TABLE>
     </TD></TR></TABLE>
  </FORM>
";


### GET SUBDIRECTORIES AND INPUT FILE NAMES...

# PROBLEMS WITH PERL'S readdir(), USE SYSTEM CALL ls INSTEAD...
# KLUDGE! THIS WON'T WORK WITH SYMBOLIC LINKS!
open( LS, "ls $input_dir|" );
@subdirs = <LS>;
close LS;
if ($#subdirs < 0) { &die_html("No directories read from $input_dir"); }
foreach $subdir (@subdirs) {
  if ($subdir =~ m/^l.*\-\>/) { &die_html("Error: will not follow link"); }
  $subdir =~ s/\s//g;	# ELIMINATE EXTRANEOUS WHITESPACE
  if ($subdir !~ /_data$/i) {
    $file_info{$subdir} = '';
    open( LS, "ls $input_dir/$subdir|" );
    @files = <LS>;
    close LS;

### NOW PARSE INPUTS AND STORE:
###   - view button
###   - edit button
###   - file name
###   - file description

    foreach $file (@files) {
      if ($file =~ m/\.inp$/i || $file =~ m/\.def$/i || $file =~ m/\.sdb$/i) {
        $file =~ s/\s//g;	# ELIMINATE EXTRANEOUS WHITESPACE
        $desc = '';
        open( FILE, "< $input_dir/$subdir/$file" ) 
          || &die_html("Cannot read file ($input_dir/$subdir/$file)");
	$cont = 1;
        while( $cont ) {
	  $_=<FILE>;
          if (/^\s*\{\+\s*file\:\s*(\S+)/i) { 
          # we could do some kind of error checking here to make sure
          #   we've got the right file.
          }
  # KLUDGE! USE APPEND_TEXT_UNTIL ?
          elsif (/^\s*\{\+\s*description\:\s*(.*)$/i) { 
            $desc = $1;
            while ($cont) {
              $desc =~ s/\s+/ /g; # SQUASH ALL WHITESPACE INTO A SINGLE SPACE
              if ($desc =~ m/^(.*)\+\}/) {
                $desc = $1;
                $cont = 0;
              } else {
                $desc .= <FILE>;	# APPEND NEXT LINE OF FILE
                if (eof) 
                  { &die_html("EOF detected while trying to parse $file"); }
              }
            } 
          }#elsif ...
        }#while ...
        close FILE;
        $file_info{$subdir} .= "<TR>
          <TD ROWSPAN=2 BGCOLOR=#777777>&nbsp</TD>
          <TD ALIGN=LEFT><FONT SIZE=+1>$file</FONT></TD>
          <TD ALIGN=RIGHT>
              <TABLE CELLSPACING=0 CELLPADDING=2 BORDER=2><TR>
              <TD><A HREF=\"javascript:cnsView(\'inputs/$subdir/$file\')\">
              <IMG SRC=\"$icon_dir/view_small.gif\" ALT=\"VIEW\"
               BORDER=0 HEIGHT=20 WIDTH=43 ALIGN=CENTER></A></TD>
              <TD><A HREF=\"javascript:cnsForm(\'$subdir/$file\',\'off\')\">
              <IMG SRC=\"$icon_dir/edit_small.gif\" ALT=\"EDIT\" 
               BORDER=0 HEIGHT=20 WIDTH=43 ALIGN=CENTER></A></TD>
              <TD><A HREF=\"javascript:cnsForm(\'$subdir/$file\',\'on\')\">
              <IMG SRC=\"$icon_dir/edit_help_small.gif\" ALT=\"EDIT+HELP\" 
               BORDER=0 HEIGHT=20 WIDTH=75 ALIGN=CENTER></A></TD>
              </TR></TABLE>
          </TD></TR>
          <TR><TD COLSPAN=2>$desc</TD></TR>\n";
      }#if ($file =~ ...
    }#foreach $file ...
  }#if ($subdir !~ /_data$/i) ...
}#foreach $subdir ...


### PRINT OUT ACCORDING TO DESCRIPTION FILE
open( DESC, "< $desc_file" );
$new = 1;
print "<HR WIDTH=\"100%\">\n";
while($line=<DESC>) {
  ($subdir, $heading) = split(/\s+/,$line,2);
  if ($subdir) {
    print "<P><A NAME=\"$subdir\"></A>\n";
    $heading =~ tr/"//d;
    if (!$heading) { $heading = $subdir; }
    if ($new) { 
      print "<CENTER><P><FONT SIZE=+3>$heading</FONT></CENTER>\n";
    } else {
      print "<BR><FONT SIZE=+2>$heading</FONT>\n";
    }
    if ($file_info{$subdir}) {
      print "<TABLE CELLSPACING=5 CELLPADDING=0 WIDTH=100%>
        $file_info{$subdir}</TABLE>\n";
      $file_info{$subdir} = '';
    }
    $new = 0;
  } else {
    print "<HR WIDTH=\"100%\">\n";
    $new = 1;
  }
}

### PRINT OUT OTHERS NOT IN DESCRIPTION FILE
foreach $other (keys %file_info) {
  if ($file_info{$other}) {
    print "<HR WIDTH=\"100%\">\n";
    print "<CENTER><P><FONT SIZE=+3>$other</FONT></CENTER>\n";
    print "<TABLE CELLSPACING=5 CELLPADDING=0>
      $file_info{$other}</TABLE>\n";
    $file_info{$other} = '';
  }
}

### PRINT BLANK TABLE TO LEAVE SOME SPACE AT END OF PAGE
print "<TABLE>\n";
foreach (1..12) { print "  <TR><TD>&nbsp;</TD></TR>\n"; }
print "</TABLE></BODY></HTML>\n";

# EOF --> <H3>An error has occurred.</H3>
