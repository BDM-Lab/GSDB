#!/usr/bin/perl
# <!--
$this_cgi = "cns_jsedit.cgi";
$modified = "Last modified 05/26/98 J. Reklaw";
#
# Written by: Jesse Reklaw and Paul Adams
#
# copyright Yale University
#

require 'cns.lib';

$sel = &init_cgi();

%nv = &parse_args($sel);
$textarea = $nv{'field'};
$heading = $nv{'heading'};
$rows = $nv{'rows'};
$cols = $nv{'cols'};
&init_html("Text Editor");
print " 
  <SCRIPT LANGUAGE=\"JavaScript\">
  <!-- Begin hiding script contents from old browsers. 
    function updateTextArea() {
      window.opener.document.forms[0]['${textarea}'].value=
        document.forms[0].extraTA.value
      window.self.close()
    } 
    function writeTextArea() {
      document.write('<TEXTAREA NAME=extraTA ROWS=$rows COLS=$cols>')
      document.write(window.opener.document.forms[0]['${textarea}'].value)
      document.write('</TEXTAREA><BR>')
    }
  // End hiding script contents. --> </SCRIPT> 
  <BODY BGCOLOR=white>
  <H3> $heading </H3>
  <FORM>
  <SCRIPT> <!-- Begin hiding script contents from old browsers. 
    writeTextArea();
  // End hiding script contents. --> </SCRIPT> 
  <INPUT TYPE=BUTTON VALUE=Done onClick=\"updateTextArea()\">
  </FORM> </BODY> </HTML>\n";

# EOF --> <H3>An error has occurred.</H3>
