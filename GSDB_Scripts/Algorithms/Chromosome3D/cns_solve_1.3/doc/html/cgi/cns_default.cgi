#!/usr/bin/perl
# <!--
$this_cgi = 'cns_default.cgi';
$modified = 'Last modified 5/29/98 P. Adams';
#
# CGI file for setting default values 
#
# Written by: Paul Adams and Jesse Reklaw
#
# copyright Yale University
#

require 'cns.lib';

##### PROGRAM FLOW OF CONTROL #####

&init_cgi();
$arglen=$ENV{'CONTENT_LENGTH'};

# go through the input if something is passed on stdin
if ( $arglen >> 0 ) {

    while ($line = <STDIN>) {
	next if ($line =~ /^\s*$/);
	if ($line =~ /^\s*\-+(\d+)\s*$/) {
	    $file_id = $1;
	} 
	elsif ($line =~ /name="default_file"/) {
	    %input = &get_vars_from_file;
	    &write_form_page (%input);
	}
	elsif ($line =~ /default_values=/) {
	    %input = &parse_args ($line);
	    foreach $name (keys %input) {
		if ( $name eq save_default_to_file ) {
		    print "Content-type: application/octet-stream name=project.def\n";
		    print "Content-disposition: attachment; filename=project.def\n\n";
		    &save_to_file ('stdout', %input);
		} 
		elsif ( $name eq continue_edit ) {
		    $def_file = $tmp_dir . 
                                "/${full_date}_${remote_host}_defaults";
		    &save_to_file ($def_file, %input);
		    &write_frame($def_file);
		}
	    }
	}
    }
}

# if no stdin then we are just meant to print out the standard form

else { 
  %input = {};
  &write_form_page(%input); 
}

exit;

###########################################################################

sub get_vars_from_file {
    local (%nv);
    while ( <STDIN> ) {
        ### if we have any variable(s) just on one line...
	if (/^\s*\{===>\}\s*(\S.*)\;\s*$/) {
	    @varpairs = split(';',$1);
	    foreach $varpair (@varpairs) {
		($lhs,$rhs) = split('=',$varpair,2);
                if ($rhs =~ /^\s*\"(.*)\"\s*$/) { $rhs = $1; }
		$nv{$lhs} .= $rhs;		
	    }
	}
        ### if we have a multi-line variable...
	elsif (/^\s*\{===>\}\s*([^\s=]+)\s*=\s*(\S.*)$/) {
	    ($lhs,$rest) = ($1,$2);
	    ($rhs,$remainder) = &append_text_until( STDIN, '\;', $rest );
            if ($rhs =~ /^\s*\"(.*)\"\s*$/) { $rhs = $1; }
	    $nv{$lhs} .= $rhs;		
	}
    }
    return %nv;
}

sub write_form_page {
  local (%input) = @_;
  local($jscode,$def_str,$a); 

  if ($default ne "") { $def_str = "?default=$default"; }
  $jscode = "    function Cancel() {
      parent.text.location.href = \"$cns_input\";
      parent.menu.location.href = \"$cns_menu\";
    }\n";
  $jscode .= &jscode_confirmedReset();
  $jscode .= &jscode_assertReal();
  &init_html('CNSsolve Default Setup',$jscode);

  print "<BODY BGCOLOR=$CNS_BGCOLOR>
    <HR WIDTH='100%'>
    <CENTER><P><FONT SIZE=+3>Read defaults from file</FONT></P></CENTER>
    <FORM NAME='userUpload' METHOD=POST ACTION=\"$cns_default$def_str\"
     ENCTYPE='multipart/form-data'>
    <TABLE> <TR>
      <TD> <INPUT NAME=default_file TYPE=FILE SIZE=30></TD>
      <TD> <INPUT TYPE=SUBMIT VALUE='Read file'></TD>
    </TR> </TABLE>
    </FORM>

    <HR WIDTH='100%'>
    <CENTER><P><FONT SIZE=+3>Current defaults</FONT></P></CENTER>
    <CENTER>
    <FORM NAME='inputFields' METHOD=POST ACTION=\"$cns_default$def_str\">
    <INPUT TYPE=HIDDEN NAME='default_values' VALUE=''>
    <TABLE CELLPADDING=3 BORDER=2 BGCOLOR=#99AA99>
    <TR>
      <TD ALIGN=RIGHT><FONT SIZE=+1>spacegroup</FONT></TD>
      <TD COLSPAN=3><INPUT NAME=sg VALUE=\"$input{sg}\"></TD>
    </TR>
    <TR>
      <TD ALIGN=RIGHT ROWSPAN=2><FONT SIZE=+1>unit cell</FONT></TD>
    ";

  foreach $a ('a','b','c') {
    print "      <TD ALIGN=CENTER>$a<BR><INPUT NAME=$a VALUE='$input{$a}' 
        SIZE=8 onChange=\"assertReal(this,'$input{$a}')\"></TD>\n";
  }

  print "        </TR><TR>\n";

  foreach $a ('alpha','beta','gamma') {
    print "      <TD ALIGN=CENTER>$a<BR><INPUT NAME=$a VALUE='$input{$a}' 
        SIZE=8 onChange=\"assertReal(this,'$input{$a}')\"></TD>\n";
  }

  print "    </TR>
    <TR>
      <TD ALIGN=RIGHT><FONT SIZE=+1>anomalous file</FONT></TD>
      <TD COLSPAN=3><INPUT NAME=anom_library 
          VALUE=\"$input{anom_library}\" SIZE=30></TD>
    </TR>
    </TABLE>
    </CENTER>

    <HR WIDTH='100%'>
    <TABLE>
    <TR>
      <TD><INPUT TYPE=SUBMIT NAME='continue_edit'
           VALUE='Start editing files'></TD>
      <TD><INPUT TYPE=SUBMIT NAME='save_default_to_file' 
           VALUE='Save defaults to file'></TD>
      <TD><INPUT TYPE=BUTTON VALUE='Reset' 
           onClick=\"confirmedReset(\'inputFields\')\"></TD>
      <TD><INPUT TYPE=BUTTON VALUE='Cancel' 
           onClick=\"Cancel()\"></TD>
    </TABLE> </FORM> </BODY> </HTML>\n";
}

sub save_to_file {
    local ($out_file, %nv) = @_;

    if ( $out_file ne 'stdout' ) {
	open (SAVEOUT,">&STDOUT");
	open (STDOUT,">$out_file");
    }

    print "{+ description: defaults file written on $full_date +}\n";
    print '
{====================== crystallographic data ========================}

{* space group *}
{* use International Table conventions with subscripts substituted
     by parenthesis *}
';
print "{===>} sg=\"$nv{sg}\";\n";
print '
{* unit cell parameters in Angstroms and degrees *}
{+ table: rows=1 "cell" cols=6 "a" "b" "c" "alpha" "beta" "gamma" +}
';
print "{===>} a=$nv{a};\n";
print "{===>} b=$nv{b};\n";
print "{===>} c=$nv{c};\n";
print "{===>} alpha=$nv{alpha};\n";
print "{===>} beta=$nv{beta};\n";
print "{===>} gamma=$nv{gamma};\n";
print '
{* anomalous f\' f\'\' library file *}
{* If a file is not specified, no anomalous contribution will be included *}
{+ choice: "CNS_XRAYLIB:anom_cu.lib" "CNS_XRAYLIB:anom_mo.lib" "" user_file +}
';
print "{===>} anom_library=\"$nv{anom_library}\";\n";
print '
{===========================================================================}
';

    if ( $out_file ne 'stdout' ) {
        close (STDOUT);
        open (STDOUT,">&SAVEOUT");
    }

}

sub write_frame {
  local ($default_fname) = @_;
  local($jscode,$def_str); 

  $def_str = "?default=$default_fname";
  $jscode = "    function Load_Pages() {
      parent.menu.location.href = '${cns_menu}${def_str}';
      parent.text.location.href = '${cns_input}${def_str}';
    }\n";
  &init_html('',$jscode);
  print "<BODY BGCOLOR=$CNS_BGCOLOR onLoad=\"Load_Pages()\"></BODY></HTML>\n";
}

# EOF --> <H3>An error has occurred.</H3>
