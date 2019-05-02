#!/usr/bin/perl
# <!--
$this_cgi = 'cns_menu.cgi';
$modified = 'Last modified 5/11/98 J. Reklaw';
#
# CGI file for generating input file menu
#
# Written by: Paul Adams and Jesse Reklaw
#
# copyright Yale University
#

require 'cns.lib';

##### PROGRAM FLOW OF CONTROL #####

$argstr=&init_cgi();
$arglen=$ENV{'CONTENT_LENGTH'};

# go through stdin if nothing is passed in the QUERY_STRING

if ( ! $argstr && ( $arglen >> 0 ) ) {
  
    while ($line = <STDIN>) {
	next if ($line =~ /^\s*$/);
	if ($line =~ /^\s*\-+(\d+)\s*$/) {
	    $file_id = $1;
	} 
    }
}

# take the arguments from the QUERY-STRING

else {
    %input = &parse_args ($argstr);
    $menu = $input{menu};
    $default = $input{default};
    &write_menu_page ($menu,$default);
}

exit;

###########################################################################

sub write_menu_page {
  local ($menu,$default) = @_;
  local($entry,@menulist,$def_str,$def_toggle);

  if ($default ne "") {
    $def_str = "\n       + '?default=${default}'\n       ";
    $def_toggle = 'on';
  } else {
    $def_str = '';
    $def_toggle = 'off';
  }

  open (DESC,"< $desc_file");
  while ( $entry = <DESC> ) {
    if ( $entry =~ /(\w+)\s+(\w*)\s*/ ) { 
      push(@menulist,$1); 
    }
  }
  close (DESC);

  if ( $menu eq "disabled" ) {
    $jscode = "      function Goto_Anchor(anchor) { }\n";
  } else {
    $jscode = "      function Goto_Anchor(anchor) {
        parent.text.location.href = '${cns_input}' $def_str + '#' + anchor
      }\n";
  }
  $jscode .= "      function Set_Defaults() {
        parent.text.location.href = '${cns_default}'
        parent.menu.location.href = '$cns_menu?menu=disabled'
      }
      function Reread_Frame() {
        parent.text.location.href = '$cns_input'
        parent.menu.location.href = '$cns_menu'
      }\n";

  &init_html('CNSsolve Menu',$jscode);

  print "<BODY BGCOLOR=$CNS_BGCOLOR>
    <CENTER><TABLE CELLSPACING=0 CELLPADDING=0>
    <TR ALIGN=CENTER VALIGN=CENTER>
      <TD> <IMG SRC=\"$icon_dir/input_files/input_files_title.gif\"
            HEIGHT=28 WIDTH=156 ALIGN=ABSCENTER ALT=\"Input Files\">
      </TD> </TR>
    <TR> 
      <TD> <A HREF=\"javascript:Goto_Anchor('own-file')\">
           <IMG SRC=\"$icon_dir/input_files/edit_file.gif\"
            BORDER=0 HEIGHT=28 WIDTH=156 ALIGN=ABSCENTER ALT=\"Edit File\"></A>
      </TD> </TR>\n";

    foreach $item (@menulist) {
      print "    <TR>\n      <TD>
           <A HREF=\"javascript:Goto_Anchor('${item}')\">
           <IMG SRC=\"$icon_dir/input_files/$item.gif\"
	    BORDER=0 HEIGHT=28 WIDTH=156 ALIGN=ABSCENTER ALT=\"$item\"> </A>
      </TD> </TR>\n";
    }

  print "
    <TR>
      <TD> <IMG SRC=\"$icon_dir/green.gif\"
            BORDER=0 HEIGHT=28 WIDTH=156 ALIGN=ABSCENTER> 
      </TD> </TR>
    <TR>
      <TD> <IMG SRC=\"$icon_dir/input_files/default_${def_toggle}_title.gif\"
            HEIGHT=28 WIDTH=156 ALIGN=ABSCENTER ALT=\"Defaults ${def_toggle}\">
      </TD> </TR>
    <TR>
      <TD> <A HREF=\"javascript:Set_Defaults()\">
           <IMG SRC=\"$icon_dir/input_files/default_set.gif\"
            BORDER=0 HEIGHT=28 WIDTH=156 ALIGN=ABSCENTER ALT=\"Set Defaults\"></A>
      </TD> </TR>
    <TR>
      <TD> <A HREF=\"javascript:Reread_Frame()\">
           <IMG SRC=\"$icon_dir/input_files/default_remove.gif\"
            BORDER=0 HEIGHT=28 WIDTH=156 ALIGN=ABSCENTER ALT=\"Remove Defaults\"> </A>
      </TD> </TR>
    </TABLE></CENTER> </BODY> </HTML>\n";

}

# EOF --> <H3>An error has occurred.</H3>
