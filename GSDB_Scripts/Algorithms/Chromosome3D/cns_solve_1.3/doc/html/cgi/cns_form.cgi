#!/usr/bin/perl
# <!--
$this_cgi = 'cns_form.cgi';
$modified = 'Last modified 98may26 J. Reklaw';
#
# CGI file for parsing CNSsolve task file, creating form page and applying form
# to original input file
#
# Written by: Jesse Reklaw and Paul Adams
#
# copyright Yale University
#

require 'cns.lib';
 
# note to Paul -- I've marked everything that still needs a fix with "KLUDGE"
# DO WE EVER CONSIDER EMPTY FORMS? (KLUDGE)

########## GLOBAL DEFINITIONS ##########

  $NUMBERED = 'numbered';	# user keyword for numbered rows and cols

# FORMATTING PREFERENCES:

  $COLOR='#698B69';	# color of the form section title blocks (medium green)
  $INFOCOL='#DFDFDF';	# color of information header (gray)
  $ALIGN='RIGHT';	# align documentation right
  $DOCMAX=300;	# if documentation > DOCMAX then fill whole row (was 200)
  $A_FEW=3;	# if choices <= A_FEW then selection list else radio buttons
  $RADIO_LONG=45; # maximum total number of characters for all radio buttons
  $BIG_SIZE=30;	# default size for big string fields (filenames) -- was 50
  $STR_SIZE=10;		# default size for string fields (was 15)
  $NUM_SIZE=8;		# default size for number fields
  $SIZE_INC=5;		# extra space for sizes larger than max
  $TEXTAREA_COLS= 60;
  $MIN_TEXTAREA_ROWS= 2;
  $AUX_TEXTAREA_ROWS= 15;
  $COL1=''; $COL2='COLSPAN=2'; $COL3='COLSPAN=3';
  $SPACER1='<IMG SRC="' . "$icon_dir/whiteblock.gif" . '" HEIGHT=1 WIDTH=100>';
  $BORDER = 1;		# border width around tables?
  $INFOWIN_ICON = '=';	# text printed on buttons that call up infoWindow

# CHECK TO SEE IF BROWSER CAN HANDLE NESTED TABLES. NOT SURE WHAT ALL OF 
# THESE STRINGS COULD BE, PARTICULARLY THE MSIE ONES.

  $browser = $ENV{'HTTP_USER_AGENT'};
  if ($browser=~m/Mozilla\/(\d+)\./i && $1 >= 3) { $ok_to_nest = 1; }
  elsif ($browser=~m/\(\s*compatible\s*\;\s*MSIE\s*(\d+)\./i) {	# KLUDGE
    if ($1 >= 2) { $ok_to_nest = 1; }
    else { $ok_to_nest = 0; }
  } else { $ok_to_nest = 0; }
  if (!$ok_to_nest) { $browser_warning = "<H3><FONT COLOR=RED>
    The CNSsolve Forms interface has not been tested on your browser
    ($browser). Netscape3.0 or better is recommended.</FONT></H3>\n"; }
  else { $browser_warning = ''; }

##############################


##### PROGRAM FLOW OF CONTROL #####

$argstr = &init_cgi();

# For backwards compatibility, if there is no equals sign ('=') 
#   in the $argstr, then the $argstr is the $inp_file.
#   NOTE: the assumption is that no filename can have '=' in it.
if ($argstr !~ m/\=/) { 
  $action   = 'inp2form';
  $inp_file = $argstr; 
  $def_file = '';
}
else {
  %argnv = &parse_args($argstr);
  $action   = $argnv{'action'};
  $inp_file = $argnv{'file'};
  if (!$inp_file) { $inp_file = $argnv{$INP_FILE}; }
  $def_file = $argnv{'default'};
  if ($argnv{'cns_help'} =~ /on/i) { $print_info_buttons=1; }
}
if ( $inp_file !~ /^\/\w/ && $inp_file) { # shortcut for standard inp files
  $inp_file = "$input_dir/" . $inp_file; 
}
# dont open anything with an absolute path
elsif ( $action == 'inp2form' && $inp_file =~ /^\/\w/ ) {
    &die_html( "cannot open file $inp_file" );
}
# remove special characters from user-input pathname
$inp_file =~ s/[\~\`\!\$\&\;]//g;
$inp_file =~ s/\.\./\./g;

# remove any ^M and newlines at end of file names
$inp_file =~ s/\cM//g;
$inp_file =~ s/\n$//;
$def_file =~ s/\cM//g;
$def_file =~ s/\n$//;

undef( %argnv );
undef( $argstr );

if ($action eq 'form2inp') { 
  &form2inp();		# previous program is now a subroutine
}

else {			# action should be 'inp2form'	
  ### If no inp file pathname is given, 
  ###   then we need to parse it from stdin and save a temporary copy.
  if (!$inp_file) {

    # make temporary directory if not present. assumes /usr/tmp is open acccess
    if ( ! opendir TMPDIR, "$tmp_dir" ) {
      if ( system ("mkdir -p $tmp_dir") != 0 ) 
        { &die_html("could not make temporary directory ($tmp_dir)"); }
    }

    # get rid of out of date temporary files - 
    #   3 days for files, 1 day for empty dirs
    system 
 ("find $tmp_dir -type f -atime +3 '!' -name access.log -exec rm -f '{}' ';'");
    system 
 ("find $tmp_dir -type d -atime +1 depth -exec rmdir '{}' ';' 2>/dev/null");

    $inp_file = $tmp_dir . "/${full_date}_${remote_host}";
    while ($line = <STDIN>) {
      next if ($line =~ /^\s*$/);			# skip blank lines
      if ($line =~ /^\s*\-+(\d+)\s*$/) {		# if delimiter line for Netscape 
          $file_id = $1;				# ... get file ID
# modification ATB 06/03/08
      } elsif ($line =~ /-----------------------------/)  {  # if delimiter line for Internet Explorer 
          $file_id = $line;
          $file_id =~ s/-//g;  
          $file_id =~ s/\cM//g;
          $file_id =~ s/\n$//;
      } elsif ($line =~ /------WebKitFormBoundary/)  {  # if delimiter line for Safari...
          $file_id = $line;
          $file_id =~ s/WebKitFormBoundary//g;
          $file_id =~ s/-//g;
          $file_id =~ s/\cM//g;
          $file_id =~ s/\n$//;
# end modification
      } elsif ($line =~ /name="cns_file"/) {	# then get entire file...
	  $fname = &get_file( $file_id, $line, $inp_file );
	  last;
      } elsif ($line =~ /cns_edit_file=/) {
	  %input=parse_args($line);
	  $line_out = $input{cns_edit_file};
	  $line_out = unescape($line_out);
	  $line_out =~ s/\cM//g;
	  $line_out =~ s/\n\n$/\n/;
	  open (OUT_FILE, "> $inp_file");
	  print OUT_FILE "$line_out";
	  close OUT_FILE;         
	  if ( $input{cns_edit_def} ) {
	      if (!$def_file) {
		  $def_file = $inp_file . "_default";
	      }
	      $line_out = $input{cns_edit_def};
	      $line_out = unescape($line_out);
	      $line_out =~ s/\cM//g;
	      $line_out =~ s/\n\n$/\n/;
	      open (OUT_FILE, "> $def_file");
	      print OUT_FILE "$line_out";
	      close OUT_FILE;         
	  }
	  if ( $input{cns_edit_fname} ) {
	      $fname = $input{cns_edit_fname};
	      $fname =~ s/\cM//g;
	      $fname =~ s/\n$//;
	  }
	  else {
	      $fname = 'cns_edit';
	  }
      }
    }
    if ($fname !~ /^ERR:/) {
      &log( $tmplog, "$remote_host $fname R" );	# log the transaction
   } else { &die_html( "no file received ($fname)" ); }	
  } 

  else {
    ### $fname was not acquired from client; get from $inp_file .
    $fname = $inp_file;
    &log( $tmplog, "$remote_host $fname L" );	# log the transaction
  }
  @dirs_and_file = split(/[\/\\]+/,$fname);
  $fname = $dirs_and_file[$#dirs_and_file];
  undef( @dirs_and_file );
  if ($fname =~ /^\s*$/) { &die_html( "no filename given ($inp_file)" ); }	

  ### Print header, JavaScript, initial body & form tags, etc.
  $jscode  = &jscode_assertReal();
  $jscode .= &jscode_confirmedReset();
  $jscode .= &jscode_editWindow();
  $jscode .= &jscode_cnsInfo();
  &init_html("Edit $fname",$jscode);
  print "
    <BODY BGCOLOR=white> 
    $browser_warning
    $noscript_warning
    <H1><CENTER>$fname</CENTER></H1>
    <FORM ACTION=\"${this_cgi}?action=form2inp\" NAME=f1 METHOD=POST>
    <INPUT TYPE=HIDDEN NAME=\"$INP_FILE\" VALUE=\"$inp_file\">
    <INPUT TYPE=HIDDEN NAME=\"$TRUE_FNAME\" VALUE=\"$fname\">\n";

  ### Read, parse and set defaults
  if ($def_file) {
    print "<CENTER><TABLE><TR><TD BGCOLOR=#FF0000><FONT SIZE=+1>\n";
    print "Using defaults</FONT></TD></TR></TABLE></CENTER>\n";
    open( INP, "< $def_file" ) || &die_html( "cannot read $def_file" );
    while (<INP>) {
      $current_line_number++;
      if (/^\s*\{===>\}\s*(.*)$/) {
	@vars = &get_vars(INP,$1);
        foreach $var (@vars) { 
	  ($lhs,$rhs) = split('.cns_var_is_set_to.',$var);
	  if ( defined($default_value{$lhs}) ) {
	      &die_html( "error in defaults file at line $current_line_number:
                          parameter $lhs is multiply defined" );
	  }
	  $default_value{$lhs} = $rhs; 
	}
      }
    }
    close INP;
  }

  ### Read and parse '.inp' file.
  @form = ('info');		# initialize @form for the information header
  open( INP, "< $inp_file" ) || &die_html( "cannot read $inp_file" );
  while ($line=<INP>) { $current_line_number++; process_line($line); }
  &output_form();		# output the remainder 
  close INP;
  print '
    <INPUT TYPE=SUBMIT NAME="submit_view" VALUE="View updated file"> 
    <INPUT TYPE=SUBMIT NAME="submit_download" VALUE="Save updated file">
    <INPUT TYPE=BUTTON VALUE="Reset values" onClick="confirmedReset(\'f1\')">';
  print "\n</FORM>";
  print "<DIV ALIGN=right>";
  print "<B><U>CNSsolve HTML interface</U>&nbsp;&nbsp;</B><BR>";
  print "Copyright &copy; Yale University&nbsp;&nbsp;</DIV>\n";
  print "</BODY></HTML>\n";
}
# END FLOW OF CONTROL


#################### PARSING SUBROUTINES ###################
sub process_line { local($line) = @_;
  local($str,$remainder,$numrows,$rownames,$numcols,$colnames,$type);
  $_ = $line;
  if (/^\s*$/) {				# blank lines break doc's
    if ($form[$#form - 1] eq 'doc') { push(@form, 'break'); } 
  }
  elsif (/^remarks$/i) { }			# skip remarks
  elsif (/^\s*\{\-(.*)$/) {			# skip guidelines
    ($str,$remainder) = &append_text_until( INP, '\-\}', $1 ); 
    &process_line($remainder);
  }

  ### if we start a new heading...
  elsif (/^\s*\{====+\s*([^=]*)\s*=+===\}\s*$/) {
    $str = $1;
#    $str =~ tr/=//d;			# remove any spurious '=' chars
    &output_form(); 			# output previous table section
    @form = ('form',$str);		# begin new table section
  }

  ### if we begin documentation...
  elsif (/^\s*\{\*(.*)$/) {	
    ($str,$remainder) = &append_text_until( INP, '\*\}', $1 );
#    &process_line($remainder);
    $str =~ s/\s+/ /g;
    if ($str ne ' ' && $str ne '') {	# only accept non-blank lines
      $str =~ s/\&/\&amp\;/g;		# make HTML-friendly...
      $str =~ s/</\&lt\;/g;
      $str =~ s/>/\&gt\;/g;
      push(@form, 'doc', $str);
    } else { push(@form, 'break' ); }
  }

  ### if we have a predefined option list...
  elsif (/^\s*\{\+\s*choice:(.*)$/) {
    ($str,$remainder) = &append_text_until( INP, '\+\}', $1 );
    $str =~ s/choice://g;	# delete 'choice:' from list-string
    $str =~ s/\n/ /g;
    push(@form, 'varchoice', $str);
    &process_line($remainder);
  }

  ### if we have a table...
  elsif (/^\s*\{\+\s*table:(.*)$/) {
    ($str,$remainder) = &append_text_until( INP, '\+\}', $1 );
    $str =~ s/\n/ /g;
    if ($str =~ /\s*rows\s*=\s*(\d+)(.*)cols\s*=\s*(\d+)(.*)$/i) {
      ($numrows,$rownames,$numcols,$colnames) = ($1,$2,$3,$4);
    } elsif ($str =~ /\s*cols\s*=\s*(\d+)(.*)rows\s*=\s*(\d+)(.*)$/i) {
      ($numcols,$colnames,$numrows,$rownames) = ($1,$2,$3,$4);
    } else { $numrows=-1; $numcols=-1; }	# BAD TABLE PARSE
    push(@form, 'vartable', $numrows, $numcols, $rownames, $colnames);
    &process_line($remainder);
  }

  ### if we have special documentation...
  elsif (/^\s*\{\+\s*([^:]+):\s*(.*)$/) {
    $type = $1;
    ($str,$remainder) = &append_text_until( INP, '\+\}', $2 );
    if ($str =~ m/^\s*(\S.*)$/) { $str = $1; } # clear beginning spaces
    if ($str =~ m/^(.*\S)\s*$/) { $str = $1; } # clear ending spaces         
    $str =~ s/\&/\&amp\;/g;		# make HTML-friendly...
    $str =~ s/</\&lt\;/g;
    $str =~ s/>/\&gt\;/g;
    if ($type =~ /^description/i) {
      push(@form, 'desc', $str);
    } elsif ($type =~ /^directory/i) {
      $familyName = $str;
    } elsif ($type =~ /^author/i) {
      push(@form, 'auth', $str);
    } elsif ($type =~ /^copyright/i) {
      push(@form, 'copy', $str);
    } elsif ($type =~ /^reference/i) {
      push(@form, 'ref', $str);
    } elsif ($type =~ /^comment/i) {
      push(@form, 'com', $str);
    } elsif ($type =~ /^list/i) {
      push(@form, 'doc', "<PRE>$str</PRE>");
    } else { push(@form, "ERR$str" ) ;}
#    &process_line($remainder);
  }

  ### if we have any variable(s) ...
  elsif (/^\s*\{===>\}\s*(.*)$/) {
    @vars = &get_vars(INP,$1);
    foreach $var (@vars) { 
      ($lhs,$rhs) = split('.cns_var_is_set_to.',$var);
      if ( defined($varstore{$lhs}) ) {
	  &die_html( "error in input file at line $current_line_number:
                      parameter $lhs is multiply defined" );
      }
      else {$varstore{$lhs} = 1}
      # NOTE THAT DEFAULT REPLACEMENT WILL ONLY WORK ON '_paren' AND
      # '_quote' VARIABLES IF THE DEFAULT VALUE WAS ALSO ORIGINALLY
      # IN PARENS OR QUOTED.
      if ($default_value{$lhs} ne '') { $rhs = $default_value{$lhs}; }
      # if $lhs is a textarea or is marked for parenthesis then it must
      # be a multi-line variable, PDA
      if ($lhs =~ m/$PARENFLAG$/ || $lhs =~ m/^$TEXTAREAROOT\d+$/) 
        { $type = 'varlong'; }
      # if $rhs contains a newline it must be a multi-line variable, PDA
      elsif ( $rhs =~ m/\n/ )
        { $type = 'varlong'; }
      else  { $type = 'var'; }
      if ($lhs =~ m/$QUOTEFLAG$/) { $size = $STR_SIZE; }
      else  { $size = $NUM_SIZE; }
      if (length($rhs) > $size) { $size =length($rhs)+$SIZE_INC; }
      push(@form, $type, $size, $lhs, $rhs); 
    }
  }

  ### default
  else {}
}

sub output_form { 
#print join(' | ',@form); 	# UNCOMMENT TO DEBUG FORM PARSING
  local ($desc_str,$com_str,$auth_str,$ref_str,$tabcol);
  local ($i,$item,$doc_str,$lhs,$rhs,$q,$size,$cho_len,$max_cho_len);
  local ($choice_str,$other,@choices,$choice,@lines,$lines);

  ## I. HEADER INFORMATION

  if ($form[0] eq 'info') {
    $desc_str=''; $com_str=''; $auth_str=''; $ref_str='';
    while ($form[++$i]) { 
      if    ($form[$i] eq 'desc') { $desc_str .= " $form[++$i]"; }
      elsif ($form[$i] eq 'com')  { $com_str  .= $form[++$i]; }
      elsif ($form[$i] eq 'auth') { $auth_str .= " $form[++$i]"; }
      elsif ($form[$i] eq 'copy') { $copy_str .= " $form[++$i]"; }
      elsif ($form[$i] eq 'ref')  { $ref_str  .= "  <LI>$form[++$i]\n"; }
    }
    print "<TABLE WIDTH=\"100%\" CELLPADDING=10><TR><TD BGCOLOR=$INFOCOL>\n";
    if ($desc_str) { print "<BIG>$desc_str</BIG><BR>\n"; }
    if ($com_str)  { print "<PRE>$com_str</PRE>\n"; } 
    if ($auth_str) 
      {print "<BR><BIG><B>Authors:</B></BIG>$auth_str<BR>\n";}
    if ($copy_str) 
      {print "<BR><BIG><B>Copyright &copy; </B></BIG>$copy_str<BR>\n";}
    if ($ref_str)  { print "<H3>References</H3><UL>\n$ref_str  </UL>\n"; }
    print "</TD></TR></TABLE>\n";

  ## II. DATA FORM

  } elsif ($form[0] eq 'form' && $form[2]) { # if !$form[2] then form is blank

  ## A. HEADING

    print "<TABLE BORDER=$BORDER WIDTH=\"100\%\">\n";
    if ($form[1]) {
      print "  <TR><TD $COL3 BGCOLOR=$COLOR><CENTER>\n";
      print "    <FONT FACE=helvetica COLOR=white><BIG><B>$form[1]\n";
      print "    </B></BIG></FONT></CENTER></TD></TR>\n";
    } else {
      print "  <TR><TD $COL3 BGCOLOR=$COLOR><CENTER>\n";
      print "    <FONT SIZE=4 FACE=helvetica COLOR=white>&nbsp;\n";
      print "    </FONT></CENTER></TD></TR>\n";
    }
    print "<TR><TD>$SPACER1<TD>$SPACER1<TD>$SPACER1</TR>\n";

    $i = 2;
    while ($item=$form[$i++]) {


  ## B. DOCUMENTATION

      if ($item eq 'doc') {
        $doc_str = $form[$i++];
        $this_doc = $doc_str;
        if ($doc_str=~m/file\(?s?\)?\s*$/) { $bigstr=1; } else { $bigstr=0; }
        $doc_str = "    <B>$doc_str</B><BR>\n";
        while ($form[$i] eq 'doc') { 
          $doc_str .= "    <I>$form[$i+1]</I>\n"; 
          $i+=2;
        }
        if ($form[$i] eq 'var' || $form[$i] eq 'varchoice') { 
          if ($form[$i+1] < $BIG_SIZE && !$bigstr) {
            print "  <TR><TD $COL2 ALIGN=$ALIGN>\n$doc_str";
            print "    </TD><TD $COL1>\n";
            $tabcol=1;
          } else {
            if (length($doc_str) > $DOCMAX) {
              print "  <TR><TD $COL3>\n$doc_str\n<BR>$SPACER1\n";
              $tabcol=3;
            } else {
              print "  <TR><TD $COL1 ALIGN=$ALIGN>\n$doc_str";
              print "    </TD><TD $COL2>\n";
              $tabcol=2;
            }
          }
        } else {
          print "  <TR><TD $COL3>\n$doc_str\n";	# for vartables only ???
          $tabcol=3;
        }


  ## C. CHOICES AND TABLE FORMATTING

      } elsif ($item eq 'varchoice') {
        $choice_str = $form[$i++];
        while ($form[$i] eq 'break') { $i++; }
        if ($form[$i] !~ /^var/) { print "ERROR!</TD></TR>\n"; }

      } elsif ($item eq 'vartable') {
        if ($ok_to_nest) {
          print "<BLOCKQUOTE><TABLE BORDER=$BORDER WIDTH=\"95%\">\n";
        } else {
          print "</TD></TR></TABLE>\n";		# break for inner table
          print "<TABLE BORDER=$BORDER WIDTH=\"100%\">\n";
        }
        ($numrows, $numcols) = ($form[$i++], $form[$i++]);
        ($rownames, $colnames) = ($form[$i++], $form[$i++]);
        if ($rownames=~m/$NUMBERED/i) 
          { $rownames = '"' . join('" "',1..$numrows) . '"' ; }
        if ($colnames=~m/$NUMBERED/i) 
          { $colnames = '"' . join('" "',1..$numcols) . '"' ; }
        $nested_table=1;	# /set flag that we're in a nested table
        if ($colnames=~m/^\s+$/) { $colnames = ''; }
        if ($rownames=~m/^\s+$/) { $rownames = ''; }
        $colnames =~ s/\n/ /g;
        $rownames =~ s/\n/ /g;

        if ($colnames) {
          foreach $i (0..$numcols-1) {
            if ($colnames=~m/^\s*"([^"]*)"(.*)$/) {
              ($name,$colnames) = ($1, $2);
            } else { $name = ''; }
            $colnames[$i] = $name;
          }
          if ($rownames) { print "  <TR><TH>&nbsp;</TH>"; }
          else { print "  <TR>"; }
          foreach $name (@colnames) { print "<TH>$name</TH>"; }  
          print "    </TR>\n";
        }
        $colcount = 1; # colcount is 1-based

        if ($rownames) {
          foreach $i (0..$numrows-1) {
            if ($rownames=~m/^\s*"([^"]*)"(.*)$/) {
              ($name,$rownames) = ($1, $2);
            } else { $name = ''; }
            $rownames[$i] = $name;
          }
          print "  <TR ALIGN=CENTER><TH>$rownames[0]</TH><TD>\n"; 
          $rownames = 1;
        } else { print "  <TR ALIGN=CENTER><TD>\n"; }
	$rowcount = 0; # rowcount is 0-based


  ## D. VARIABLES

      } elsif ($item eq 'var' || $item eq 'varlong') {
        ($size,$lhs,$rhs) = ($form[$i++], $form[$i++], $form[$i++]);
        if ($bigstr) { $size=$BIG_SIZE; } #bigstr flag set in B.(documentation)
	# hardwired names that will be recognized as being files - PDA
        if ( $lhs =~ /_infile/ || $lhs =~ /_outfile/ ||
             $lhs =~ /_library/ || $lhs =~ /output_root/ ) { $size=$BIG_SIZE; }

  ## D. 1) ANY VARIABLE REQUIRING MORE THAN 1 LINE (E.G.: TEXTAREA)

        if ($item eq 'varlong') {
          @lines = split(/\n/,$rhs);
          $lines = $MIN_TEXTAREA_ROWS;
          if ($lines < ($#lines + 1)) { $lines = ($#lines + 1); }
          print "    </TD></TR>\n";	# always break row for textarea
          print "  <TR><TD $COL3><TEXTAREA NAME=$lhs ROWS=$lines ";
          print "COLS=$TEXTAREA_COLS>\n$rhs</TEXTAREA>\n";
          print "    <INPUT TYPE=BUTTON VALUE='More Lines' ";
#          ($doc_safe = $this_doc) =~ s/'/\'/g;
          $doc_safe = &escape($this_doc);
          print "onClick=\"editWindow('${lhs}','${doc_safe}',";
          print   "$AUX_TEXTAREA_ROWS,$TEXTAREA_COLS)\">\n";
          &print_cnsInfo_button($lhs);
          undef( @lines );
	  undef( $lines );
	  undef( $rhs );
	  undef( $lhs );
	  undef( $doc_safe );

  ## D. 3) CHOICE VARIABLE

        } elsif ($choice_str) {

        ## Split up choices, grouping characters enclosed in double-quotes as
        ##   a single choice, or non-whitespace strings as single choices.
        ##   Example: {+ choice: "this" "that" "" user_file +} 
        ##   becomes a list with three elements: (this,that,,),
        ##   where the third element is null and user_file indicates the 
        ##   <<other>> field should be given.

          while ($choice_str !~ /^\s*$/) {
            if ($choice_str=~/^\s*"\s*([^"]*)"(.*)$/) {
              ($choice, $choice_str) = ($1, $2);
            } elsif ($choice_str=~/^\s*(\S+)(.*)$/) {
              ($choice, $choice_str) = ($1, $2);
            } # should be no 'else'
            if ($choice =~ m/^$OTHER/) { $other=1; }
            else {
              if ($choice =~ m/^\s*$/) { $choice = $NONE; }
              push( @choices, $choice );
            }
          }

  ## D. 3a) RADIO BUTTONS

          $ch_len=0;
          foreach $choice (@choices) {
            $ch_len += length($choice);
          }
	  # if we have few choices...
          if ( $#choices < $A_FEW && ! $other && $ch_len < $RADIO_LONG ) {
            print "<TABLE><TR>\n";
            foreach $choice (@choices) {	# print radio buttons
              print "    <TD><INPUT TYPE=RADIO NAME=$lhs VALUE=\"$choice\"";
              if ( $choice eq $rhs || $rhs =~ /^\s*$/ && $choice eq $NONE )
                { print ' CHECKED'; }
              print ">$choice</TD>\n";
            }
            if ( $print_info_buttons ) {
              print "";
            }
            else {
              print "</TR></TABLE>\n";
            }

  ## D. 3b) SELECTION LIST

          } else {				# else print selection list
            print "    <SELECT NAME=$lhs>\n";
            $size = $STR_SIZE;
            $selected = 0;
            foreach $choice (@choices) {
              if ( $choice eq $rhs || $rhs =~ /^\s*$/ && $choice eq $NONE ) {
                print '    <OPTION SELECTED>';
                $selected = 1;
              }
              else { print '    <OPTION>'; }
              print "$choice</OPTION>\n";
              if (length($choice) > $size) { $size = length($choice); }
              }
            print "    </SELECT>\n";
            if ($other) {
              if ($size > $STR_SIZE || $tabcol == 1) { 
                 print "    <BR>\n"; 
                 $size = $size + $SIZE_INC;
              }
              print "    Other: ";
              if ( $selected ) {
                print "<INPUT TYPE=CHECKBOX NAME=$lhs$OTHERFLAG VALUE=1>\n";
                print "           <INPUT TYPE=TEXT NAME=\"$lhs$OTHERVAL\" ";
                print "VALUE=\"\" SIZE=$size>\n";
              }
              else {
                print "<INPUT TYPE=CHECKBOX CHECKED NAME=$lhs$OTHERFLAG VALUE=1>\n";
                print "           <INPUT TYPE=TEXT NAME=\"$lhs$OTHERVAL\" ";
                print "VALUE=\"$rhs\" SIZE=$size>\n";
              }
            }
          }
          if ($#choices < $A_FEW && !$other && $ch_len < $RADIO_LONG &&
              $print_info_buttons) {
            print "<TD>";
            &print_cnsInfo_button($lhs);
            print "</TD></TR></TABLE>\n";
          }
          else {
            &print_cnsInfo_button($lhs);
          }
          $choice_str=''; $other=0; undef( @choices ); 

  ## D. 4) NORMAL VARIABLE

        } else { &output_var($size,$lhs,$rhs); }
        while ($form[$i] eq 'break') { $i++; }
        if (!$nested_table) {
          if ($form[$i] eq 'doc') { print "    </TD></TR>\n"; }
          elsif ($form[$i] eq 'var') { print "    <BR>\n"; }
        } else {
          print "    </TD>";
          if ($form[$i] eq 'doc') { 
            foreach (1..$numcols-$colcount) { print "<TD>&nbsp;</TD>"; }
            print "</TR>\n";
            &cleanup_nested_table();
          } else { 
            $colcount++;
            if ($colcount <= $numcols) {	# colcount is 1-based
              print "<TD>\n";
            } else {
              $rowcount++;
              $colcount=1;
              print "</TR>\n";
              if ($rowcount >= $numrows) {	# rowcount is 0-based
                &cleanup_nested_table();
              } else {
                if ($rownames) {
                  $r = $rownames[$rowcount];
                  if ($r=~m/^\s*$/) { $r = '&nbsp;'; }
                  print "  <TR ALIGN=CENTER><TH>$r</TH><TD>"; 
                } else { print "  <TR ALIGN=CENTER><TD>"; }
              }
            } 
          }
        }


  ## E. OTHER

      } else {}			# NO DEFAULT FOR TABLE BLOCK (may be 'break')
    }
    print "</TABLE>\n";
  } else {}				# NO FORM PROCESSED (e.g.: first form)
}

sub cleanup_nested_table {
  print "</TABLE>\n";		# end nested table
  if ($ok_to_nest) {
    print "</BLOCKQUOTE></TD></TR>\n";
  } else {
    print "<TABLE BORDER=$BORDER WIDTH=\"100\%\">\n"; # continue original table
    print "<TR><TD>$SPACER1<TD>$SPACER1<TD>$SPACER1</TR>\n";
  }
  $nested_table=0;		# nested_table flag set in C. (table formatting)
  undef( $numrows );
  undef( $rowcount );
  undef( $rownames );
  undef( @rownames );
  undef( $numcols );
  undef( $colcount );
  undef( $colnames );
  undef( @colnames );
}

sub output_var { local( $size, $lhs, $rhs ) = @_;
  if ($lhs=~m/$QUOTEFLAG$/ || $rhs=~m/\s/)	# if quoted or contains spaces
  { $js_err = ''; $rhs = '"' . $rhs . '"'; } 
  else {$js_err="onChange=\"assertReal( this, $rhs)\"";}
  print "    <INPUT NAME=$lhs VALUE=$rhs SIZE=$size $js_err>\n"; 
  &print_cnsInfo_button($lhs);
}

sub print_cnsInfo_button { local( $lhs ) = @_;
  if ($lhs !~ /^$TEXTAREAROOT\d+$/ && $print_info_buttons) {
    if ($lhs =~ m/^(.*)$PARENFLAG$/) { $lhs = $1; }
    elsif ($lhs =~ m/^(.*)$QUOTEFLAG$/) { $lhs = $1; }
    print "      <INPUT TYPE=BUTTON VALUE='${INFOWIN_ICON}' 
      onClick=\"cnsInfo('${lhs}','${familyName}')\">\n";
  }
}


##### FORM2INP #####

sub form2inp {
  local($buf,%formval,$lhs,$rhs,@lines,@newlines,$line);
  local($inp_file,@vars,$prelen,$indent);

  $output="";

  read(STDIN, $buf, $ENV{'CONTENT_LENGTH'});  
  %formval = &parse_args( $buf );
  foreach $lhs (keys %formval) {
    if ( $lhs eq submit_view ) {
      $output .= "Content-type: text/plain\n\n";
    } elsif ( $lhs eq submit_download ) {
      $output .= "Content-type: application/octet-stream; name=$formval{$TRUE_FNAME}\n";
      $output .= "Content-disposition: attachment; filename=$formval{$TRUE_FNAME}\n\n";
    } else {
      if ($formval{"$lhs$OTHERFLAG"}) { $rhs=$formval{"$lhs$OTHERVAL"}; }
      else { $rhs = $formval{$lhs}; }
      if ($rhs eq $NONE || $rhs =~ m/^\s+$/) { $rhs = ''; }
      else {
        @lines = split(/\n/,$rhs);
        @newlines = ();
        foreach $line (@lines) {
          if ($line =~ m/^\s*(\S.*)$/) { $line = $1; } # clear beg spaces
          if ($line =~ m/^(.*\S)\s*$/) { $line = $1; } # clear end spaces
          if ($line !~ m/^\s*$/) {push(@newlines,$line);} # keep non-blank lines
        }
        $rhs = join("\n",@newlines);
      }
      $formval{$lhs} = $rhs;
    }
  }

  ($inp_file = $formval{$INP_FILE}) || &die_html("No input file given.");
  open(INP,$inp_file) || &die_html("Cannot read source file ($inp_file).");

  while ( $line=<INP> ) {
    if ($line =~ m/^\s*\{===>\}\s*(.*)$/) {
      @vars = &get_vars(INP,$1);
      $output .= '{===>}';
      foreach $var (@vars) { 
        ($lhs,$rhs) = split('.cns_var_is_set_to.',$var);
        $rhs = $formval{$lhs};
        if ($lhs =~ m/^$TEXTAREAROOT(\d+)$/) {
          $rhs = "\n$rhs";
          $rhs =~ s/\n/\n       /g;
          $output .= "$rhs\n{<===}";
        } else {
          if ($lhs =~ m/^(.*)$PARENFLAG$/) 
            { $lhs = $1; 
	      if ( $rhs !~ /^\s*\([^()]*\)\s*$/ ) {
		  $rhs = "($rhs)";
	      }
	      $prelen=9; }
          elsif ($lhs =~ m/^(.*)$QUOTEFLAG$/) 
            { $lhs = $1;
	      if ( $rhs !~ /^\s*\"[^"]*\"\s*$/ ) {
		  $rhs = "\"$rhs\"";
	      }
	      $prelen=9; }
          else { 
            if ( $rhs =~ m/^\s*$/ ) {
              &die_html("Input field <I>$lhs</I> is blank, this must have a value.<BR>
                         <FORM><INPUT TYPE=\"button\" 
                                      VALUE=\"Return to editing file\"
                                      onClick=\"history.back()\"></FORM>");
            } 
            else {
              $prelen=8;
            }
          }
          $indent = ' ' x ($prelen + length($lhs));
          $rhs =~ s/\n/\n$indent/g;
          $output .= " ${lhs}=$rhs;";
        }
      }
      $output .= "\n";
    } 
    else { 
	$output .= "$line";
    }
  }
  print $output;
  close INP;
}

# EOF --> <H3>An error has occurred.</H3>
