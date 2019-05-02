#!/usr/bin/perl
# <!--
$this_cgi = 'cns_save.cgi';
$modified = 'Last modified 12/08/98 Paul Adams';
#
# CGI file for saving CNSsolve files
#
# Written by: Jesse Reklaw and Paul Adams
#
# copyright Yale University
#

require 'cns.lib';
 
##### PROGRAM FLOW OF CONTROL #####

$argstr = &init_cgi();

# For backwards compatibility, if there is no equals sign ('=')
#   in the $argstr, then the $argstr is the $inp_file.
#   NOTE: the assumption is that no filename can have '=' in it.
#   KLUDGE! for security, some sanity checking should be done on filenames.
if ($argstr !~ m/\=/) {
    $inp_file = $argstr;
}
else {
    %argnv = &parse_args($argstr);
    $inp_file = $argnv{'file'};
    if (!$inp_file) { $inp_file = $argnv{$INP_FILE}; }
}

if ( $inp_file !~ /^\/\w/ && $inp_file) { # shortcut for no full path
    $inp_file = "$www_root/$cns_version/" . $inp_file; 
}
else {
    &die_html( "cannot open file $inp_file" );
}

undef( %argnv );
undef( $argstr );

open (IN_FILE, "< $inp_file") || &die_html( "could not open file $inp_file" );

print "Content-type: application/octet-stream; name=$inp_file\n";
print "Content-disposition: attachment; filename=$inp_file\n\n";

while ( <IN_FILE> ) {
    print "$_";
}

close (IN_FILE);

# EOF --> <H3>An error has occurred.</H3>
