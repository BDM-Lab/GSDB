#! /usr/bin/perl -w

print "<html>\n";
print "<body>\n";

while (<>) {
  next if (/^#/);
  if (/vars:\s*(.*)/) {
    $vars = $1;
    $vars =~ s/\&/\&amp;/g;
    $vars =~ s/</\&lt;/g;
    $vars =~ s/>/\&gt;/g;
    $vars =~ s/"/\&quot;/g;
    print "<hr>\n";
    print "<strong>$vars</strong>\n";
    print "<p>\n";
  }
  elsif (/info:\s*(.*)/) {
    print "$1\n";
  }
  else {
    print;
  }
}

print "<hr>\n";
print "</body>\n";
print "</html>\n";

__END__
