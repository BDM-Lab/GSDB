#!/usr/bin/perl

use lib "/var/www/html/3dgenome/GSDB/MIME-Lite-2.117/lib/";

use MIME::Lite;

my $from = "MULTICOM-CLUSTER <MULTICOM-CLUSTER\@missouri.edu>";
my $to = "oeow39@mail.missouri.edu";
my $smtp = "localhost";


$date = `date`; 
chomp $date;
$email_content="[$date]: Remind of RR category for novel: Start running on lotus server under !!!!\n\n";
$email_content .="[$date]: Start running <perl /storage/casp13/DNCON2.0/RR_targets/scripts/multicom-novel.pl >\n\n";


	$msg1 = MIME::Lite->new(
			From    => "$from",
			To      => "$to",
			Subject => "Remind of RR category for novel: Start running on lotus server!!!!",
			#Cc     => 'tianqiwu@mail.missouri.edu',
			#Cc     => 'aomqc@mail.missouri.edu',
			Type    => 'multipart/mixed'
	) or die ("cannot create MIME object!");

	$msg1->attach(
		Type     =>'TEXT',
		Data     =>"$email_content"
	); 

	$msg1->send or print "Error sending email, MIME::Lite->send failed: $!\n";



