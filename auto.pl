#!/usr/bin/perl -w

use strict;
use Getopt::Std;

my $doquiet = "./quiet.pl";

my %opts;

getopts('vcCt:', \%opts);

my $quiet = !defined $opts{'v'};
my $doconf = defined $opts{'c'};
my $docheck = defined $opts{'C'};
if ( $docheck ) {
  $doconf = 1;
}

my $tee = "";
if ( defined $opts{'t'} ) {
  system("echo \"make check and distcheck\" > $opts{'t'}") == 0 or die "aborted";
  $tee = "2>&1 | tee -a " . $opts{'t'} if defined $opts{'t'};
}

system("aclocal; autoconf; libtoolize --automake; autoheader; automake  --foreign --add-missing") == 0 or die "aborted";
if ( $quiet ) {
  my $files = `find . -name Makefile.in`;
  $files =~ s/\n/ /gs;
  system("$doquiet $files") == 0 or die "aborted";
}
if ( $doconf ) {
  system("./configure $tee") == 0 or die "aborted";
}
if ( $docheck ) {
  system("(time make check) $tee") == 0 or die "aborted";
  system("(time make distcheck) $tee") == 0 or die "aborted";
}
