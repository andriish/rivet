#!/usr/bin/perl -pi

if ( /\@LIBTOOL\@/ ) {
  print "ifdef SHOWCOMMAND\n";
  print;
  print "else\n";
  s/\@LIBTOOL\@/\@LIBTOOL\@ --quiet/;
  print;
  print "endif\n";
  $_ = "";
}
elsif ( /^\s+\$\(CXXLINK\)/ ) {
  print "ifdef SHOWCOMMAND\n";
  print;
  print "else\n\t\@echo \"linking \$(subdir)/\$\@ ...\"\n";
  s/\$\(CXXLINK\)/\@\$\(CXXLINK\)/;
  print;
  print "endif\n";
  $_ = "";
}
elsif ( /^\s+\$\(LINK\)/ ) {
  print "ifdef SHOWCOMMAND\n";
  print;
  print "else\n\t\@echo \"linking \$(subdir)/\$\@ ...\"\n";
  s/\$\(LINK\)/\@\$\(LINK\)/;
  print;
  print "endif\n";
  $_ = "";
}
elsif ( /^\s+\$\(F77COMPILE\)/ ) {
  print "ifdef SHOWCOMMAND\n";
  print;
  print "else\n\t\@echo \"compiling \$(subdir)/\$\@ ...\"\n";
  s/\$\(F77COMPILE\)/\@\$\(F77COMPILE\)/;
  print;
  print "endif\n";
  $_ = "";
}
elsif ( /^\s+\$\(LTF77COMPILE\)/ ) {
  print "ifdef SHOWCOMMAND\n";
  print;
  print "else\n\t\@echo \"compiling \$(subdir)/\$\@ ...\"\n";
  s/\$\(LTF77COMPILE\)/\@\$\(LTF77COMPILE\)/;
  print;
  print "endif\n";
  $_ = "";
}
elsif ( /\.cc\.lo/ || /\.cc\.o/ || /\.cc\.obj/ ) {
  $chunk = "else\n\t\@echo \"compiling \$(subdir)/\$\@ ...\"\n";
  print;
  print "ifdef SHOWCOMMAND\n";
  $_ = "";
}
elsif ( $_ =~ /^\s*$/ && $chunk ) {
  $chunk .= "endif\n";
  $chunk =~ s/if \$\(LTCXXCOMPILE\)/\@if \$\(LTCXXCOMPILE\)/g;
  $chunk =~ s/if \$\(CXXCOMPILE\)/\@if \$\(CXXCOMPILE\)/g;
  print $chunk;
  $chunk = "";
}
elsif ( $chunk ) {
  $chunk .= $_;
}
