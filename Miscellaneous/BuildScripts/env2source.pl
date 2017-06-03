#!/usr/bin/perl
# A small wrapper around calls like
#   ./env2source.pl FOO=BAR baz2=123
# which outputs the variables suitable for sourcing in a shellscript.
# Use this to convert command line arguments to shell/environment variables, ie.
#   VARS=$(./env2source.pl $@) && eval $VARS || { echo "Failure parsing $@"; exit -1; }
# Sven K, 2017-06-03

my $ret=0;
foreach (@ARGV) {
  if(m/([A-Za-z][A-Za-z0-9-_]*)=([^\s]+)\s*/)
  {print "$1=\"$2\"\n";}
  else {print STDERR "Malformed command line argument: $_\n"; $ret=1;}
}
exit $ret;
