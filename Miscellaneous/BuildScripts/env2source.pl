#!/usr/bin/perl
# A small wrapper around calls like
#   ./env2source.pl FOO=BAR baz2=123
# which outputs the variables suitable for sourcing in a shellscript.
# Use this to convert command line arguments to shell/environment variables, ie.
#   VARS=$(./env2source.pl $@) && eval $VARS || { echo "Failure parsing $@"; exit -1; }
#
# You can extend/control the behaviour with a single meta argument --mapper=reducechar
# having mapper one of the quantitiy mapping functions (as indicated below) and the
# reducing (joining) character by default a newline.
# Useful examples are: --keys=,  or --values=\;
#
# Output keys are always sorted by name (as the internal order is lost in the hash).
# Sven K, 2017-06-03

use strict;
my $ret=0, my %vals, my %map, my %reduce;
my $mapper="source", my $reducechar="NL"; # defaults, can be changed with argument

# available mappers:
$map{'source'} =  sub { return "$_[0]=\"$_[1]\""; }; # code which can be sourced by shells
$map{'export'} = sub { return "export " . $map{'source'}(@_); }; # exporting variables
$map{'keys'} = sub { return $_[0]; }; # only the keys
$map{'values'} = sub { return $_[1]; }; # only the values

foreach (@ARGV) {
  # control argument
  if(m/^--([a-z]+)(?:=([^\s+k]+))?$/i) { $mapper=$1; $reducechar=$2 if($2);
    if(not exists $map{$mapper}) { print STDERR "Mapper $mapper does not exist.\n"; exit 2; }
    next; }
  # data argument
  if(m/^([A-Za-z][A-Za-z0-9-_]*)=([^\s]+)\s*$/) 	{ $vals{$1} = $2; }
  else {print STDERR "Malformed command line argument: $_\n"; $ret=1;}
}

# service: convert control characters
$reducechar = "\n" if("nl" eq lc $reducechar);
$reducechar = "\0" if("zero" eq lc $reducechar); # zero byte (like in find -1)

#use Data::Dumper; print Dumper(%vals);

my @lines = map { $map{$mapper}($_, $vals{$_}); } sort keys %vals;
print join($reducechar, @lines);

# for control characters, print them also at the end
print $reducechar if($reducechar eq "\n" or $reducechar eq "\0");

exit $ret;
