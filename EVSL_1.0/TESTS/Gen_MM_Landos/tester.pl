#!/usr/bin/perl -lw
use strict;

my $regex = "number of eigenvalues found"; #Regex to match
my %data = ( 
  'RMSE1' => 0.3, 
  'RMSE2' => 0.3
); 
my $nev = 0;

while (my $line = <>) {
  chomp($line);
  my @x = split(/ /,$line);
  if($line =~ $regex) {
    $nev += $x[5]
  }
}
if ($nev == 707) {
  print('./MMRLan eigenvalue count Passed');
  exit 0;
}
else {
  print('./MMRLan eigenvalue count Failed');
  exit -1;
}

