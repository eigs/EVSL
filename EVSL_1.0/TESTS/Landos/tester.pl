#!/usr/bin/perl -lw
use strict;

my $regex = "RMSE"; #Regex to match
my %data = ( 
  'RMSE' => 0.3, 
  'RMSE1' => 0.3, 
  'RMSE2' => 0.3
); 
my $flag = 0;

while (my $line = <>) {
  chomp($line);
  my @x = split(/ /,$line);
  if($x[1] =~ $regex) {
    if ($x[2] < $data{$x[1]}) {
      print $line .  " Passed"; 
    }
    else {
      print $line .  "Failed";
      $flag = 1
    }
  }
}
exit $flag;
