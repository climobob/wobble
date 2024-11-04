#!/usr/bin/perl

#Scan orbit (Horizons) file for julian day and earth-sun distance, plus rdot

while ($line=<STDIN>) {
  chop($line);
  @words = split(/,/,$line);
  #$words = split($line,/,/);
  #if (($words[0] - int($words[0]) ) == 0) {
    print $words[1] - 2400000,"	",$words[8],"	",$words[9],"\n";
  #}
}
# 1940-Jan-01 00:00, 2429629.500000000, , , 18 44 38.80,-23 03 00.8,  0.000000000000,  0.0000000, .983275809925875, -0.0143278,       n.a.,      n.a.,

