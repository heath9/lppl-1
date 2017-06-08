#!/usr/bin/perl -w
use strict;

use IO::Handle;
use constant PI => 4 * atan2(1, 1);

my ($n, $B, $T, $m, $C, $w, $p, $variance, $tau, @seed, $intensity, $timefile, 
    $freqfile, $i, $B_base, $T_base, $m_base, $C_base, $w_base, $p_base,
    $n1, $alpha);

sub safe_unlink {
  my $filename;

  $filename = shift;
  if (-e $filename) {
    unlink $filename;
  }
}


sub genpar {
  my ($filename, $m, $C);

  $B = shift;
  $T = shift;
  $m = shift;
  $C = shift;
  $w = shift;
  $p = shift;

  $filename = "base.par";
  safe_unlink($filename);
  open(PAR, "> $filename");
  print PAR "A=100\nB=$B\nT=$T\nm=$m\nC=$C\nw=$w\np=$p\n";
  close PAR;

  return $filename;
}


sub launch {
  my ($x, $i, $filename, $basename, $parname, $intensity, $m, $C);

  $parname = genpar @_;

  $B = shift;
  $T = shift;
  $m = shift;
  $C = shift;
  $w = shift;
  $p = shift;
  $variance = shift;
  $tau = shift;
  $basename = shift;

  print "$basename ...\n";
  for ($i = 0; $i <= $#seed; $i++) {
    $x = $seed[$i];
    $filename = $basename . "t" . $i;
    $timefile = $filename . ".plt";
    $freqfile = $filename . "-dft.plt";
    safe_unlink($filename);
    safe_unlink("tmp0.plt");
    safe_unlink("tmp1.plt");
    system ("gen-lppl -l --nopenny -p $parname -N 0 --diffuse=$variance --relax=$tau --seed=$x | cut -f 2 | head -25000 > tmp0.plt");
    system ("tac tmp0.plt >  tmp1.plt");
    system ("tail -n +2 tmp0.plt | head -n -1 >> tmp1.plt");
    safe_unlink("tmp0.plt");
    system ("nl -n ln tmp1.plt > $timefile");
    safe_unlink("tmp1.plt");
    system ("cut -f 2 $timefile | fft -n 49998 > $freqfile");
  }
}


sub gnuplotto {
  my ($basename, $filename, $gptname); 

  $basename = shift;

  $gptname = $basename . ".gpt";
  unlink $gptname;
  open (PLOT, "> $gptname");

  print PLOT "set logscale x\n";
  print PLOT "set logscale y\n";
  print PLOT "set grid x\n";
  print PLOT "set grid y\n";

  print PLOT "set yrange[1e-10:]\n";
  print PLOT "set xrange [:.5]\n";

  print PLOT "set key top right\n";

  print PLOT "n2=49998\n";

  print PLOT "plot ";
  $filename = $basename . "t0-dft.plt";
  print PLOT "\"$filename\" using (\$1/n2):(\$2*sqrt(n2)) title \"$basename\" with points ," ;
  $filename = $basename . "-mt0-dft.plt";
  print PLOT "\"$filename\" using (\$1/n2):(\$2*sqrt(n2)) title \"$basename C=0\" with points ," ;
  $filename = $basename . "-Ct0-dft.plt";
  print PLOT "\"$filename\" using (\$1/n2):(\$2*sqrt(n2)) title \"$basename m=0\" with points\n" ;
  print PLOT "pause -1 \"Hit Enter to continue\"";
  close(PLOT);

  system "gnuplot $gptname";

  print "Hit Enter again to continue, q to exit\n";
  $_ = <>;
  chomp;
  exit if ($_ eq "q");
}


# For each case, compare with: 
# - C=0 (pure power law)
# - m=0 (pure log-periodic oscillations)
sub trilaunch {
  my ($basename, $m, $C);

  $B = shift;
  $T = shift;
  $m = shift;
  $C = shift;
  $w = shift;
  $p = shift;
  $variance = shift;
  $tau = shift;
  $basename = shift;

  launch($B, $T, $m, $C, $w, $p, $variance, $tau, $basename);
  launch($B, $T, 0., $C, $w, $p, $variance, $tau, $basename . "-C");
  launch($B, $T, $m, 0., $w, $p, $variance, $tau, $basename . "-m");

  gnuplotto($basename);
}


# @seed = (663411404, 44046727, 672331627, 363612029, 329747984);
@seed = (663411404);

$B = $B_base = 0.01;
$T = $T_base = 26000;
$m = $m_base = 0.7;
$C = $C_base = 0.05;
$w = $w_base = 2 * PI;
$p = $p_base = 0.;
$variance = 0.;
$tau = 1.;

trilaunch($B, $T, $m, $C, $w, $p, $variance, $tau, "base");

for ($m = .1; $m <= 1.; $m += 0.1) {
  trilaunch($B, $T, $m, $C, $w, $p, $variance, $tau, "m" . $m);
}
$m = $m_base;

for ($T = 26500; $T <= 28000; $T += 500) {
  trilaunch($B, $T, $m, $C, $w, $p, $variance, $tau, "T" . $T);
}
$T = $T_base;

for ($C = 0.01; $C <= 0.04; $C += .01) {
  trilaunch($B, $T, $m, $C, $w, $p, $variance, $tau, "C" . $C);
}
$C = .1;
trilaunch($B, $T, $m, $C, $w, $p, $variance, $tau, "C" . $C);
$C = .5;
trilaunch($B, $T, $m, $C, $w, $p, $variance, $tau, "C" . $C);
$C = $C_base;

for ($w = 7; $w <= 14; $w += 1) {
  trilaunch($B, $T, $m, $C, $w, $p, $variance, $tau, "w" . $w);
}
$w = $w_base;

for ($p = .6; $p <= PI; $p += .6) {
  trilaunch($B, $T, $m, $C, $w, $p, $variance, $tau, "p" . $p);
}
$p = $p_base;
