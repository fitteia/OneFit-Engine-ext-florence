#!/usr/bin/env perl
# Are the reference fit's errors right? Refits N data sets that differ only
# in their noise (make_data.c SEED, seeds 1..N) and prints, per parameter,
# the mean and spread of the pulls (fitted - true) / error: about 0 and 1
# when the errors are right (within about 0.14 and 0.10 for N = 50), and
# the fraction within 1 sigma (about 68%). Not part of the test: about 3 s
# per fit; the README records a run.
#
#   OFE_PATH=/opt/onefite-go ONEFITE=/opt/onefite-go/bin/onefite perl pulls.pl [N]
use strict;
use warnings;
use File::Temp qw(tempdir);
use File::Basename qw(dirname);
use Cwd qw(abs_path);
my $here;
BEGIN { $here = dirname(abs_path($0)) }
use lib $here;
use FitRef qw(%TRUE model in_path find_engine run_fit);

my $n = shift // 50;
my $root = $ENV{OFE_PATH} or die "set OFE_PATH to the engine's install root\n";
my $engine = find_engine() or die "no onefite: set ONEFITE\n";
my $cc = in_path('cc') // in_path('gcc') or die "no C compiler\n";
my $work = tempdir(CLEANUP => 1);
system($cc, "-I$root/include/ext/florence", "$here/make_data.c", "-L$root/lib",
    '-lonefit-ext-florence', '-lgfortran', '-lm', '-o', "$work/make_data") == 0
    or die "make_data.c did not build\n";

my (%pull, @bad);
for my $seed (1 .. $n) {
    my $dir = "$work/$seed";
    mkdir $dir or die $!;
    system("\"$work/make_data\" $seed > \"$dir/data.dat\"") == 0 or die "make_data $seed failed\n";
    my $r = run_fit($engine, model(2e-11, 0.03, 0.6), $dir, 'data.dat');
    if ($r->{rc} || ($r->{migrad} // '') ne 'CONVERGED' || !$r->{accurate}
        || ($r->{minos} // '') ne 'SUCCESSFUL' || keys %{ $r->{par} } != 3) {
        push @bad, $seed;
        next;
    }
    push @{ $pull{$_} }, ($r->{par}{$_}[0] - $TRUE{$_}) / $r->{par}{$_}[1] for keys %TRUE;
    printf STDERR "seed %d: chi2 %s\n", $seed, $r->{chi2};
}

printf "%d fits, %d not converged cleanly%s\n", $n, scalar @bad, @bad ? " (seeds @bad)" : '';
for my $p (qw(TR D2 S)) {
    my @x = @{ $pull{$p} // [] } or next;
    my $mean = 0;
    $mean += $_ / @x for @x;
    my $var = 0;
    $var += ($_ - $mean)**2 / (@x - 1) for @x;
    my $in1 = grep { abs($_) < 1 } @x;
    printf "%-3s pull mean %+.3f  spread %.3f  within 1 sigma %d%%\n", $p, $mean, sqrt($var), 100 * $in1 / @x;
}
