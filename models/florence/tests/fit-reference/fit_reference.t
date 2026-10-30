#!/usr/bin/env perl
# A converged Florence fit with known answers: FlorenceLS-R1.dat is R1 from
# FlorenceLS with TAURM = 2.88e-11 s, DELTA2 = 0.0236 (x 1e20 s^-2) and
# SLS = 0.7 (TAUVM = 2.86e-12 s fixed; see make_data.c), plus a fixed 1%
# Gaussian noise sequence. Fitting TAURM, DELTA2 and SLS, once from below
# and once from above the truth (20-30% off), must each end with MIGRAD
# converged, an accurate error matrix and MINOS successful, at the same
# minimum; every value must lie within 3 sigma of the truth, its MINOS
# errors close to the parabolic ones, and values, errors and chi2 must
# match the reference below (OFE-Go 0.5.5 and OFE-Raku 1.1.0 on .113,
# Florence 0.3.4: identical). TAUVM is fixed because R1 here depends on
# DELTA2 * TAUVM only (correlation -1.000 when both are free). That the
# errors are right, not only consistent, is pulls.pl's question (README).
#
# Needs an installed OneFit Engine with this extension:
#   ONEFITE=/opt/onefite-go/bin/onefite prove -v fit_reference.t
# (default: onefite-go, else onefite, from PATH). Skipped without one.
# OFE_PATH=<install root> also rebuilds the data with make_data.c and
# compares it with the file.
use strict;
use warnings;
use Test::More;
use File::Temp qw(tempdir);
use File::Basename qw(dirname);
use File::Copy qw(copy);
use Cwd qw(abs_path);
my $here;
BEGIN { $here = dirname(abs_path($0)) }
use lib $here;
use FitRef qw(%TRUE model in_path find_engine run_fit);

my $data = 'FlorenceLS-R1.dat';
# name => [value, parabolic error]; MINOS from the reference run
my %ref = (
    TR => [0.28804E-10, 0.50520E-12],
    D2 => [0.23611E-01, 0.65299E-04],
    S  => [0.70143,     0.67653E-02],
);
my $ref_chi2 = 38.33097;
my %starts = (
    'from below' => [2e-11, 0.03,  0.6],
    'from above' => [3.7e-11, 0.018, 0.84],
);

my $engine = find_engine();
plan skip_all => 'no onefite (set ONEFITE=/path/to/onefite)' unless $engine && -x $engine;
diag "engine: $engine";

for my $start (sort keys %starts) {
    my $work = tempdir(CLEANUP => 1);
    copy("$here/$data", "$work/$data") or die "$data: $!";
    my $r = run_fit($engine, model(@{ $starts{$start} }), $work, $data);
    is($r->{rc}, 0, "$start: the fit runs") or diag $r->{output};
    ok($r->{log}, "$start: it writes fit1.log") or next;
    is($r->{migrad}, 'CONVERGED', "$start: MIGRAD converges");
    ok($r->{accurate}, "$start: ...with an accurate error matrix") or diag "MINUIT: " . ($r->{matrix} // 'none');
    is($r->{minos}, 'SUCCESSFUL', "$start: MINOS succeeds");
    my $chi2 = $r->{chi2};
    ok(defined $chi2 && abs($chi2 - $ref_chi2) < 1e-3, "$start: chi2 " . ($chi2 // 'none') . " is the reference $ref_chi2");
    ok(defined $chi2 && $chi2 / (31 - 3) < 2, "$start: chi2 per degree of freedom below 2");
    for my $p (qw(TR D2 S)) {
        my ($v, $err, $neg, $pos) = @{ $r->{par}{$p} // [] };
        ok(defined $pos, "$start: $p in the MINOS table") or next;
        ok(abs($v - $TRUE{$p}) < 3 * $err, sprintf('%s: %s = %g +- %g, within 3 sigma of the true %g', $start, $p, $v, $err, $TRUE{$p}));
        ok(abs(-$neg / $err - 1) < 0.2 && abs($pos / $err - 1) < 0.2, "$start: $p MINOS errors within 20% of the parabolic one");
        my ($rv, $re) = @{ $ref{$p} };
        ok(abs($v - $rv) < 0.05 * $re && abs($err / $re - 1) < 0.01, "$start: $p value and error are the reference ones");
    }
}

SKIP: {
    my $root = $ENV{OFE_PATH};
    my $lib = $root && "$root/lib/libonefit-ext-florence.a";
    skip 'OFE_PATH not set, or no Florence library in it: data not rebuilt', 1 unless $lib && -f $lib;
    my $cc = in_path('cc') // in_path('gcc');
    skip 'no C compiler: data not rebuilt', 1 unless $cc;
    my $work = tempdir(CLEANUP => 1);
    my $built = system($cc, "-I$root/include/ext/florence", "$here/make_data.c", "-L$root/lib",
        '-lonefit-ext-florence', '-lgfortran', '-lm', '-o', "$work/make_data") == 0;
    skip 'make_data.c did not build here', 1 unless $built;
    my @new = map { [split] } split /\n/, qx{"$work/make_data"};
    my @old = map { [split] } split /\n/, do { local $/; open my $fh, '<', "$here/$data" or die $!; <$fh> };
    my $same = @new == @old;
    for my $i (0 .. $#old) {
        for my $j (0 .. 2) {
            $same &&= abs($new[$i][$j] / $old[$i][$j] - 1) < 1e-5;
        }
    }
    ok($same, "make_data.c rebuilds $data");
}

done_testing();
