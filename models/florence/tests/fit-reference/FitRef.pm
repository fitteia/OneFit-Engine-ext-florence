package FitRef;
# The reference fit of fit_reference.t and pulls.pl: the model, running it
# with an installed onefite, and reading MINUIT's results from fit1.log.
use strict;
use warnings;
use File::Find;
use Exporter qw(import);
our @EXPORT_OK = qw(%TRUE model in_path find_engine run_fit);

# make_data.c's parameters
our %TRUE = (TR => 2.88e-11, D2 => 0.0236, S => 0.7);

# FlorenceLS with TAURM (TR), DELTA2 (D2) and SLS (S) free from the given
# starting values, TAUVM fixed at its true value
sub model {
    my ($tr, $d2, $s) = @_;
    return "R1(f[1e4<1e9],TR:$tr\[1e-12<1e-9],D2:$d2\[1e-4<1],S:$s\[0<1])="
        . 'FlorenceLS(0,f,0.5,2.675222e8,3.5,D2,TR,2.86e-12,1e-7,2,0,0,2.003,2.003,2.003,'
        . '0,0,0,0.001,3.1,0,0,0.748,12,0,0,S,5e-12,1)';
}

sub in_path {
    my ($p) = @_;
    for my $d (split /:/, $ENV{PATH} // '') {
        return "$d/$p" if -x "$d/$p";
    }
    return;
}

# ONEFITE, else onefite-go, else onefite from PATH
sub find_engine { return $ENV{ONEFITE} // in_path('onefite-go') // in_path('onefite') }

# run_fit(ENGINE, MODEL, DIR, DATA): fits DIR/DATA in DIR/wf. Returns
# { rc, output, log, migrad, matrix, accurate, minos, chi2, par => { NAME =>
# [value, parabolic error, MINOS negative, MINOS positive] } } - the last
# MIGRAD and the MINOS of fit1.log. matrix is what MINUIT printed for the
# error matrix; accurate says whether it is: MINUIT prints either "ERROR
# MATRIX ACCURATE" or, when it has an estimate, "ERROR MATRIX UNCERTAINTY=
# x%", which it uses for an accurate and an approximate matrix alike
# (mnprin.F) - then the matrix is accurate when the parameter table MINUIT
# prints after that MIGRAD does not head its error column APPROXIMATE (or
# CURRENT GUESS; mnprin.F, ISW(2) < 3) and MINUIT added no "ERR MATRIX
# APPROXIMATE"/"NOT POS-DEF" (mnmatu.F). The MINOS table says PARABOLIC
# either way.
sub run_fit {
    my ($engine, $model, $dir, $data) = @_;
    my $out = qx{cd "$dir" && "$engine" fit '$model' "$data" --work-folder="$dir/wf" 2>&1};
    my %r = (rc => $?, output => $out, par => {});
    find(sub { $r{log} //= $File::Find::name if $_ eq 'fit1.log' }, "$dir/wf") if -d "$dir/wf";
    return \%r unless $r{log};
    my $text = do { local $/; open my $fh, '<:raw', $r{log} or die "$r{log}: $!"; <$fh> };
    my @migrad = $text =~ /FROM MIGRAD\s+STATUS=(\S+).*\n.*?(ERROR MATRIX \S+|ERR MATRIX \S+|NO ERROR MATRIX)/g;
    @r{qw(migrad matrix)} = @migrad[-2, -1] if @migrad;
    my ($after) = $text =~ /.*FROM MIGRAD\s+STATUS=\S+(.*)/s;
    $r{accurate} = defined $r{matrix} && ($r{matrix} eq 'ERROR MATRIX ACCURATE'
        || ($r{matrix} eq 'ERROR MATRIX UNCERTAINTY='
            && $after !~ /ERR MATRIX (?:APPROXIMATE|NOT POS-DEF)/
            && ($after =~ /^\s*EXT PARAMETER(.*)$/m)[0] !~ /APPROXIMATE|CURRENT GUESS/)) ? 1 : 0;
    @r{qw(chi2 minos)} = $text =~ /FCN=\s*(\S+)\s+FROM MINOS\s+STATUS=(\S+)/;
    my ($table) = $text =~ /MINOS ERRORS.*?\n(.*?)\*{10}/s;
    for (split /\n/, $table // '') {
        $r{par}{$1} = [$2, $3, $4, $5] if /^\s*\d+\s+(TR|D2|S)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
    }
    return \%r;
}

1;
