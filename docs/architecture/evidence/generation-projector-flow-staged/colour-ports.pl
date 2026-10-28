use strict;
use warnings;
use JSON::PP;
local $/;
open my $input, '<', (shift // 'tests/resources/graphs/qqx_aaa_subtracted.dot') or die $!;
my $source = <$input>;
my @certificates;
for my $graph (split /(?=digraph\s)/, $source) {
    next unless $graph =~ /^digraph\s+(\S+)/;
    my $name = $1;
    my ($projector) = $graph =~ /^projector="([^"]*)"/m;
    next unless defined $projector && $projector =~ /spenso::cof/;
    my %ports;
    while ($graph =~ /(?:^projector|\bnum)="([^"]*)"/mg) {
        my $factor = $1;
        $factor =~ s/gammalooprs::hedge\((\d+),1\)/hedge_$1/g;
        $factor =~ s{spenso::dind\(spenso::cof\(3,(hedge_\d+)\)\)}{
            $ports{$1}{dual}++; 'DUAL'
        }ge;
        while ($factor =~ /spenso::cof\(3,(hedge_\d+)\)/g) {
            $ports{$1}{fundamental}++;
        }
    }
    for my $port (sort keys %ports) {
        die "$name $port does not have exactly one fundamental and one dual\n"
            unless ($ports{$port}{fundamental} // 0) == 1 && ($ports{$port}{dual} // 0) == 1;
    }
    die "$name missing physical external ports\n" unless exists $ports{hedge_1} && exists $ports{hedge_2};
    push @certificates, { graph => $name, exact_incidence => \%ports };
}
die 'expected eight explicit colour projectors' unless @certificates == 8;
print JSON::PP->new->canonical->pretty->encode({
    certificate => 'exact fundamental/dual port incidence; no numerator expansion',
    scope => 'exact authored factor incidence in native model-one slots; runtime admission checked separately',
    graphs => \@certificates,
});
