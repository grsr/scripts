# find proxy SNPs for a set of lead SNPs

use strict;
use warnings;

use Data::Dumper;

my $INCLUDE_POS = 1;

my $tag_dir = "/lustre/scratch113/teams/soranzo/users/vi1/uk10k/clump_distance/ALL";

my ($snp_file) = @ARGV;

die "Usage: $0 <snps>\n" unless $snp_file;

open my $SNPS, "<$snp_file";

my %lead_snps;

while (<$SNPS>) {
    chomp;
    my ($chr, $start, $end, $id, $pval) = split /\t/;
    $lead_snps{$id} = [$chr, $start, $end, $pval];
}

exit if keys %lead_snps == 0;

close $SNPS;

# get tags

my %tags;

for my $chr (1..22, 'X') {

    my $TAGS;

    my $tf = "${tag_dir}/TAGS.r02.chr${chr}.txt";

    if (-e $tf) {
        open $TAGS, "<$tf";
    }
    elsif (-e $tf.'.gz') {
        $tf .= '.gz';
        open $TAGS, "gunzip -c $tf |";
    }
    else {
        die "Can't find a tag file for chr $chr\n";
    }


    while (<$TAGS>) {
        # SNP  CHR         BP NTAG       LEFT      RIGHT   KBSPAN TAGS
        # rs181691356   21    9411245    0    9411245    9411245        0 NONE
        # chr21:9411318   21    9411318    2    9411318    9575754  164.436 chr21:9412269|rs140668083
        s/^\s+//;
        next if /^SNP/;
        my @cols = split /\s+/; 
        my $id = $cols[0];
        if ($lead_snps{$id}) {
            $tags{$id} = $cols[7];
        }
    }

    close $TAGS;
}

for my $lead (keys %lead_snps) {
    my @tags = ($lead);
    if ($tags{$lead} && $tags{$lead} ne 'NONE') {
        push @tags, split /\|/, $tags{$lead};
    }
    print join("\n", @tags), "\n";
}

