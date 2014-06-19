use strict;
use warnings;

use Data::Dumper;

my $INCLUDE_POS = 1;

#my $pval_dir = "/lustre/scratch113/projects/uk10k/users/jh21/scratch3/Meta-Fixed";
my $pval_dir = "/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/Replication/MetaAnalysis/4-way_20131008/options-qt-m-gc-gco/";
my $tag_dir = "/lustre/scratch113/teams/soranzo/users/vi1/uk10k/clump_distance/ALL";

my ($trait, $chr) = @ARGV;

die "Usage: $0 <trait> <chr>\n" unless $trait && $chr;

# get p values

open my $PVAL, "gunzip -c ${pval_dir}/${trait}.out.gz |" 
    or die "Can't open association file for $trait\n";

my %pvals;
my %extra;

while (<$PVAL>) {
    next if /^chromosome/;
    chomp;
    # chromosome      position        rs_number       reference_allele        other_allele    eaf     beta    se      beta_95L        beta_95U        z       p-value _-log10_p-value q_statistic     q_p-value       i2      n_studies       n_samples     effects
    # 1       28590   chr1:28590      TTGG    T       0.899517        0.072994        0.063377        -0.051226       0.197213        1.151732        0.253235        0.596476        1.748555        0.626193        0.000000        4       9968 +-++
    my @cols = split /\s+/;
    my $chrom = $cols[0];
    my $pos = $cols[1];
    my $id = $cols[2];
    my $eaf = $cols[5];
    my $beta = $cols[6];
    my $se = $cols[7];
    my $pval = $cols[11];

    # filter variants with MAF < 1%
    my $maf = $eaf < 0.5 ? $eaf : 1 - $eaf;
    next if $maf < 0.01;

    # convert chr 23 to X
    $chrom = $chrom eq 23 ? 'X' : $chrom;

    # we're only processing a single chromosome
    next unless $chrom eq $chr;

    $pvals{$id} = $pval;
    $extra{$id} = [$pos, $beta, $se];
}

close $PVAL;

warn "Read in ",scalar(keys %pvals), " p values\n";

# get tags

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
    die "Can't find a tag file for $trait chr $chr\n";
}

my %tags;

while (<$TAGS>) {
    # SNP  CHR         BP NTAG       LEFT      RIGHT   KBSPAN TAGS
    # rs181691356   21    9411245    0    9411245    9411245        0 NONE
    # chr21:9411318   21    9411318    2    9411318    9575754  164.436 chr21:9412269|rs140668083
    s/^\s+//;
    next if /^SNP/;
    my @cols = split /\s+/; 
    my $id = $cols[0];
    # don't worry about NONE values, as no SNP will have this as an ID
    $tags{$id} = $cols[7];
}

close $TAGS;

warn "Read in tags for ",scalar(keys %tags), " variants\n";

# find independent SNPs

my %seen;

# loop over the SNPs with the strongest association first

my $ind = 0;

#print join("\t", qw(ID CHR POS BETA SE P)), "\n";

for my $snp (sort {$pvals{$a} <=> $pvals{$b}} keys %pvals) {
    
    # skip this SNP if we've already seen it or one of its proxies
    next if $seen{$snp};

    # also skip if it's not in the tags file (as we can't be sure it is independent)
    next unless $tags{$snp};

    # if we get here the SNP is independent, so print the ID and the p value
    my ($pos, $beta, $se) = @{ $extra{ $snp } };
    print join("\t", $snp, $chr, $pos, $beta, $se, $pvals{$snp}), "\n";
    
    # add this SNP's proxies to the seen list
    if ($tags{$snp}) {
        for my $tag (split /\|/, $tags{$snp}) {
            $seen{$tag} = 1;
        }
    }
    else {
        warn "$snp not found in the tags file\n";
    }
    
    $ind++;
}

warn "Found $ind independent SNPs\n";

