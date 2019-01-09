#!/usr/bin/perl -w
use strict;
use warnings;
use FindBin qw($Bin);
#use Config::IniFiles;
use File::Basename;
use Getopt::Long;
use Cwd qw(abs_path getcwd);

my $trinity="/data02/software/trinity/trinityrnaseq-Trinity-v2.8.2/";
my $R="/data02/software/R/R-v3.5.0/bin/";

my ($sample,$contrasts,$outdir,$dispersion,$matrix,$foldchange,$FDR);
$outdir||=getcwd;
$outdir=abs_path($outdir);
$FDR||=0.001;
$foldchange||=2;
$dispersion||=0.1;
GetOptions(
    "g:s"=>\$sample,
    "v:s"=>\$contrasts,
    "o:s"=>\$outdir,
    "m:s"=>\$matrix,
    "d:s"=>\$dispersion,
    "C:s"=>\$foldchange,
    "P:s"=>\$FDR,
           );
if(!$sample||!$contrasts||!$matrix)
{
    &usage();
}

sub usage{
    print qq{
This script will run DGE.
usage:
perl $0 -g group.txt -v vs.txt -m matrix.txt -d 0.00001 -c 2 -f 0.05

options:
-g          tab-delimited text file indicating biological replicate relationships.
                cond_A  sample1
                cond_A  sample2
                cond_B  sample3
                cond_B  sample4
-v          tab-delimited text file containing the pairs of sample comparisons to perform.
                cond_A  cond_B
                cond_Y  cond_Z
-m          matrix of raw read counts (not normalized!)
-P          p-value cutoff for FDR  (default: 0.001)
-C          min abs(log2(a/b)) fold change (default: 2(meaning 2^(2) or 4-fold).
-d          edgeR dispersion value :0.4 for human data,0.1 for data on genetically identical model organisms or 0.01 for technical replicates.
-o          output directory

Email:fanyucai1\@126.com
2018.8
    };
    exit;
}
system "mkdir -p $outdir";
$matrix=abs_path($matrix);
$sample=abs_path($sample);
$contrasts=abs_path($contrasts);
$outdir=abs_path($outdir);
`export PATH=$R:\$PATH && perl $trinity/Analysis/DifferentialExpression/run_DE_analysis.pl  --matrix $matrix --method edgeR --dispersion $dispersion --samples_file $sample --contrasts $contrasts --output $outdir`;
`export PATH=$R:\$PATH && cd $outdir && perl $trinity/Analysis/DifferentialExpression/run_TMM_normalization_write_FPKM_matrix.pl --matrix $matrix --just_do_TMM_scaling`;
`export PATH=$R:\$PATH && cd $outdir && perl $trinity/Analysis/DifferentialExpression/analyze_diff_expr.pl --samples $sample --matrix $matrix.TMM_rescaled -P $FDR -C $foldchange --output $outdir`;


open(VS1,"$contrasts");
while(<VS1>)
{
	chomp;
	my @array=split;
	my $string;
	$string.=basename ($matrix);
	$string.=".$array[0]\_vs_$array[1].";
	$string.="edgeR.DE_results";
	`cd $outdir && perl $Bin/plot_Volcano_MA.pl -i $string -o $outdir -vs $array[0].$array[1]`;
}
