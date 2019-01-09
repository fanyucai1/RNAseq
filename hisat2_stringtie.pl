#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin);
use Getopt::Long;
use Cwd;
use Cwd 'abs_path';

my $hisat2="/local_data1/software/hisat2/hisat2-2.1.0/";
my $stringtie="/local_data1/software/stringtie/stringtie-1.3.4d.Linux_x86_64";
my $samtools="/local_data1/software/samtools/samtools-1.9/";
my $R="/local_data1/software/R/R-v3.4.0/bin/";
my $gffcompare="/local_data1/software/GffCompare/gffcompare/";
my $qsub="/home/fanyucai/software/qsub/qsub-sge.pl";
my $gffread="/local_data1/software/gffread/gffread-0.9.12.Linux_x86_64/";
my $preDE="/local_data1/work/fanyucai/script/prepDE.py";

my (@pe1,@pe2,@prefix,$outdir,$ref,$gtf,$lib,$readlength,$k);
$outdir||=getcwd;
$readlength||=75;
$k||="all"；
GetOptions(
    "pe1:s{1,}"=>\@pe1,       
    "pe2:s{1,}"=>\@pe2,
    "p:s{1,}"=>\@prefix,
    "o:s"=>\$outdir,
    "r:s"=>\$ref,
    "l:s"=>\$lib,
    "g:s"=>\$gtf,
    "len:s"=>\$readlength,
    "k:s"=>\$k,
           );

sub usage{
    print qq{
This script will run the hisat2 and stringtie.
usage:
perl $0 -pe1 s1.1.fq(.gz) s2.1.fq(.gz) -pe2 s1.2.fq(.gz) s2.2.fq(.gz) -o $outdir -p s1 s2 -r reference.fna -g $gtf -l FR

options:
-pe1        5 reads several split by space:force
-pe2        3 reads several split by space:force
-p          prefix output several split by space:forcre
-r          hisat2 index reference fasta sequence:force
-l          strand-specific,choices=["RF", "FR", "unstranded"]:force
-g          gtf file
-o          output directory
-len        the average read length[default: 75]
-k          only quality known transcript(known,k) or all denovo trancript(all,a),defualt:a

Email:fanyucai1\@126.com
2018.8
    };
    exit;
}

if(!@pe1|| !@pe2 ||!@prefix ||!$readlength)
{
    &usage();
}

system "mkdir -p $outdir";
$gtf=abs_path($gtf);
$outdir=abs_path($outdir);
#1: mapping using hisat2,convert sam to bam and sort the bam file
open(SH1,">$outdir/hisat2.sh");
for(my $i=0;$i<=$#pe1;$i++)
{
    system "mkdir -p $outdir/$prefix[$i]";
    $pe1[$i]=abs_path($pe1[$i]);
    $pe2[$i]=abs_path($pe2[$i]);
    print SH1 "$hisat2/hisat2 -x $ref -p 20 -1 $pe1[$i] -2 $pe2[$i] --rna-strandness $lib -S $outdir/$prefix[$i]/$prefix[$i].sam && ";
    print SH1 "$samtools/samtools view -@ 20 -uS $outdir/$prefix[$i]/$prefix[$i].sam|$samtools/samtools sort -@ 20 -o $outdir/$prefix[$i]/$prefix[$i].bam\n";
}
if(! -e "$outdir/$prefix[$#prefix]/$prefix[$#prefix].bam")
{
    `perl $qsub $outdir/hisat2.sh`;
}

#2: using stringtie to assemble RNA-Seq alignments into potential transcripts
open(SH2,">$outdir/stringtie.sh");
open(LIST,">$outdir/merge.list");
for(my $i=0;$i<=$#pe2;$i++)
{
    if($lib=~/RF/)
    {
        print SH2 "$stringtie/stringtie $outdir/$prefix[$i]/$prefix[$i].bam --rf -p 20 -G $gtf -l $prefix[$i] -m 100 -o $outdir/$prefix[$i]/$prefix[$i].out\n";
    }
    elsif($lib=~/FR/)
    {
        print SH2 "$stringtie/stringtie $outdir/$prefix[$i]/$prefix[$i].bam --fr -p 20 -G $gtf -l $prefix[$i] -m 100 -o $outdir/$prefix[$i]/$prefix[$i].out\n";
    }
    else
    {
      print SH2 "$stringtie/stringtie $outdir/$prefix[$i]/$prefix[$i].bam -p 20 -G $gtf -l $prefix[$i] -m 100 -o $outdir/$prefix[$i]/$prefix[$i].out\n";
    }
    print LIST "$outdir/$prefix[$i]/$prefix[$i].out\n";
}
if(! -e "$outdir/$prefix[$#prefix]/$prefix[$#prefix].out")
{
    `perl $qsub $outdir/stringtie.sh`;
}

#3: Merged GTF and compare use cuffcompare to compare gtf file
`$stringtie/stringtie --merge -G $gtf -p 20 -o $outdir/stringtie_merge.gtf $outdir/merge.list`;
`cd $outdir && $gffcompare/gffcompare -R -r $gtf -o strtcmp stringtie_merge.gtf`;

#4:Extracting transcript sequences
if(! -e "$ref.fai")
{
    `$samtools faidx $ref`;
}
if(! -e "$outdir/transcripts.fa")
{
    `$gffread/gffread -w $outdir/transcripts.fa -g $ref stringtie_merge.gtf`;
}

#5:Estimate transcript abundances and create table counts for Ballgown
open(SH3,">$outdir/ballgown.sh");
open(GTF,">$outdir/gtf.list");
my $gtftmp;
if($k=~"a")
{
    $gtftmp=$outdir/stringtie_merge.gtf;
}
else
{
    $gtftmp=$gtf;
}
for(my $i=0;$i<=$#pe2;$i++)
{
    if($lib=~/RF/)
    {
        print SH3 "$stringtie/stringtie $outdir/$prefix[$i]/$prefix[$i].bam --rf -e -B -p 20 -F 0.01 -m 100 -A $outdir/$prefix[$i]/$prefix[$i].tab -G $gtftmp -o $outdir/$prefix[$i]/$prefix[$i].gtf\n";
    }
    elsif($lib=~/FR/)
    {
        print SH3 "$stringtie/stringtie $outdir/$prefix[$i]/$prefix[$i].bam --fr -e -B -p 20 -F 0.01 -m 100 -A $outdir/$prefix[$i]/$prefix[$i].tab -G $gtftmp -o $outdir/$prefix[$i]/$prefix[$i].gtf\n";
    }
    else
    {
        print SH3 "$stringtie/stringtie $outdir/$prefix[$i]/$prefix[$i].bam -e -B -p 20 -F 0.01 -m 100 -A $outdir/$prefix[$i]/$prefix[$i].tab -G $gtftmp -o $outdir/$prefix[$i]/$prefix[$i].gtf\n";
    }
    print GTF "$prefix[$i]\t$outdir/$prefix[$i]/$prefix[$i].gtf\n";
}
if(! -e "$outdir/$prefix[$#prefix]/$prefix[$#prefix].gtf")
{
    `perl $qsub $outdir/ballgown.sh`;
}
`python /local_data1/work/fanyucai/script/prepDE.py -l $readlength -i $outdir/gtf.list -g $outdir/gene_count_matrix.csv -t $outdir/transcript_count_matrix.csv`;
`sed -i "s:,:\t:g" $outdir/gene_count_matrix.csv`;
`sed -i "s:,:\t:g" $outdir/transcript_count_matrix.csv`;







