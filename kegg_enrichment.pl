#!/usr/bin/perl -w
use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Long;
use Cwd;
use Cwd 'abs_path';

my $kobas="/data02/software/kobas/kobas-3.0/";
my $blastplus="/data02/software/blast+/ncbi-blast-2.7.1+/bin";
my $R="/data02/software/R/R-v3.5.0/";
my $python="/data02/software/python/Python-v2.7.15/bin/";
my $pcre="/data02/software/pcre/pcre-v8.42/lib/";
my $env="export PATH=$R/bin:$python:\$PATH";
my $lib="export LD_LIBRARY_PATH=$R/lib/R/lib:$pcre:\$LD_LIBRARY_PATH";

my($input,$type,$outdir,$species,$prefix,$FDR);
$species||="hsa";
$outdir||=getcwd;
$FDR||=0.1;
my $usr=`whoami`;
chomp $usr;
my $dir="/home/$usr/.kobasrc";
system "echo '
[DEFAULT]
kobas_home =$kobas
blast_home =$blastplus
[KOBAS]
kobasdb =$kobas/sqlite3/
gmt = $kobas/gmt/
grn = $kobas/grn/
model =$kobas/model/

[BLAST]
blastp = $blastplus/blastp
blastx = $blastplus//blastx
blastdb = $kobas/seq_pep/
' >$dir ";

GetOptions(
    "i:s"=>\$input,
    "t:s"=>\$type,
    "o:s"=>\$outdir,
    "p:s"=>\$prefix,
    "sp:s"=>\$species,
           );
sub usage{
    print qq{
This script will run enrichment using kobas.
options:
-i          input data file
-t          input type
            (fasta:pro, fasta:nuc, blastout:xml,blastout:tab,
             id:ncbigi, id:uniprot, id:ensembl,id:ncbigene), default fasta:pro
-o          output directory(defualt:$outdir)
-sp         species abbreviation(defualt:hsa)
-p          prefix of output

Email:fanyucai1\@126.com
2018.9.4
    };
    exit;
}
if(!$input||!$type||!$prefix||!$species)
{
    &usage();
}
$input=abs_path($input);
$outdir=abs_path($outdir);

my $par=" -k $kobas -v $blastplus -y $kobas/sep_pep/ -q $kobas/sqlite3/ -p $blastplus/blastp -x $blastplus/blastx ";
##############1:annotate#############################
system "echo '$env && $lib && python $kobas/scripts/annotate.py -i $input -t $type -n 30 -s $species -o $outdir/$prefix.annotate $par'>$outdir/anno.sh";
`sh $outdir/anno.sh`;

##############2:identify############################
system "echo '$env && $lib &&  python $kobas/scripts/identify.py -f $outdir/$prefix.annotate -b $species -d K -o $outdir/$prefix.kegg.identify $par'>$outdir/identify.sh";
system "echo '$env && $lib &&  python $kobas/scripts/identify.py -f $outdir/$prefix.annotate -b $species -d G -o $outdir/$prefix.go.identify $par'>>$outdir/identify.sh";
`sh $outdir/identify.sh`;

#########extract significant pathways#########
open(KEGG,"$outdir/$prefix.kegg.identify");
open(OUT1,">$outdir/$prefix.kegg.sig.tsv");
print OUT1 "#Term\tInput_number\tEnrichment Factor\tFDR\tGeneID\n";
while(<KEGG>)
{
    chomp;
    if($_!~/^#/)
    {
        if($_=~/[A-Z]/)
        {
            my @array=split("\t",$_);
            if($array[6]<=$FDR)
            {
                my $en=$array[3]/$array[4];
                print OUT1 "$array[0]\t$array[3]\t$en\t$array[6]\t$array[7]\n";
            }
        }
    }
}

open(GO,"$outdir/$prefix.go.identify");
open(OUT2,">$outdir/$prefix.go.sig.tsv");
print OUT2 "#Term\tInput_number\tEnrichment Factor\tFDR\tGeneID\n";
while(<GO>)
{
    chomp;
    if($_!~/^#/)
    {
        if($_=~/[A-Z]/)
        {
            my @array=split("\t",$_);
            if($array[6]<=$FDR)
            {
                my $en=$array[3]/$array[4];
                print OUT2 "$array[0]\t$array[3]\t$en\t$array[6]\t$array[7]\n";
            }
        }
    }
}