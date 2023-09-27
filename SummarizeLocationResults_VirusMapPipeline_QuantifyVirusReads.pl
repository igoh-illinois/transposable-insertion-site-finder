#!/usr/bin/perl

# ************
# Author: J. Cristobal Vera
# email: jcvera@illinois.edu

#description: parse individual sample results from the Run_VirusPipeline_QuantifyReads.pl pipeline and combine virus insertion
#bin and overlap info into a summary table. Intended as an alternative virus insert location summary method for multiple isolate samples, but it
#can be used for deep sequencing samples as well (although less useful)

#use strict;
use Cwd;
use File::Spec;

# defaults, initializations, and constants
my $help = "\n\nSummarizeLocationResults_VirusMapPipeline_QuantifyVirusReads:\nparse individual sample results from the Run_VirusMapPipeline_QuantifyReads.pl pipeline into a summary table.\n".
            "\t-d  Option: Input directory. Required.\n".
            "\t-o  Option: Output file. Default=STDOUT.\n".
            "\t-s  Option: Suggested tie-breaker bin position to report. Default=970501.\n".
            "\t-v  Option: Virus header. Required.\n".
            "\n************\nAuthor: J. Cristobal Vera, email: jcvera\@illinois.edu\n"; 
my $usage = "\nSummarizeLocationResults_VirusMapPipeline_QuantifyVirusReads.pl -d [Input Dir] -o [Output File] -s [Tie-Breaker Bin Position] -v [Virus Header]\n";
my ($indir,$outfile);
my $outfh = 'STDOUT';
my $i = 0;  ##total file lines parsed
my $n = 0;  ##total files parsed
my $tiepos = 970501;  ##suggested default bin start position to use/report in case of tie during frequency sorting
my $vheader = '';   ##virus header; virus contig name
my (%bins,%quals,%overlaps);

#process command line custom script tags
my %cmds = ReturnCmds(@ARGV);
die "\n$help\n$usage\n" if ($cmds{'h'});
if ($cmds{'o'}) {
  $outfh = 'out';
  $outfile = $cmds{'o'};
  open ($outfh, ">$outfile") or die "Cannot create $outfile: $!\n";
}
$indir = $cmds{'d'} if ($cmds{'d'});
$tiepos = $cmds{'s'} if ($cmds{'s'});
$vheader = $cmds{'v'} if ($cmds{'v'});

##make absolute paths
$indir = File::Spec->rel2abs($indir);

#get all sample files
opendir(INDIR,$indir) or die "Can't open directory: $indir: $!\n";
my @inputfiles = grep {m/_Quantification\.tsv$/ && -f "$indir/$_"} readdir(INDIR);
close(INDIR);
opendir(INDIR,$indir) or die "Can't open directory: $indir: $!\n";
my @qualfiles = grep {m/_genome_results\.txt$/ && -f "$indir/$_"} readdir(INDIR);
close(INDIR);
opendir(INDIR,$indir) or die "Can't open directory: $indir: $!\n";
my @overlapfiles = grep {m/_Overlaps\.tsv$/ && -f "$indir/$_"} readdir(INDIR);
close(INDIR);
@inputfiles = sort @inputfiles;
@qualfiles = sort @qualfiles;
@overlapfiles = sort @overlapfiles;

##parse qualimap files
foreach my $qualfile (@qualfiles){
  my $t = 0;
  my $samplename = $qualfile;
  $samplename =~ s/\/([^\/]+)$/$1/;
  $samplename =~ s/_genome_results\.txt$//;
  $samplename =~ s/\.txt$//;
  open (IN, "<$indir/$qualfile") or die "Cannot open $indir/$qualfile: $!\n";
  $n++;
  while (my $line = <IN>) {
    chomp $line;
    next if ($line eq '');
    if (!$t and $line =~ m/number of mapped reads = ([0-9,]+)/){
      $quals{$samplename}{'mapped'} = $1;
    }
    elsif ($line =~ m/mean coverageData = ([0-9.,]+)X/){
      $quals{$samplename}{'totalcov'} = $1;
    }
    elsif (!$t and $line eq '>>>>>>> Coverage per contig'){
      $t = 1;
    }
    elsif ($t){
      my @line = split /\t/,$line;
      if ($line[1]){
        $quals{$samplename}{$line[1]} = $line[4];
      }
    }
  }
  close(IN);
}
print STDERR "\nTotal Qualimap files parsed: $n\n\n";
$n = 0;

##parse virus overlap files
foreach my $overlapfile (@overlapfiles){
  my $samplename = $overlapfile;
  $samplename =~ s/\/([^\/]+)$/$1/;
  $samplename =~ s/_hostSplits_Overlaps\.tsv$//;
  $samplename =~ s/\.tsv$//;
  
  open (IN, "<$indir/$overlapfile") or die "Cannot open $indir/$overlapfile: $!\n";
  $n++;
  my $top = 0;
  while (my $line = <IN>){
    chomp $line;
    next if ($line eq '');
    if ($line =~ m/^Chromosome/){
      next;
    }
    else{
      $i++;
      my ($chr,$range,$treads) = (split /\t/,$line);
      if ($treads > $top){
        my ($start,$end) = (split /-/,$range);
        $overlaps{$samplename}{2} = $overlaps{$samplename}{1} if (exists $overlaps{$samplename}{1});
        $overlaps{$samplename}{1} = "$chr\t$start\t$treads";
        $top = $treads;
      }
    }
  }
  close(IN);
}
print STDERR "\nTotal Overlap files parsed: $n\n\n";
$n = $i = 0;

##parse virus bin files
foreach my $inputfile (@inputfiles){
  my $samplename = $inputfile;
  $samplename =~ s/\/([^\/]+)$/$1/;
  $samplename =~ s/_virusSplitPairs_Quantification\.tsv$//;
  $samplename =~ s/\.tsv$//;
    
  open (IN, "<$indir/$inputfile") or die "Cannot open $indir/$inputfile: $!\n";
  $n++;
  my ($loc);
  my $x = my $y = my $z = my $top1 = my $top2 = my $top3 = my $top4 = 0;
  while (my $line = <IN>) {
    chomp $line;
    next if ($line eq '');
    if ($line =~ m/^For split reads/){
      $loc = 1;
      next;
    }
    elsif ($line =~ m/^For discordant read mates/){
      $loc = 2;
      next;
    }
    elsif ($line =~ m/VirusPosition/){
      next;
    }
    elsif ($loc >= 1){
      $i++;
      my ($readtype,$chr,$bin,$start,$freq) = (split /\t/,$line);
      $x++ if ($loc == 1 and $readtype ne 'center');
      $y++ if ($loc == 2 and $readtype ne 'center');
      $z++ if ($readtype eq 'center');
      if ($loc == 1 and $readtype eq 'onBorder' and $freq > $top1){
        $bins{$samplename}{'splits2'}{$readtype} = $bins{$samplename}{'splits'}{$readtype} if (exists $bins{$samplename}{'splits'}{$readtype});
        $bins{$samplename}{'splits'}{$readtype} = "$chr\t$start\t$freq";
        $top1 = $freq;
      }
      elsif ($loc == 1 and $readtype eq 'onBorder' and $freq == $top1 and $start == $tiepos){
        $bins{$samplename}{'splits2'}{$readtype} = $bins{$samplename}{'splits'}{$readtype} if (exists $bins{$samplename}{'splits'}{$readtype});
        $bins{$samplename}{'splits'}{$readtype} = "$chr\t$start\t$freq";
      }
      elsif ($loc == 1 and $readtype eq 'nearBorder' and $freq > $top2){
        $bins{$samplename}{'splits2'}{$readtype} = $bins{$samplename}{'splits'}{$readtype} if (exists $bins{$samplename}{'splits'}{$readtype});
        $bins{$samplename}{'splits'}{$readtype} = "$chr\t$start\t$freq";
        $top2 = $freq;
      }
      elsif ($loc == 1 and $readtype eq 'nearBorder' and $freq == $top2 and $start == $tiepos){
        $bins{$samplename}{'splits2'}{$readtype} = $bins{$samplename}{'splits'}{$readtype} if (exists $bins{$samplename}{'splits'}{$readtype});
        $bins{$samplename}{'splits'}{$readtype} = "$chr\t$start\t$freq";
      }
      elsif ($loc == 2 and $readtype eq 'onBorder' and $freq > $top3){
        $bins{$samplename}{'paired2'}{$readtype} = $bins{$samplename}{'paired'}{$readtype} if (exists $bins{$samplename}{'paired'}{$readtype});
        $bins{$samplename}{'paired'}{$readtype} = "$chr\t$start\t$freq";
        $top3 = $freq;
      }
      elsif ($loc == 2 and $readtype eq 'onBorder' and $freq == $top3 and $start == $tiepos){
        $bins{$samplename}{'paired2'}{$readtype} = $bins{$samplename}{'paired'}{$readtype} if (exists $bins{$samplename}{'paired'}{$readtype});
        $bins{$samplename}{'paired'}{$readtype} = "$chr\t$start\t$freq";
      }
      elsif ($loc == 2 and $readtype eq 'nearBorder' and $freq > $top4){
        $bins{$samplename}{'paired2'}{$readtype} = $bins{$samplename}{'paired'}{$readtype} if (exists $bins{$samplename}{'paired'}{$readtype});
        $bins{$samplename}{'paired'}{$readtype} = "$chr\t$start\t$freq";
        $top4 = $freq;
      }
      elsif ($loc == 2 and $readtype eq 'nearBorder' and $freq == $top4 and $start == $tiepos){
        $bins{$samplename}{'paired2'}{$readtype} = $bins{$samplename}{'paired'}{$readtype} if (exists $bins{$samplename}{'paired'}{$readtype});
        $bins{$samplename}{'paired'}{$readtype} = "$chr\t$start\t$freq";
      }
    }
  }
  close(IN);
  $bins{$samplename}{'splitbins'} = $x;
  $bins{$samplename}{'pairedbins'} = $y;
  print STDERR "\nTotal sample $samplename:\n\tsplit border positions/bins parsed: $x\n\tpaired border positions/bins parsed: $y\n\tcenter positions/bins: $z\n";
}
print STDERR "\nTotal Quantification files found: $n\nTotal border positions/bins parsed: $i\n\n";

##print out summary
print $outfh "SampleName\tMappedReads\tTotalCoverage\tVirusCoverage\tContig_Overlaps\tStart_Overlaps\tCoverage_Overlaps\tSplitBins\tContig_SplitsOnBorder\tBinStart_SplitsOnBorder\tFrequency_SplitsOnBorder\tContig_SplitsNearBorder\tBinStart_SplitsNearBorder\tFrequency_SplitsNearBorder\t";
print $outfh "PairedBins\tContig_ReadsOnBorder\tBinStart_PairedOnBorder\tFrequency_PairedOnBorder\tContig_PairedNearBorder\tBinStart_PairedNearBorder\tFrequency_PairedNearBorder\tContig_SplitsOnBorder_2nd\tBinStart_SplitsOnBorder_2nd\tFrequency_SplitsOnBorder_2nd\t";
print $outfh "Contig_SplitsNearBorder_2nd\tBinStart_SplitsNearBorder_2nd\tFrequency_SplitsNearBorder_2nd\tContig_ReadsOnBorder_2nd\tBinStart_PairedOnBorder_2nd\tFrequency_PairedOnBorder_2nd\tContig_PairedNearBorder_2nd\tBinStart_PairedNearBorder_2nd\tFrequency_PairedNearBorder_2nd\n";
foreach my $sample (sort keys %bins){
  if ($bins{$sample}{'splitbins'}){
    print $outfh "$sample";
    print $outfh "\t$quals{$sample}{'mapped'}\t$quals{$sample}{'totalcov'}\t$quals{$sample}{$vheader}";
    print $outfh "\t$overlaps{$sample}{'1'}" if (exists $overlaps{$sample}{'1'});
    print $outfh "\t\t\t" if (!exists $overlaps{$sample}{'1'});
    print $outfh "\t$bins{$sample}{'splitbins'}";
    print $outfh "\t$bins{$sample}{'splits'}{'onBorder'}" if (exists $bins{$sample}{'splits'}{'onBorder'});
    print $outfh "\t"x4 if (!exists $bins{$sample}{'splits'}{'onBorder'});
    print $outfh "\t$bins{$sample}{'splits'}{'nearBorder'}" if (exists $bins{$sample}{'splits'}{'nearBorder'});
    print $outfh "\t"x3 if (!exists $bins{$sample}{'splits'}{'nearBorder'});
  }
  else{
    print $outfh "$sample\t$quals{$sample}{'mapped'}\t$quals{$sample}{'totalcov'}\t$quals{$sample}{$vheader}\t$bins{$sample}{'splitbins'}\t".("\t"x5);
  }
  if ($bins{$sample}{'pairedbins'}){
    print $outfh "\t$bins{$sample}{'pairedbins'}";
    print $outfh "\t$bins{$sample}{'paired'}{'onBorder'}" if (exists $bins{$sample}{'paired'}{'onBorder'});
    print $outfh "\t"x3 if (!exists $bins{$sample}{'paired'}{'onBorder'});
    print $outfh "\t$bins{$sample}{'paired'}{'nearBorder'}" if (exists $bins{$sample}{'paired'}{'nearBorder'});
    print $outfh "\t"x3 if (!exists $bins{$sample}{'paired'}{'nearBorder'});
  }
  else{
    print $outfh "\t$bins{$sample}{'pairedbins'}\t".("\t"x5);
  }
  if (exists $bins{$sample}{'splits2'}{'onBorder'}){
    print $outfh "\t$bins{$sample}{'splits2'}{'onBorder'}";
  }
  else{
    print $outfh "\t"x3;
  }
  if (exists $bins{$sample}{'splits2'}{'nearBorder'}){
    print $outfh "\t$bins{$sample}{'splits2'}{'nearBorder'}";
  }
  else{
    print $outfh "\t"x3;
  }
  if (exists $bins{$sample}{'paired2'}{'onBorder'}){
    print $outfh "\t$bins{$sample}{'paired2'}{'onBorder'}";
  }
  else{
    print $outfh "\t"x3;
  }
  if (exists $bins{$sample}{'paired2'}{'nearBorder'}){
    print $outfh "\t$bins{$sample}{'paired2'}{'nearBorder'}\n";
  }
  else{
    print $outfh ("\t"x3)."\n";
  }
}
print STDERR "\nFinished!\n\n";


sub ReturnCmds{
  my (@cmds) = @_;
  my ($opt);
  my %cmds;
  foreach my $cmd (@cmds){
    if (!$opt and $cmd =~ m/^-([a-zA-Z])/) {
      $opt = $1;
    }
    elsif ($opt and $cmd =~ m/^-([a-zA-Z])/){
      $cmds{$opt} = 1;
      $opt = $1;
    }
    elsif ($opt){
      $cmds{$opt} = $cmd;
      $opt = '';
    }
  }
  $cmds{$opt} = 1 if ($opt);
  return %cmds;
}
