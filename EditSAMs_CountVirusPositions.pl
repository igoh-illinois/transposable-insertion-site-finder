#!/usr/bin/perl

#use strict;
# ************
# Author: J. Cristobal Vera
# email: jvera8888@gmail.com

##Description: parse virus/host alignments (discordant mates and split reads spanning the virus/host insertion site) in SAM file
##and return quantification in host by bin. Basically asks: where do discordant mates and split reads from different parts of the
##inserted phage go to the host?  SAM file should only contain virus alignments (i.e. filter the SAM, see below).

##NOTES:
##-multiple splits to a read are currently ignored (except the first in list, which is usually the longest)
##-filter SAM so that it contains only virus genome alignments with discordant mates and/or split reads spanning the virus/host insertion site
##--e.g.
##   **first concatenate the host and virus genomes into a single multi-fasta file and generate the [host/virus ref] using the bwa index command, then...
##   bwa mem -t 8 [host/virus ref] [trimfile1].fastq [trimfile2].fastq | sambamba view -F "not unmapped" -S -f bam /dev/stdin >[bamout].bam
##   sambamba sort -t 8 -m 64G -F "ref_name =~ /[Virus Header]/ and not [SA] == null and not [SA] =~ /[Virus Header]/ and [XA] == null" -o [bamout].sortfilter_splits.bam [bamout.bam] && sambamba view -h [bamout].sortfilter_splits.bam >[bamout].sortfilter_splits.sam
##   sambamba sort -t 8 -m 64G -F "ref_name =~ /[Virus Header]/ and mate_ref_name =~ /^[Host Header Base]/ and [XA] == null and [SA] == null" -o [bamout].sortfilter_mates.bam [bamout].bam | sambamba view -h [bamout].sortfilter_mates.bam >[bamout].sortfilter_mates.sam
##   perl CombineSAMs.pl -i [bamout].sortfilter_splits.sam -j [bamout].sortfilter_mates.sam -o [bamout].sortfilter_SplitsAndMates.sam
##   **run this script on [bamout].sortfilter_SplitsAndMates.sam!


# defaults, initializations, and constants
my $help = "\n\nEditSAMs_CountVirusPositions.\nDescription: parse virus/host alignments (discordant mates and split reads spanning the virus/host insertion site) in SAM file and return quantification in host by bin. Basically asks: where do discordant mates and split reads from different parts of the inserted phage go to the host?\n".
            "\t-i  Option: Input file. Required.\n".
            "\t-o  Option: Output file. Default=STDOUT.\n".
            "\t-c  Option: virus contig/chromosome header. Default='NC_008717.1'.\n".
            "\t-b  Option: host chromosome bin size. Default=1000.\n".
            "\n************\nAuthor: J. Cristobal Vera, email: jvera8888\@gmail.com\n"; 
my $usage = "\nEditSAMs_CountVirusPositions.pl -i [Input File] -o [Output File] -c [Virus Header] -b [Host Bin Size]\n";
my $outfh = my $infh = 'STDOUT';
my $i = my $j = my $x = my $y = my $z = 0;
my $virus_header = 'NC_008717.1';
my $virus_end_len = 500;  ##length from end considered cannonical; depends on read lengths; used for distance from virus ends and for binning sizes
my $bin_size = 1000;
my $virus_len;
my %contig_lens;
my %contig_bins;
my %virus_dist_splits;
my %virus_dist_mates;


#process command line custom script tags
my %cmds = ReturnCmds(@ARGV);
die "\n$help\n$usage\n" if ($cmds{'h'});
if ($cmds{'i'}) {
  $infh = 'IN';
  open ($infh, "<$cmds{'i'}") or die "Cannot open $cmds{'i'}: $!\n";
}
if ($cmds{'o'}) {
  $outfh = 'OUT';
  open ($outfh, ">$cmds{'o'}") or die "Cannot open $cmds{'o'}: $!\n";
}
$virus_header = $cmds{'c'} if ($cmds{'c'});
$bin_size = $cmds{'b'} if ($cmds{'b'});

while (my $line = <$infh>){
  chomp $line;
  next if ($line eq '');
  $i += 1;
  if ($line =~ m/^\@SQ/){
    my ($tag,$chr,$len) = (split /\t/,$line);
    $chr =~ s/^SN://;
    $len =~ s/^LN://;
    if ($chr eq $virus_header){
      $virus_len = $len;
    }
    else{
      $contig_lens{$chr} = $len;
      ##create bin categories for each chromosome/contig
      my $contig_bins = int($len/$bin_size)+1; ##get number of bins on chromosome plus 1 for remainder
      my $start = 1;
      #my $stop = $bin_size;
      for (my $c = 1;$c <= $contig_bins;$c++){
        $contig_bins{$chr}{$c} = $start; ##add start positions of bins to bin hash
        $start += $bin_size;
      }
    }
  }
  elsif ($line =~ m/^\@/){
    ##do nothing
  }
  else{
    $j += 1;
    my @tags = split /\t/,$line;
    my ($read,$flag,$ref1,$pos1,$mapq,$cigar,$ref2,$pos2,$tlen,$seq,$qual) = (splice @tags,0,11);  ##leave tags in @tags
    die "\nError: read is not within virus contig:\n\tvirus: $virus_header\n\tref: $ref1\n\n" if ($ref1 ne $virus_header);  ##check that all reads in SAM are mapped to the virus
    my $len1 = ParseCigar('M',$cigar);
    my $read_virus_pos = CheckVirusPos($pos1,$len1,$virus_len,$virus_end_len);
    my $split = CheckSplit(@tags);
    my $mate = 0;
    $mate = 1 if ($ref2 ne '=');
    die "\nError: read alignment is neither discordant nor split across contigs at $j:\nline:$line\n" if (!$mate and !$split); ##check that all reads are split or paired across virus/host boundry
    if ($split){
      my ($ref3,$pos3,$len3) = (split /\t/,$split);
      my $chr_pos_bin = GetBin($ref3,$pos3,$contig_lens{$ref3},\%contig_bins);
      my $keys = keys %{$contig_bins{$ref3}};
      die "\nError: binning error1: $line\n\n" if (!$chr_pos_bin);
      $virus_dist_splits{$read_virus_pos}{$ref3}{$chr_pos_bin} += 1;
      $x += 1;
    }
    if ($mate){
      my $chr_pos_bin = GetBin($ref2,$pos2,$contig_lens{$ref2},\%contig_bins);
      die "\nError: binning error2: $line\n\n" if (!$chr_pos_bin);
      $virus_dist_mates{$read_virus_pos}{$ref2}{$chr_pos_bin} += 1;
      $y += 1;
    }
    if ($mate and $split){
      $z += 1;
    }
  }
}
##print out
print $outfh "For split reads (host/virus):\nVirusPositionType\tChromosome\tBin\tStartPosition\tFrequency\n";
foreach my $virus_pos_type (keys %virus_dist_splits){
  foreach my $chr (sort keys %{$virus_dist_splits{$virus_pos_type}}){
    foreach my $bin (sort keys %{$virus_dist_splits{$virus_pos_type}{$chr}}){
      print $outfh "$virus_pos_type\t$chr\t$bin\t$contig_bins{$chr}{$bin}\t$virus_dist_splits{$virus_pos_type}{$chr}{$bin}\n";
    }
  }
}
print $outfh "\nFor discordant read mates (host/virus):\nVirusPositionType\tChromosome\tBin\tStartPosition\tFrequency\n";
foreach my $virus_pos_type (keys %virus_dist_mates){
  foreach my $chr (sort keys %{$virus_dist_mates{$virus_pos_type}}){
    foreach my $bin (sort keys %{$virus_dist_mates{$virus_pos_type}{$chr}}){
      print $outfh "$virus_pos_type\t$chr\t$bin\t$contig_bins{$chr}{$bin}\t$virus_dist_mates{$virus_pos_type}{$chr}{$bin}\n";
    }
  }
}

print STDERR "\nTotal lines read: $i\n";
print STDERR "Total alignments: $j\n";
print STDERR "Total discordant read mates (virus/host): $y\n";
print STDERR "Total split reads (virus/host): $x\n";
print STDERR "Total reads that are both discordant and split: $z\n";

##assign host chromosome bin/position
sub GetBin{
  my ($chr,$chr_pos,$contig_len,$contig_bins) = @_;
    my $n = 0;
    foreach my $bin (keys %{$contig_bins{$chr}}){
      my $start = $contig_bins->{$chr}{$bin};
      my $stop = $contig_len;
      $stop = $contig_bins->{$chr}{($bin+1)} if (exists $contig_bins->{$chr}{($bin+1)});
      return $bin if ($chr_pos >= $start and $chr_pos <= ($stop-1));
      $n++;
    }
    #for (my $c = 1;$c <= $contig_bins;$c++){
    #  for (my $q = 1;$q <= $bin_size;$q++){
    #    $n++;
    #    return $p if ($chr_pos == $n);
    #  }
    #}
    #print "c:$chr\tcl:$contig_len\tcp:$chr_pos\tsa:$start\tst:$stop\n";
    return 0;
}

##check for split read tag and return parsed alignment info
sub CheckSplit{
  my (@tags) = @_;
  foreach my $tag (@tags){
    if ($tag =~ m/^SA/){
      my @SA = split /;/,$tag;
      if ($SA[0] =~ m/SA:Z:([^,]+)\,([0-9]+)\,[+-]\,([0-9SHM]+)/){
        my ($chr,$start,$cigar) = ($1,$2,$3);
        my $len = ParseCigar('M',$cigar);
        return "$chr\t$start\t$len";
      }
    }
  }
  return 0;
}

##parse SAM file CIGAR string and return specified part (either M, S, or H)
sub ParseCigar{
  my ($part,$cigar) = @_;
  if ($cigar =~ m/[SHM]{0,1}([0-9]+)$part/){
    return $1;
  }
  return 0;
}

sub CheckVirusPos{
  my ($pos,$len,$virus_len,$virus_end_len) = @_;
  ##check if read start is on or very close to virus/host border
  if ($pos <= 5 or (($pos+$len) >= ($virus_len-5))){
    return "onBorder";
  }
  ##check if read start is close to virus/host border (but not on it)
  elsif ($pos <= $virus_end_len or (($pos+$len) >= ($virus_len-$virus_end_len))){
    return "nearBorder";
  }
  else{
    return "center";
  }
}

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
