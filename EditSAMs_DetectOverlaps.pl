#!/usr/bin/perl

#use strict;
# ************
# Author: J. Cristobal Vera
# email: jcvera@illinois.edu

##Description: parse virus/host alignments (split reads spanning the virus/host insertion site) in SAM file
##and returns position and depth of overlap regions. Overlap regions represent the duplicated host sequences at the exact Mu phage insertion site.
##SAM file should only contain host split read alignments that span the virus insertion junctions/borders (i.e. filter the SAM, see below).

##NOTES:
##-multiple splits to a read are currently ignored (except the first in list, which is usually the longest)
##-this script currently filters overlaps only for read support, not directionality of read splits. it is recommended that overlaps be post-filtered for at least one left and one right split occurence
##-filter SAM so that it contains only unique HOST genome alignments with split reads spanning the virus/host
## insertion site
##--e.g.
##   **first concatenate the host & virus genomes into a single multi-fasta file and generate the [host/virus ref] using the bwa index command, then...
##   bwa mem -t 8 [host/virus ref] [trimfile1].fastq [trimfile2].fastq | sambamba view -S -f bam /dev/stdin >[bamout].bam
##   sambamba sort -t 8 -m 64G -F "not unmapped and not ref_name =~ /[Virus Header]/ and [SA] =~ /[Virus Header]/ and [XA] == null" -o [bamout].sortfilter_HostSplits.bam [bamout.bam] && sambamba view -h [bamout].sortfilter_HostSplits.bam >[bamout].sortfilter_Hostplits.sam
##   **run this script on [bamout].sortfilter_HostSplits.sam


# defaults, initializations, and constants
my $help = "\n\nEditSAMs_DetectOverlaps.\nDescription: parse virus/host alignments (split reads spanning the virus/host insertion site) in SAM file and returns location and depth of overlap regions.\n".
            "\t-i  Option: Input file. Default=STDIN.\n".
            "\t-o  Option: Output file. Default=STDOUT.\n".
            "\t-c  Option: virus contig/chromosome header. Default='NC_008717.1'.\n".
            "\t-g  Option: genome/reference tab-delimited sequence file. Required.\n".
            "\t-f  Option: minimum read support filter value. Default=2.\n".
            "\n************\nAuthor: J. Cristobal Vera, email: jcvera\@illinois.edu\n";
my $usage = "\nEditSAMs_DetectOverlaps.pl -i [Input File] -o [Output File] -c [Virus Header] -g [Genome File] -f [Min Reads]\n";
my $outfh = my $infh = 'STDOUT';
my $i = my $j = my $x = my $y = my $z = my $m = my $n = my $c = my $d = my $e = my $f = 0;
my $virus_header = 'NC_008717.1';
my $virus_end_len = 500;  ##used for distance from virus ends; depends on read lengths;
my $overlap_len = 5;      ##overlap region length to check; MU phage = 5bp
my $virus_len;            ##virus genome length
my $reffile;
my $min_reads = 2;        ##minimum reads filter value
my %seqs;                 ##5bp overlap seqs
my %contig_lens;          ##list of chromosome lengths
my %contig_seqs;          ##contig sequences
my %splitreads;           ##list of splitread names with orientation values
my %total_splitends;      ##split read IDs: IDs are location (chromosome & overlap position) of reads that are split across virus/host junction
my %rightend_splitends;   ##split reads with split occuring on right end of read
my %leftend_splitends;    ##split reads with split occuring on left end of read
my %total_con_splitends;  ##split reads from concordant alignments (i.e. host/virus split occurs where expected on read)
my %rightend_con_splitends;##concordant split reads with split occuring on right end of read
my %leftend_con_splitends;##concordant split reads with split occuring on left end of read
my %total_center_splitends;##split reads originating from center of virus
my %total_con_center_splitends;##concordant split reads originating from center of virus

###sam file bits
#my %bits = (
#    1 => '0x0001',      #read was/is a segment of a pair (before mapping)
#    2 => '0x0002',      #both reads are mapped 'properly' (either within reasonable distance or in proper orientation or both, depending on aligner)
#    4 => '0x0004',      #read is unmapped
#    8 => '0x0008',      #read mate is unmapped
#    16 => '0x0010',     #read is reverse complimented (it's minus strand)
#    32 => '0x0020',     #read mate is reverse complimented (it's minus strand)
#    64 => '0x0040',     #first read in fragment
#    128 => '0x0080',    #second read in fragment
#    256 => '0x0100',    #secondary alignment
#    512 => '0x0200',    #failed read
#    1024 => '0x0400',   #duplicate fragment
#    2048 => '0x0800',   #supplementary alignment
#  );

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
if ($cmds{'g'}) {
  open (REF, "<$cmds{'g'}") or die "Cannot open $cmds{'g'}: $!\n";
  while (my $line = <REF>){
    chomp $line;
    $i++;
    my @line = split /\t/,$line;
    $line[0] =~ s/^([^ ]+) .+/$1/;
    $contig_seqs{$line[0]} = $line[1];
  }
  close(REF);
  print STDERR "\nContigs found: $i\n\n";
  $i = 0;
}
$virus_header = $cmds{'c'} if ($cmds{'c'});
$min_reads = $cmds{'f'} if ($cmds{'f'});
die "\nPlease select minimum read support threshold of 2 or higher.\n\n" if ($min_reads <= 1);

while (my $line = <$infh>){
  chomp $line;
  next if ($line eq '');
  $i += 1;
  ##parse headers
  if ($line =~ m/^\@SQ/){
    my ($tag,$chr,$len) = (split /\t/,$line);
    $chr =~ s/^SN://;
    $len =~ s/^LN://;
    if ($chr eq $virus_header){
      $virus_len = $len;
    }
    else{
      ##get contig names and lengths
      $contig_lens{$chr} = $len;
      die "\nError: cannot find contig $chr in genome file.\n\n" if (!exists $contig_seqs{$chr});
    }
  }
  elsif ($line =~ m/^\@/){
    ##do nothing
  }
  ##parse alignments
  else{
    $j += 1;
    my @tags = split /\t/,$line;
    my ($read,$flag,$ref1,$pos1,$mapq,$cigar,$ref2,$pos2,$tlen,$seq,$qual) = (splice @tags,0,11);  ##leave tags in @tags
    ##critical error check, all reads must be on the host, not the virus
    die "\nError: read is not on a host contig. please filter alignments correctly:\n\tvirus: $virus_header\n\tref: $ref1\n\n" if ($ref1 eq $virus_header);
    
    ##get read length and start/stop positions on host chunk
    my ($len1,$cutend1);
    ($pos1,$len1,$cutend1) = GetMatchPart($pos1,$cigar);
    
    ##check to make sure read on host is clipped (it has to be a split read)
    die "\nError: host read has no clipping. please filter alignments correctly:\n\talignment row: $j\n\tRead: $read\n\tLine: $line\n\n" if (!$pos1);
    ##read chunk should not contain deletions
    print STDERR "\nWarning: host read $read chunk at alignment $j contains deletion, skipping read:\n\tcigar: $cigar\n\tline: $line\n\n" if ($pos1 eq 'd');
    $d++ if ($pos1 eq 'd');
    next if ($pos1 eq 'd');
    ##read chunk should not contain insertion in cigar
    print STDERR "\nWarning: host read $read chunk at alignment $j contains insertion, skipping read:\n\tline: $line\n\n" if ($pos1 eq 'i');
    $f++ if ($pos1 eq 'i');
    next if ($pos1 eq 'i');
    ##read should only be split on one end or we can't determine where the virus starts
    print STDERR "\nWarning: skipping read $read at alignment $j: HOST read chunk is clipped on both ends:\nLine: $line\n\n" if (!$cutend1);
    $n++ if (!$cutend1);
    next if (!$cutend1);
    
    ##get split tag and match chunk info for virus chunk
    my ($ref3,$pos3,$len3,$cutend3) = CheckSplit(@tags);
    
    ##critical error: check for no splitting
    die "\nError: read $read at alignment $j is not split: $line\n\n" if (!$ref3);
    ##check for complex splitting
    print STDERR "\nWarning: skipping read $read at alignment $j: not split across virus/host boundry or split is too complex:\nline:$line\n" if ($ref3 ne $virus_header);
    $z++ if ($ref3 ne $virus_header);
    next if ($ref3 ne $virus_header);
    ##read should only be split on one end or we can't determine where the host starts
    print STDERR "\nWarning: skipping read $read at alignment $j: VIRUS read chunk is clipped on both ends:\nLine: $line\n\n" if (!$cutend3);
    $n++ if (!$cutend3);
    next if (!$cutend3);
    ##virus read chunk should not contain deletion in cigar
    print STDERR "\nWarning: virus read $read chunk at alignment $j contains deletion, skipping read:\n\tline: $line\n\n" if ($pos3 eq 'd');
    $d++ if ($pos1 eq 'd');
    next if ($pos1 eq 'd');
    ##virus read chunk should not contain insertion in cigar
    print STDERR "\nWarning: virus read $read chunk at alignment $j contains insertion, skipping read:\n\tline: $line\n\n" if ($pos3 eq 'i');
    $f++ if ($pos1 eq 'i');
    next if ($pos1 eq 'i');
    
    ##check whether or not we've already parsed this read's mate.  if so, we skip it (they are probably overlapped anyways, so this doesn't represent new information)
    print STDERR "\nWarning: skipping read mate for read $read at alignment $j: read mate is also split and probably overlapped:\nLine: $line\n\n" if (exists $splitreads{$read});
    $m++ if (exists $splitreads{$read});
    next if (exists $splitreads{$read});
    ##check for bit parsing errors: read has both 1st and 2nd read bit set, meaning it's pretty messed up or perhaps it's a singleton?  not sure... but either way get rid of it
    print STDERR "\nWarning: skipping read $read at alignment $j due to bit error, possibly not Illumina read: $line\n\n" if ($flag & 64 and $flag & 128);
    $c++ if ($flag & 64 and $flag & 128);
    next if ($flag & 64 and $flag & 128);
    
    ##check for read orientation info. we want unique reads with orientations relative to the virus that make sense
    $splitreads{$read} = 0;
    $splitreads{$read} = '1to2v' if ($ref2 eq $virus_header and $flag & 64); ##this read is first in pair and its mate maps to virus
    $splitreads{$read} = '2to1v' if ($ref2 eq $virus_header and $flag & 128); ##this read is second in pair and its mate maps to virus
    $splitreads{$read} = '1to2h' if ($ref2 eq '=' and $flag & 64); ##this read is first in pair and its mate maps to host
    $splitreads{$read} = '2to1h' if ($ref2 eq '=' and $flag & 128); ##this read is second in pair and its mate maps to host
    
    ##double check to see if we've determined the reads mate orientation
    print STDERR "\nWarning: unknown read pair configuration for read $read at alignment $j, possible incorrect filtering?:\nLine: $line\n" if ($splitreads{$read} == 0);
    $e++  if ($splitreads{$read} == 0);
    
    ##where does split read piece go on virus genome? onBorder, nearBorder, or center?
    my $read_virus_pos = CheckVirusPos($pos3,$len3,$cutend3,$virus_len,$virus_end_len);
    die "\nError: center clipped read can't have its virus position checked.\n" if (!$read_virus_pos); ##this shouldn't ever happen, but...
    
    ##get coordinates of cut end (5bp) of read
    ##also test for concordant or canonical read pair and read split orientation and alignment relative to virus insertion site
    my ($endid,$sense,$concordant) = GetVirusEnd($ref1,$pos1,$len1,$cutend1,$overlap_len,$read_virus_pos,$splitreads{$read});
    
    ##get 5bp overlap region sequence
     if (!exists $seqs{$endid}){
      $seqs{$endid} = GetSeq($endid,\%contig_seqs);
      die "\nError: could not parse overlap sequence for overlap: $endid; line: $i\n\n" if (!$seqs{$endid});
     }
    
    #bin virus ends
    $total_splitends{$endid} += 1;
    $total_center_splitends{$endid} += 1 if ($read_virus_pos eq 'center');
    $rightend_splitends{$endid} += 1 if ($cutend1 eq 'right');
    $leftend_splitends{$endid} += 1 if ($cutend1 eq 'left');
    if ($concordant){
      $total_con_splitends{$endid}{$sense} += 1;
      $total_con_center_splitends{$endid}{$sense} += 1 if ($read_virus_pos eq 'center');
      $rightend_con_splitends{$endid}{$sense} += 1 if ($cutend1 eq 'right');
      $leftend_con_splitends{$endid}{$sense} += 1 if ($cutend1 eq 'left');
      $y++;
    }
    $x++;
  }
}

##print out overlaps
print $outfh "Chromosome\tOverlapRange\tOverlapSequence\tTotalSplits\tTotalSplits_VirusCenter\tRightSplits\tLeftSplits\tTotalConcordantSplits\tTotalConcordantSplits_VirusCenter\t";
print $outfh "TotalConcordantSenseVirusSplits\tTotalConcordantAntisenseVirusSplits\tRightConcordantSplits\tLeftConcordantSplits\n";
foreach my $endid (sort keys %total_splitends){
  if ($total_splitends{$endid} >= $min_reads){
    my ($chr,$start,$end) = (split /:/,$endid);
    print $outfh "$chr\t$start-$end\t$seqs{$endid}\t$total_splitends{$endid}";
    if (exists $total_center_splitends{$endid}){
      print $outfh "\t$total_center_splitends{$endid}";
    }
    else{
      print $outfh "\t0";
    }
    if (exists $rightend_splitends{$endid}){
      print $outfh "\t$rightend_splitends{$endid}";
    }
    else{
      print $outfh "\t0";
    }
    if (exists $leftend_splitends{$endid}){
      print $outfh "\t$leftend_splitends{$endid}";
    }
    else{
      print $outfh "\t0";
    }
    if (exists $total_con_splitends{$endid}{'plus'} or exists $total_con_splitends{$endid}{'minus'}){
      my $total = $total_con_splitends{$endid}{'plus'} + $total_con_splitends{$endid}{'minus'};
      print $outfh "\t$total";
    }
    else{
      print $outfh "\t0";
    }
    if (exists $total_con_center_splitends{$endid}{'plus'} or exists $total_con_center_splitends{$endid}{'minus'}){
      my $total = $total_con_center_splitends{$endid}{'plus'} + $total_con_center_splitends{$endid}{'minus'};
      print $outfh "\t$total";
    }
    else{
      print $outfh "\t0";
    }
    if (exists $total_con_splitends{$endid}{'plus'}){
      print $outfh "\t$total_con_splitends{$endid}{'plus'}";
    }
    else{
      print $outfh "\t0";
    }
    if (exists $total_con_splitends{$endid}{'minus'}){
      print $outfh "\t$total_con_splitends{$endid}{'minus'}";
    }
    else{
      print $outfh "\t0";
    }
    if (exists $rightend_con_splitends{$endid}{'plus'} or exists $rightend_con_splitends{$endid}{'minus'}){
      my $total = $rightend_con_splitends{$endid}{'plus'} + $rightend_con_splitends{$endid}{'minus'};
      print $outfh "\t$total";
    }
    else{
      print $outfh "\t0";
    }
    if (exists $leftend_con_splitends{$endid}{'plus'} or exists $leftend_con_splitends{$endid}{'minus'}){
      my $total = $leftend_con_splitends{$endid}{'plus'} + $leftend_con_splitends{$endid}{'minus'};
      print $outfh "\t$total";
    }
    else{
      print $outfh "\t0";
    }
    print $outfh "\n";
  }
}


print STDERR "\nTotal lines read: $i\n";
print STDERR "Total alignments: $j\n";
print STDERR "Total split reads (virus/host): $x\n";
print STDERR "Total concordantly aligned split reads: $y\n";
print STDERR "Total reads with unknown pair orientation: $e\n";
print STDERR "Reads skipped due to split complexity: $z\n";
print STDERR "Reads skipped due to mate overlap (duplicate splitting): $m\n";
print STDERR "Reads skipped due to double-sided clipping: $n\n";
print STDERR "Reads skipped due to deletion: $d\n";
print STDERR "Reads skipped due to insertion: $f\n";
print STDERR "Reads skipped due to bit error: $c\n\n";

##get 5bp overlap sequence
sub GetSeq{
  my ($seqid,$seqs) = @_;
  my ($chr,$start,$stop) = (split /:/,$seqid);
  if ($seqs->{$chr}){
    my $seq = substr $seqs->{$chr},($start-1),($stop-($start-1));
    return 0 if (length($seq) != 5);
    return $seq;
  }
  return 0;
}

##get virus overlap region info and test for concordance
##assumes filtering was done correctly
sub GetVirusEnd{
  my ($ref,$pos,$len,$cutend,$overlaplen,$readviruspos,$readpair) = @_;
  my ($endid,$concordance,$overlap_start,$overlap_end,$sense);
  $overlap_start = ($pos+$len) - $overlaplen if ($cutend eq 'right');
  $overlap_start = $pos if ($cutend eq 'left');
  $overlap_end = $pos+$len-1 if ($cutend eq 'right');
  $overlap_end = $pos+$overlaplen-1 if ($cutend eq 'left');
  
  ##test for virus insertion orientation and split read concordance
  ##most standard read configuration around virus insertion:
  ##  --1st read = sense & 2nd = antisense;
  ##  --1st read is split with 'left' chunk on host & 'right' chunk on virus;
  ##  --'left' chunk is split on the right side
  ##  --'right' chunk is on 'left' side of virus sequence with cut end facing to the left (i.e. virus is sense orientation to host)
  ##  --2nd read is on virus
  if ($readpair eq '1to2v' and $readviruspos eq 'onLeftBorder1' and $cutend eq 'right'){
    $endid = "$ref:$overlap_start:$overlap_end";
    $sense = 'plus';
    $concordance = 1;
  }
  ##same as above but virus is antisense relative to host
  elsif ($readpair eq '1to2v' and $readviruspos eq 'onRightBorder1' and $cutend eq 'right'){
    $endid = "$ref:$overlap_start:$overlap_end";
    $sense = 'minus';
    $concordance = 1;
  }
  ##second read is on host & virus is sense orientation
  elsif ($readpair eq '1to2h' and $readviruspos eq 'onRightBorder1' and $cutend eq 'left'){
    $endid = "$ref:$overlap_start:$overlap_end";
    $sense = 'plus';
    $concordance = 1;
  }
  ##second read is on host & virus is antisense orientation
  elsif ($readpair eq '1to2h' and $readviruspos eq 'onLeftBorder1' and $cutend eq 'left'){
    $endid = "$ref:$overlap_start:$overlap_end";
    $sense = 'minus';
    $concordance = 1;
  }
  ##first read is on virus (second read is split) & virus is sense orientation
  elsif ($readpair eq '2to1v' and $readviruspos eq 'onRightBorder1' and $cutend eq 'left'){
    $endid = "$ref:$overlap_start:$overlap_end";
    $sense = 'plus';
    $concordance = 1;
  }
  ##first read is on virus & virus is antisense
  elsif ($readpair eq '2to1v' and $readviruspos eq 'onLeftBorder1' and $cutend eq 'left'){
    $endid = "$ref:$overlap_start:$overlap_end";
    $sense = 'minus';
    $concordance = 1;
  }
  ##first read is on host & virus is sense
  elsif ($readpair eq '2to1h' and $readviruspos eq 'onRightBorder1' and $cutend eq 'right'){
    $endid = "$ref:$overlap_start:$overlap_end";
    $sense = 'plus';
    $concordance = 1;
  }
  ##first read is on host & virus is antisense
  elsif ($readpair eq '2to1h' and $readviruspos eq 'onLeftBorder1' and $cutend eq 'right'){
    $endid = "$ref:$overlap_start:$overlap_end";
    $sense = 'minus';
    $concordance = 1;
  }
  else{
    $endid = "$ref:$overlap_start:$overlap_end";
    $sense = 'na';
    $concordance = 0;
  }
  return ($endid,$sense,$concordance);
}

##check for split read tag and return parsed alignment info
sub CheckSplit{
  my (@tags) = @_;
  foreach my $tag (@tags){
    if ($tag =~ m/^SA/){
      my @SA = split /;/,$tag;
      if ($SA[0] =~ m/SA:Z:([^,]+)\,([0-9]+)\,[+-]\,([0-9SHMDI]+)/){
        my ($chr,$start,$cigar) = ($1,$2,$3);
        my ($len,$cutend);
        ($start,$len,$cutend) = GetMatchPart($start,$cigar);
        return ($chr,$start,$len,$cutend);
      }
    }
  }
  return 0;
}

##parse SAM file CIGAR string and return specified part (either M, S, or H), AND modified read start-stop positions
sub GetMatchPart{
  my ($start,$cigar) = @_;
  my ($sta,$len,$cutend);
  if ($cigar =~ m/^([0-9]+)[HS]([0-9]+)M([0-9]+)[HS]$/){
    $sta = $start;
    $len = $2;
    $cutend = 0;
    return ($sta,$len,$cutend);
  }
  elsif ($cigar =~ m/^([0-9]+)[HS]([0-9]+)M$/){
    $sta = $start;
    $len = $2;
    $cutend = 'left';
    return ($sta,$len,$cutend);
  }
  elsif ($cigar =~ m/^([0-9]+)M([0-9]+)[HS]$/){
    $sta = $start;
    $len = $1;
    $cutend = 'right';
    return ($sta,$len,$cutend);
  }
  elsif ($cigar =~ m/D/){
    return ('d',0,0);
  }
  elsif ($cigar =~ m/I/){
    return ('i',0,0);
  }
  return 0;
}

sub CheckVirusPos{
  my ($pos,$len,$cutend,$virus_len,$virus_end_len) = @_;
  
  ##check if read start is on or very close to virus/host border
  if ($pos <= 5){
    return "onLeftBorder0" if ($cutend eq 'right');
    return "onLeftBorder1" if ($cutend eq 'left');
    return 0;
  }
  elsif (($pos+$len) >= ($virus_len-5)){
    return "onRightBorder0" if ($cutend eq 'left');
    return "onRightBorder1" if ($cutend eq 'right');
    return 0;
  }
  ##check if read start is close to virus/host border (but not on it)
  elsif ($pos <= $virus_end_len){
    return "nearLeftBorder0" if ($cutend eq 'right');
    return "nearLeftBorder1" if ($cutend eq 'left');
    return 0;
  }
  elsif (($pos+$len) >= ($virus_len-$virus_end_len)){
    return "nearRightBorder0" if ($cutend eq 'left');
    return "nearRightBorder1" if ($cutend eq 'right');
    return 0;
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
