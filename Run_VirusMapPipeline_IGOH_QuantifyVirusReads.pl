#!/usr/bin/perl

#use strict;
use Cwd;
use File::Spec;
use File::Basename;

# Code by:      Cris Vera
# Group:        IGOH Theme, IGB, University of Illinois
# Date:         06-24-2019
# Contact:      jcvera@illinois.edu

##description: runs pipeline to quantify virus reads within isolate or deep coverage samples. Uses an unmodifed virus genome reference with either
##unmodifed or mod2 modified host genome reference (i.e. mod2=remove prophage sequence if it exists).
###++Runs the main quantification script (EditSAMs_CountVirusPositions) directly on SAM file of each sample for generation of broad overview sample figures
###++Also runs a virus locator script on all the samples (summary file output) to locate virus insertion sites within 500bp bins for multiple
###  isolate samples. This can be run on deep coverage samples as well, but will return less useful results, since it's intended to detect virus
###  insertion site movement between isolate samples.
###++Finally, runs the overlap detection script on each sample to list out exact virus insertion locations (for mu-like phages where there's an
###  insertion overlap region). This is intended for use on deep coverage samples, but it can be run on isolates with a mu-like phage as well.
##Notes:
####-to prep a Host/Virus DB:
####-1-if the Host reference contains the prophage, remove it using the PrepGenome script.  I tag references with modified host genomes with '_mod2'
####   e.g. UnFastaFormat.pl -i HostGenomeFile.fasta | PrepGenome.pl -r 1:5000:25000 | FastaFormat.pl -o HostGenomeFile.mod2.fasta
####-2-concatenate the Host reference with the Phage reference
####   e.g. cat HostGenomeFile.mod2.fasta PhageGenomeFile.fasta >HostAndPhageGenomeFile.fasta
####-3-run BWA index
####   e.g. bwa index HostAndPhageGenomeFile.fasta
####-4-add to this script (sorry, this requires too many manual decisions to easily automate)

###CONFIG FILE NAME###
my $config_file = 'WGSPipeline_REF_CONFIG.tsv';
###

###CONFIG FILE PATH###
my $config_dir = Cwd::abs_path(dirname(__FILE__));
###

##version
my $version = v1.0;

# defaults, initializations, and constants
my $help = "\n\nRun_VirusMapPipeline_QuantifyVirusReads $version. Runs pipeline to quantify virus reads within isolate or deep coverage samples.\n".
            "\t-d  Option: specify input directory.  Required.\n".
            "\t-e  Option: specify output directory. Default='out/'.\n".
            "\t-t  Option: specify number of threads. Default=8.\n".
            "\t-p  Option: specify project name (summary file base). Default='PA'.\n".
            "\t-b  Option: specify quantification bin size. Default=1000.\n".
            "\t-f  Option: specify minimum overlap coverage. Default=2.\n".
            "\t-V  Option: specify Hybrid Host/Virus DB to use. Default=1.\n";
my $usage = "\nRun_VirusMapPipeline_QuantifyVirusReads.pl -d [Input Dir] -e [Output Dir] -p [Project Name] -v [Host/Virus BWA DB] -b [Bin Size] -f [Min Overlap Coverage]\n";
my ($indir,$outdir,$project);
my ($refdir,$scriptdir) = GetDirs("$config_dir/$config_file");
my $date = GetDate();
my $i = my $j = 0;
my $virus_output = 'VirusPipeline_QuantifyVirus';
my $linkdir = 'quantificationlinks';
my $procs = 8;
my $ref = 1;
my $bordersize = 1200;  ##virus modification border size
my $binsize = 1000;     ##quantification bin size
my $minreads = 2;       ##minimum split read overlap coverage


#process command line custom script tags
my %cmds = ReturnCmds(@ARGV);
if ($cmds{'h'}){
  $help .= GetRefList("$config_dir/$config_file");
  $help .= "\n************\nAuthor:\t\tCris Vera\nGroup:\t\tIGOH Theme, IGB, University of Illinois\nContact:\tjcvera\@illinois.edu\n";  
  die "\n$help\n$usage\n";
}
$indir = $cmds{'d'} if ($cmds{'d'});
$outdir = $cmds{'e'} if ($cmds{'e'});
$ref = $cmds{'V'} if ($cmds{'V'});
$project = $cmds{'p'} if ($cmds{'p'});
$binsize = $cmds{'b'} if ($cmds{'b'});
$minreads = $cmds{'f'} if ($cmds{'f'});
die "\nError: must specify input directory.\n\n" if (!$indir);
die "\nError: must specify output directory.\n\n" if (!$outdir);
die "\nError: must specify project base.\n\n" if (!$project);

###tool directories###
my %tools = GetConfigs('TOOL',"$config_dir/$config_file");
my $sambamba = $tools{'Sambamba'};
my $qualimap = $tools{'Qualimap'};
my $bwa = "$tools{'BWA'} mem";

###script directories###
my $quant = "perl $scriptdir/EditSAMs_CountVirusPositions.pl";
my $overlap = "perl $scriptdir/EditSAMs_DetectOverlaps.pl";
my $prep_virus = "perl $scriptdir/PrepGenome.pl";
my $unfasta = "perl $scriptdir/UnFastaFormat.pl";
my $fasta = "perl $scriptdir/FastaFormat.pl";
my $combine = "perl $scriptdir/CombineSAMs.pl";
my $table = "perl $scriptdir/SummarizeLocationResults_VirusMapPipeline_QuantifyVirusReads.pl";

###reference genomes###
my %refs = GetConfigs('HREF',"$config_dir/$config_file");
my $hostvirus_db = my $hostvirus_db_txt = $refs{$ref}{'BWA'};
my $refdes = $refs{$ref}{'Refname'};
my $refname = $refs{$ref}{'Ref'};
my $virus_header = $refs{$ref}{'VirusHeader'};
$hostvirus_db_txt .= '.txt';
$virus_output .= "_$refname";
$linkdir .= "_$refname";

##make absolute paths
$indir = File::Spec->rel2abs($indir);
$outdir = File::Spec->rel2abs($outdir);

##logfile start
my $logfile = "$project\_VirusPipeline.log";
if (!-e "$outdir/$logfile"){
  open ('LOG', ">$outdir/$logfile") or die "Cannot create $outdir/$logfile: $!\n";
  print LOG "Run date: $date\n";
  print LOG "Threads: $procs\n";
  print LOG "\nProject name: $project\n";
  print LOG "Input Dir: $indir\n";
  print LOG "Output Dir: $outdir\n";
  print LOG "\nHybrid reference genome: $refdes\n";
}
else{
  open ('LOG', ">>$outdir/$logfile") or die "Cannot write to $outdir/$logfile: $!\n";
  print LOG "\n\n**************\nRe-Running Virus Pipeline\n**************\n";
  print LOG "Run date: $date\n";
  print LOG "Threads: $procs\n";
  print LOG "\nProject name: $project\n";
  print LOG "Input Dir: $indir\n";
  print LOG "Output Dir: $outdir\n";
  print LOG "\nHybrid reference genome: $refdes\n";
}

#get all sample file names to run script on
opendir(INDIR,$indir) or die "Can't open directory: $indir: $!\n";
my @samplefiles = grep {m/(_R1_001\.fastq$)|(_R1\.fastq$)/ && -f "$indir/$_"} readdir(INDIR);
close(INDIR);

##create link file folder
 if (!-e "$outdir/$linkdir"){
  mkdir "$outdir/$linkdir";
 }
else{
  die "\nError: link directory already exists: $outdir/$linkdir\n\n";
}

###create directories, and run pipeline commands
foreach my $samplefile1 (@samplefiles){
  $i += 1;
  my $samplefile2 = my $samplename = $samplefile1;
  $samplename =~ s/(_R1_001\.fastq)|(_R1\.fastq)//;
  $samplefile2 =~ s/_R1_001\.fastq/_R2_001.fastq/;
  $samplefile2 =~ s/_R1\.fastq/_R2.fastq/;
  my $sample = $samplename;
  $sample =~ s/_[ATGC-]+_(L00[0-9])$/_$1/;
  $sample =~ s/-pseudomonas-[ATGC-]+$//;  ##some raw read files had '-pseudomonas' added to them for some reason
  ##check for dir and run
  if (-e "$outdir/Sample_$samplename"){
    if (!-e "$outdir/Sample_$samplename/$virus_output"){
      print STDERR "\nRunning sample $i: $sample\n";
      print LOG "\nRunning sample $i: $sample\n";
      mkdir "$outdir/Sample_$samplename/$virus_output";
      
      ##main BWA alignment vs host/virus -- switch to using original BAM file in sample folder
      #my $bwa_cmd = "$bwa mem -t $procs $hostvirus_db $outdir/Sample_$samplename/$samplename\_R1_trimmed.fastq $outdir/Sample_$samplename/$samplename\_R2_trimmed.fastq 2>$outdir/Sample_$samplename/$virus_output/bwa_mem.err | $sambamba view -S -f bam /dev/stdin >$outdir/Sample_$samplename/$virus_output/$sample\_vsHostVirus_BWA-mem.bam 2>$outdir/Sample_$samplename/$virus_output/sambamba_view.err";
      #my $stats_cmd = "$bamtools stats -in $outdir/Sample_$samplename/$virus_output/$sample\_vsHostVirus_BWA-mem.bam >$outdir/Sample_$samplename/$virus_output/$sample\_vsHostVirus_BWA-mem.bam_stats 2>$outdir/Sample_$samplename/$virus_output/$sample\_vsHostVirus_BWA-mem.bam_stats.err";
      
      ##starting BAM file: $outdir/Sample_$samplename/$sample\_BWA-mem.RG.sorted.markedup.bam
      ##sort and run qualimap for coverage stats
      #my $sortfilter_cmd = "$sambamba sort -F \"not unmapped and [XA] == null\" -t $procs -m 128G --tmpdir $outdir/Sample_$samplename/temp -o $outdir/Sample_$samplename/$virus_output/$sample\_vsHostVirus_BWA-mem.sortfilter1.bam $outdir/Sample_$samplename/$virus_output/$sample\_vsHostVirus_BWA-mem.bam 2>$outdir/Sample_$samplename/$virus_output/sambamba_sortfilter1.err";
      my $sortfilter_cmd = "$sambamba sort -F \"not unmapped and [XA] == null\" -t $procs -m 128G --tmpdir $outdir/Sample_$samplename/temp -o $outdir/Sample_$samplename/$virus_output/$sample\_vsHostVirus_BWA-mem.sortfilter1.bam $outdir/Sample_$samplename/$sample\_BWA-mem.RG.sorted.markedup.bam 2>$outdir/Sample_$samplename/$virus_output/sambamba_sortfilter1.err";
      my $qualimap_cmd = "$qualimap bamqc -bam $outdir/Sample_$samplename/$virus_output/$sample\_vsHostVirus_BWA-mem.sortfilter1.bam -ip -c -outdir $outdir/Sample_$samplename/$virus_output/qualimap -outformat PDF -nt $procs --java-mem-size=128G -nw 500 -p NON-STRAND-SPECIFIC 2>$outdir/Sample_$samplename/$virus_output/$sample\_qualimap.err";
      
      ##sort for virus split reads and virus read pairs
      my $sortfilter_virussplit_cmd = "$sambamba sort -F \"ref_name =~ /$virus_header/ and not [SA] == null and not [SA] =~ /$virus_header/\" -t $procs -m 128G --tmpdir $outdir/Sample_$samplename/temp -o $outdir/Sample_$samplename/$virus_output/$sample\_vsHostVirus_BWA-mem.sortfilter_virusSplits.bam $outdir/Sample_$samplename/$virus_output/$sample\_vsHostVirus_BWA-mem.sortfilter1.bam 2>$outdir/Sample_$samplename/$virus_output/sambamba_sortfilter_virusSplits.err";
      my $sortfilter_viruspairs_cmd = "$sambamba sort -F \"ref_name =~ /$virus_header/ and not mate_ref_name =~ /$virus_header/ and [SA] == null\" -t $procs -m 128G --tmpdir $outdir/Sample_$samplename/temp -o $outdir/Sample_$samplename/$virus_output/$sample\_vsHostVirus_BWA-mem.sortfilter_virusPairs.bam $outdir/Sample_$samplename/$virus_output/$sample\_vsHostVirus_BWA-mem.sortfilter1.bam 2>$outdir/Sample_$samplename/$virus_output/sambamba_sortfilter_virusPairs.err";
      ##sort for host split reads
      my $sortfilter_hostsplit_cmd = "$sambamba sort -F \"not ref_name =~ /$virus_header/ and [SA] =~ /$virus_header/\" -t $procs -m 128G --tmpdir $outdir/Sample_$samplename/temp -o $outdir/Sample_$samplename/$virus_output/$sample\_vsHostVirus_BWA-mem.sortfilter_hostSplits.bam $outdir/Sample_$samplename/$virus_output/$sample\_vsHostVirus_BWA-mem.sortfilter1.bam 2>$outdir/Sample_$samplename/$virus_output/sambamba_sortfilter_hostSplits.err";
      
      ##make sam files
      my $sam_cmd1 = "$sambamba view -h $outdir/Sample_$samplename/$virus_output/$sample\_vsHostVirus_BWA-mem.sortfilter_virusSplits.bam >$outdir/Sample_$samplename/$virus_output/$sample\_vsHostVirus_BWA-mem.sortfilter_virusSplits.sam 2>>$outdir/Sample_$samplename/$virus_output/sam_view.err";
      my $sam_cmd2 = "$sambamba view -h $outdir/Sample_$samplename/$virus_output/$sample\_vsHostVirus_BWA-mem.sortfilter_virusPairs.bam >$outdir/Sample_$samplename/$virus_output/$sample\_vsHostVirus_BWA-mem.sortfilter_virusPairs.sam 2>>$outdir/Sample_$samplename/$virus_output/sam_view.err";
      my $sam_cmd3 = "$sambamba view -h $outdir/Sample_$samplename/$virus_output/$sample\_vsHostVirus_BWA-mem.sortfilter_hostSplits.bam >$outdir/Sample_$samplename/$virus_output/$sample\_vsHostVirus_BWA-mem.sortfilter_hostSplits.sam 2>>$outdir/Sample_$samplename/$virus_output/sam_view.err";
      
      ##combine virus sam files
      my $combine_cmd = "$combine -i $outdir/Sample_$samplename/$virus_output/$sample\_vsHostVirus_BWA-mem.sortfilter_virusSplits.sam -j $outdir/Sample_$samplename/$virus_output/$sample\_vsHostVirus_BWA-mem.sortfilter_virusPairs.sam -o $outdir/Sample_$samplename/$virus_output/$sample\_vsHostVirus_BWA-mem.sortfilter_virusSplitPairs.sam 2>$outdir/Sample_$samplename/$virus_output/combineSAMs.err";
      
      ##quantification command
      my $quant_cmd = "$quant -i $outdir/Sample_$samplename/$virus_output/$sample\_vsHostVirus_BWA-mem.sortfilter_virusSplitPairs.sam -o $outdir/Sample_$samplename/$virus_output/$sample\_vsHostVirus_BWA-mem.sortfilter_virusSplitPairs_Quantification.tsv -c $virus_header -b $binsize 2>$outdir/Sample_$samplename/$virus_output/quantification.err";
      
      ##overlap detection
      my $overlap_cmd;
      if (!-e "$hostvirus_db_txt"){
       $overlap_cmd = "$unfasta -i $hostvirus_db -o $hostvirus_db_txt 2>$outdir/Sample_$samplename/$virus_output/unfasta.err && $overlap -i $outdir/Sample_$samplename/$virus_output/$sample\_vsHostVirus_BWA-mem.sortfilter_hostSplits.sam -o $outdir/Sample_$samplename/$virus_output/$sample\_vsHostVirus_BWA-mem.sortfilter_hostSplits_Overlaps.tsv -c $virus_header -f $minreads -g $hostvirus_db_txt 2>$outdir/Sample_$samplename/$virus_output/overlap.err";
      }
      else{
       $overlap_cmd = "$overlap -i $outdir/Sample_$samplename/$virus_output/$sample\_vsHostVirus_BWA-mem.sortfilter_hostSplits.sam -o $outdir/Sample_$samplename/$virus_output/$sample\_vsHostVirus_BWA-mem.sortfilter_hostSplits_Overlaps.tsv -c $virus_header -f $minreads -g $hostvirus_db_txt 2>$outdir/Sample_$samplename/$virus_output/overlap.err";
      }
      
      ##print to log and run all commands
      print LOG "Original BWA alignment BAM file: $outdir/Sample_$samplename/$sample\_BWA-mem.RG.sorted.markedup.bam\n";
      print LOG "First Sambamba filter (unmapped & multimapped): $sortfilter_cmd\n";
      
      my $run_cmd = `$sortfilter_cmd`;
      print LOG "Qualimap: $qualimap_cmd\n";
      $run_cmd = `$qualimap_cmd`;
      sleep 15;  ##puase for qulimap to finish
      print LOG "Sambamba filter for Virus mapped split reads: $sortfilter_virussplit_cmd\n";
      $run_cmd = `$sortfilter_virussplit_cmd`;
      print LOG "Sambamba filter for Virus mapped discordant reads: $sortfilter_viruspairs_cmd\n";
      $run_cmd = `$sortfilter_viruspairs_cmd`;
      print LOG "Sambamba filter for Host mapped split reads: $sortfilter_hostsplit_cmd\n";
      $run_cmd = `$sortfilter_hostsplit_cmd`;
      print LOG "Sambamba convert 2 SAM 1: $sam_cmd1\n";
      $run_cmd = `$sam_cmd1`;
      print LOG "Sambamba convert 2 SAM 2: $sam_cmd2\n";
      $run_cmd =  `$sam_cmd2`;
      print LOG "Sambamba convert 2 SAM 3: $sam_cmd3\n";
      $run_cmd = `$sam_cmd3`;
      print LOG "Combine SAMs: $combine_cmd\n";
      $run_cmd = `$combine_cmd`;
      print LOG "Quantification of Virus mapped splits & pairs: $quant_cmd\n";
      $run_cmd =  `$quant_cmd`;
      print LOG "Find Overlaps for Host mapped splits: $overlap_cmd\n";
      $run_cmd =  `$overlap_cmd`;
      $j++;
    }
    else{
      print LOG "skipping sample $sample. It already exists.\n\n";
      print STDERR "skipping sample $sample. It already exists.\n\n";
    }
    ##make links
    print LOG "\nLink commands:\n1: ln -s -T $outdir/Sample_$samplename/$virus_output/qualimap/genome_results.txt $outdir/$linkdir/$sample\_genome_results.txt\n";
    print LOG "2: ln -s -T $outdir/Sample_$samplename/$virus_output/$sample\_vsHostVirus_BWA-mem.sortfilter_virusSplitPairs_Quantification.tsv $outdir/$linkdir/$sample\_virusSplitPairs_Quantification.tsv\n";
    print LOG "3: ln -s -T $outdir/Sample_$samplename/$virus_output/$sample\_vsHostVirus_BWA-mem.sortfilter_hostSplits_Overlaps.tsv $outdir/$linkdir/$sample\_hostSplits_Overlaps.tsv\n\n";
    my $make_link = `ln -s -T $outdir/Sample_$samplename/$virus_output/qualimap/genome_results.txt $outdir/$linkdir/$sample\_genome_results.txt 2>&1`;
    print STDERR "\nWarning: could not create link: $outdir/$linkdir/$sample\_genome_results.txt: $make_link\n\n" if ($make_link);
    $make_link = `ln -s -T $outdir/Sample_$samplename/$virus_output/$sample\_vsHostVirus_BWA-mem.sortfilter_virusSplitPairs_Quantification.tsv $outdir/$linkdir/$sample\_virusSplitPairs_Quantification.tsv 2>&1`;
    print STDERR "\nWarning: could not create link: $outdir/$linkdir/$sample\_virusSplitPairs_Quantification.tsv: $make_link\n\n" if ($make_link);
    $make_link = `ln -s -T $outdir/Sample_$samplename/$virus_output/$sample\_vsHostVirus_BWA-mem.sortfilter_hostSplits_Overlaps.tsv $outdir/$linkdir/$sample\_hostSplits_Overlaps.tsv 2>&1`;
    print STDERR "\nWarning: could not create link: $outdir/Sample_$samplename/$virus_output/$sample\_vsHostVirus_BWA-mem.sortfilter_hostSplits_Overlaps.tsv: $make_link\n\n" if ($make_link);
  }
  else{
    die "\nError: couldn't find sample dir: $outdir/Sample_$samplename\n\n";
  }
}
print STDERR "\nTotal samples found: $i\nSamples quantified: $j\n\n";
print LOG "\nTotal samples found: $i\nSamples quantified: $j\n\n";

##run summary table
my $table_cmd = "$table -d $outdir/$linkdir -o $outdir/$project\_$virus_output\_QuantLocationSummary.tsv -v $virus_header 2>$outdir/$virus_output\_QuantLocationSummary.err";
print LOG "\nGenerate Summary table: $table_cmd\n\n";
my $run = `$table_cmd`;

sub GetHeader{
  my ($unfasta,$file) = @_;
  my $run = `$unfasta -i $file -o TempFile.tmp 2>&1`;
  print STDERR "\nWarning sub GetHeader: error detected during UnFasta: $run\n\n" if ($run =~ m/error/i);
  open (TMP, "<TempFile.tmp") or die "Cannot open TempFile.tmp: $!\n";
  chomp(my @lines = <TMP>);
  close(TMP);
  my @line = split /\t/,$lines[0];
  return $line[0];
}


##get all config values sorted into integers for selection
sub GetConfigs{
  my ($reftype,$configfile) = @_;
  my %refs;
  my %refs2;
  open (CON, "<$configfile") or die "Cannot open $configfile: $!\n";
  while (my $line = <CON>){
    chomp $line;
    #if ($line =~ m/^REF\_([^\t]+)\t([^\t]+)\t(.+)/){
    #  $refs{$1}{$2} = $3;
    #}
    if ($line =~ m/^HREF\_([^\t]+)\t([^\t]+)\t(.+)/){
      $refs{$1}{$2} = $3;
    }
    elsif ($line =~ m/^$reftype\_([^\t]+)\t(.+)/){
      $refs{$1} = $2;
    }
  }
  close(CON);
  my $n = 0;
  if ($reftype eq 'HREF'){
    foreach my $ref (sort keys %refs){
      $n++;
      $refs2{$n}{'Ref'} = $ref;
      foreach my $type (keys %{$refs{$ref}}){
        $refs2{$n}{$type} = $refs{$ref}{$type};
      }
    }
    return %refs2;
  }
  else{
    return %refs;
  }
  return 0;
}

sub GetRefList{
  my ($config_file) = @_;
  my %refs;
  my $refs;
  open (CON, "<$config_file") or die "Cannot open $config_file: $!\n";
  while (my $line = <CON>){
    chomp $line;
    #if ($line =~ m/^REF_([^\t]+)\t([^\t]+)\t([^\t]+)/){
    #  $refs{$1} += 1;
    #}
    if ($line =~ m/^HREF_([^\t]+)\t([^\t]+)\t([^\t]+)/){
      $refs{$1} += 1;
    }
  }
  close(CON);
  my $n = 0;
  foreach my $ref (sort keys %refs){
    $n++;
    $refs .= "\t\t$n: $ref ($refs{$ref})\n";
  }
  return $refs;
}

sub GetDirs{
  my ($config_file) = @_;
  my ($scriptdir,$refdir);
  open (CON, "<$config_file") or die "Cannot open $config_file: $!\n";
  while (my $line = <CON>){
    chomp $line;
    if ($line =~ m/^#=REF_DIR\t(.+)/){
      $refdir = $1;
    }
    elsif ($line =~ m/^#=SCRIPT_DIR\t(.+)/){
      $scriptdir = $1;
    }
  }
  close(CON);
  return ($refdir,$scriptdir);
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
sub GetDate{
  my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
  my @months = qw(January February March April May June July August September October November December);
  $year += 1900;
  my $date = "$months[$mon] $mday, $year";
  return $date;
}

