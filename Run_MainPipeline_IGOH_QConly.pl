#!/usr/bin/perl

use strict;
use Cwd;
use File::Spec;
use File::Basename;

# Code by:      Cris Vera
# Group:        IGOH Theme, IGB, University of Illinois
# Date:         06-24-2019
# Contact:      jcvera@illinois.edu

##description: Runs the Whitaker lab WGS pipeline through initial QC and alignment steps. Establishes a directory hiearchy for output and downstream pipelines.

###CONFIG FILE NAME###
my $config_file = 'WGSPipeline_REF_CONFIG.tsv';
###

###CONFIG FILE PATH###
my $config_dir = Cwd::abs_path(dirname(__FILE__));
###

##version
my $version = v1.0;

####defaults, initializations, and constants####
my $help = "\nRun_PseudomonasPipeline_IGOH_QConly $version\nDescription: Runs Whitaker lab WGS pipeline through initial QC and alignment steps. Establishes a directory hiearchy for output and downstream pipelines.\n".          
          "\n  ###Main Options-Switches:\n".
          "\t-d  Option: specify input directory with sample files. Required.\n".
          "\t-e  Option: specify output directory. Required.\n".
          "\t-p  Option: specify project name. Required.\n".
          "\t-t  Option: specify number of threads to use (when possible). Default=8.\n".
          "\t-z  Option: raw reads are gzipped. Default=0.\n".
          "\t-u  Option: run assembly step. Defualt=0.\n".
          "\t\t 0: no assembly\n".
          "\t\t 1: Unicycler assembly, use raw reads with correction\n".
          "\t\t 2: Unicycler assembly, use trim reads with no correction\n".
          "\t\t 3: Spades assembly, use raw reads with correction\n".
          "\t\t 4: Spades assembly, use trim reads with no correction\n".
          "\t-m  Option: run ariba MLST. Default=0.\n".
          "\t\t 0: no Ariba\n".
          "\t\t 1: Ariba MLST\n".
          "\t\t 2: Ariba MLST, Magares, and VFDB\n".
          "\t-k  Option: specify sequencing library name. Default=Unknown.\n".
          "\t\t 0: Unknown\n".
          "\t\t 1: Nextera_XT_DNA_Library_Prep_Kit\n".
          "\t\t 2: TruSeq_RNA_Library_Prep_Kit\n".
          "\t\t 3: TruSeq_Stranded_mRNA_LT_Kit\n".
          "\t\t 4: TruSeq_Stranded_Total_RNA_LT_Kit\n".
          "\t\t 5: TruSeq_Small_RNA_Library_Prep_Kit\n".
          "\t\t 6: TruSeq_Chip_Library_Prep_Kit\n".
          "\t\t 7: TruSeq_DNA_Library_Prep_Kit\n".
          "\t\t 8: NEBNext_Multiplex_Small_RNA_Library_Prep_Kit\n".
          "\t\t 9: Clonetech_SMARTer_Ultra_Low_Input_RNA_Kit\n".
          "\t\t10: SureSelect_XT_Mouse_All_Exome\n".
          "\t\t11: SureSelect_XT_Human_All_Exon_V5-UTR\n".
          "\t\t12: iGenomx_Riptide_DNA_HT-RLP\n".
          "\t\t13: Kappa_DNA_HTP_Prep_Kit\n".
          "\t\t14: Kappa_PCR-free_DNA_HTP_Prep_Kit\n".
          "\t\t15: Nextera_Flex_DNA_Library_Prep_Kit\n".
          "\t-R  Option: specify reference genome/assembly. Default=1.\n";
my $usage = "\nUsage:\nRun_PseudomonasPipeline_IGOH_QConly.pl -d [Input Sample Directory] -e [Output Directory] -p [Project Name] -t [Threads] -z [Gzipped] -m [Ariba] -R [Reference] -u [Assembly] -k [Library Kit]\n";
my ($indir,$outdir,$paired,$flowcell,$logfile,$bamfiles_bwa_freebayes,$bamfiles_bowtie2_freebayes,$refsample);
my (@fastqfilelist,@bamfilelist,@ar_reports,@vf_reports);
my $date = GetDate();
my ($refdir,$scriptdir) = GetDirs("$config_dir/$config_file");
my $i = my $x = 0;
my $procs = 8;        #number of processors to shoot for
my $ref = 1;
my $project;
my $library = 0;
my $assemble = 0;
my $st = 0;
my $gzipped = 0;       #sample raw reads are gzipped

###process command line script options###
my %cmds = ReturnCmds(@ARGV);
if ($cmds{'h'}){
  $help .= GetRefList("$config_dir/$config_file");
  $help .= "\n************\nAuthor:\t\tCris Vera\nGroup:\t\tIGOH Theme, IGB, University of Illinois\nContact:\tjcvera\@illinois.edu\n";  
  die "\n$help\n$usage\n";
}
$indir = $cmds{'d'} if ($cmds{'d'});
$outdir = $cmds{'e'} if ($cmds{'e'});
$logfile = $cmds{'o'} if ($cmds{'o'});
$procs = $cmds{'t'} if ($cmds{'t'});
$project = $cmds{'p'} if ($cmds{'p'});
$assemble = $cmds{'u'} if ($cmds{'u'});
$st = $cmds{'m'} if ($cmds{'m'});
$refsample = $cmds{'r'} if ($cmds{'r'});
$ref = $cmds{'R'} if ($cmds{'R'});
$gzipped = $cmds{'z'} if ($cmds{'z'});
$library = $cmds{'k'} if ($cmds{'k'});
die "\nError: must specify input directory: $indir\n\n" if (!$indir);
die "\nError: must specify output directory.\n\n" if (!$outdir);
die "\nError: must specify project base.\n\n" if (!$project);

###tool directories###
my %tools = GetConfigs('TOOL',"$config_dir/$config_file");
my $fastqc = $tools{'FastQC'};
my $trimmomatic = $tools{'Trimmomatic'};
my $addreadgroup = "$tools{'Picard'} AddOrReplaceReadGroups";
my $fastqscreen = $tools{'FastqScreen'};
my $sambamba = $tools{'Sambamba'};
my $centrifuge = $tools{'Centrifuge'};
my $qualimap = $tools{'Qualimap'};
my $bwa = "$tools{'BWA'} mem";
my $bamtools = $tools{'Bamtools'};
my $unicycler = $tools{'Unicycler'};
my $ariba = "$tools{'Ariba'} run";
my $summary = "$tools{'Ariba'} summary";
my $spades = $tools{'Spades'};
my $multiqc = $tools{'MultiQC'};

###script directories###
my $getRG = "perl $scriptdir/GetRG_fromFASTQ.pl";
my $gatherSTs = "perl $scriptdir/GatherSTinfo.pl";

###miscellaneous###
my %misc = GetConfigs('MISC',"$config_dir/$config_file");
my $adapters = $misc{'Adapters'};
my $centrifuge_ref = $misc{'CentrifugeDB'};
my $ariba_ar_ref = $misc{'AribaAR'};
my $ariba_vf_ref = $misc{'AribaVF'};
my $trimparams = $misc{'TrimParameters'};

###reference genomes###
my %refs = GetConfigs('REF',"$config_dir/$config_file");
my $bwa_ref = $refs{$ref}{'BWA'};
my $screenconfile = $refs{$ref}{'FastqScreenConf'};
my $ariba_st_ref = $refs{$ref}{'AribaST'};
my $refname = $refs{$ref}{'Refname'};

##library kits##
my @libraries = (
  'Unknown',
  'Nextera_XT_DNA_Library_Prep_Kit',
  'TruSeq_RNA_Sample_Prep_Kit',
  'TruSeq_Stranded_mRNA_LT_Kit',
  'TruSeq_Stranded_Total_RNA_LT_Kit',
  'TruSeq_Small_RNA_Library_Prep_Kit',
  'TruSeq_Chip_Library_Prep_Kit',
  'TruSeq_DNA_Library_Prep_Kit',
  'NEBNext_Multiplex_Small_RNA_Library_Prep_Kit',
  'Clonetech_SMARTer_Ultra_Low_Input_RNA_Kit',
  'SureSelect_XT_Mouse_All_Exome',
  'SureSelect_XT_Human_All_Exon_V5-UTR',
  'iGenomx_Riptide_DNA_HT-RLP',
  'Kappa_DNA_HTP_Prep_Kit',
  'Kappa_PCR-free_DNA_HTP_Prep_Kit',
  'Nextera_Flex_DNA_Library_Prep_Kit',
);
##get library if integer
if ($library =~ m/^([0-9]+)$/){
  $library = $libraries[$1];
}

##make absolute paths
$indir = File::Spec->rel2abs($indir);
$outdir = File::Spec->rel2abs($outdir);

#create primary project directory
my $logout = "\nCreating project directory:\n";
if (!-e $outdir){
  mkdir $outdir;
  mkdir "$outdir/bamlinks";
  $logout .= "\tOut Directory created: $outdir\n";
}
else{
  $logout .= "\tSkipping out directory creation, it already exists: $outdir\n";
}

##logfile start
$logfile = "$project\_Pipeline.log";
if (!-e "$outdir/$logfile"){
  open ('LOG', ">$outdir/$logfile") or die "Cannot create $outdir/$logfile: $!\n";
  print LOG "Run date: $date\n";
  print LOG "Threads: $procs\n";
  print LOG "\nProject name: $project\n";
  print LOG "Input Dir: $indir\n";
  print LOG "Output Dir: $outdir\n";
  print LOG "\nLibrary: $library\n";
  print LOG "\nReference genome: $refname\n";
  print LOG "Assemble Reads: Unicycler\n" if ($assemble == 1 or $assemble == 2);
  print LOG "Assemble Reads: Spades\n" if ($assemble == 3 or $assemble == 4);
  print LOG "Assemble Reads: No\n" if (!$assemble);
  print LOG "Run Ariba ST: Yes\n" if ($st == 1);
  print LOG "Run Ariba ST: No\n" if (!$st);
  print LOG "Run Ariba ST, AR, & VF: Yes\n" if ($st == 2);
}
else{
  open ('LOG', ">>$outdir/$logfile") or die "Cannot write to $outdir/$logfile: $!\n";
  print LOG "\n\n**************\nRe-Running Pipeline\n**************\n";
  print LOG "Run date: $date\n";
  print LOG "Threads: $procs\n";
  print LOG "\nProject name: $project\n";
  print LOG "Input Dir: $indir\n";
  print LOG "Output Dir: $outdir\n";
  print LOG "\nLibrary: $library\n";
  print LOG "\nReference genome: $refname\n";
  print LOG "Assemble Reads: Unicycler\n" if ($assemble == 1 or $assemble == 2);
  print LOG "Assemble Reads: Spades\n" if ($assemble == 3 or $assemble == 4);
  print LOG "Assemble Reads: No\n" if (!$assemble);
  print LOG "Run Ariba ST: Yes\n" if ($st == 1);
  print LOG "Run Ariba ST: No\n" if (!$st);
  print LOG "Run Ariba ST, AR, & VF: Yes\n" if ($st == 2);
}
print LOG "$logout\n";

#get all sample file names to run script on
opendir(INDIR,$indir) or die "Can't open directory: $indir: $!\n";
my @samplefiles;
@samplefiles = grep {m/_R1_001\.fastq$/ && -f "$indir/$_"} readdir(INDIR) if (!$gzipped);
if ($gzipped){
  my @temp = grep {m/_R[12]_001\.fastq.gz$/ && -f "$indir/$_"} readdir(INDIR);
  close(INDIR);
  my $g = my $gg = 0;
  foreach my $gz (@temp){
    my $ungz = $gz;
    $ungz =~ s/\.gz$//;
    if (!-e "$indir/$ungz"){
      my $ungzip = `gunzip -cd $indir/$gz >$indir/$ungz`;
      $gg += 1;
    }
    $g += 1;
  }
  opendir(INDIR,$indir) or die "Can't open directory: $indir: $!\n";
  @samplefiles = grep {m/_R1_001\.fastq$/ && -f "$indir/$_"} readdir(INDIR);
  print STDERR "\nGzipped FASTQ files found: $g\nFASTQ files ungzipped: $gg\n\n";
  print LOG "\nGzipped FASTQ files found: $g\nFASTQ files ungzipped: $gg\n\n";
}
close(INDIR);

###create directories, and run pipeline commands
foreach my $samplefile1 (@samplefiles){
  $i += 1;
  my $run = 1;
  my $samplefile2 = my $samplename = $samplefile1;
  $samplename =~ s/_R1_001\.fastq//;
  $samplefile2 =~ s/_R1_001\.fastq/_R2_001.fastq/;
  my $sample = $samplename;
  $sample =~ s/_[ATGC-]+_(L00[0-9])$/_$1/;
  $sample =~ s/-pseudomonas-[ATGC-]+$//;  ##some raw read files had '-pseudomonas' added to them for some reason
  #create dirs
  if (!-e "$outdir/Sample_$samplename"){
    mkdir "$outdir/Sample_$samplename";
    mkdir "$outdir/Sample_$samplename/qualimap";
    mkdir "$outdir/Sample_$samplename/fastqc";
    mkdir "$outdir/Sample_$samplename/centrifuge";
    mkdir "$outdir/Sample_$samplename/fastqscreen";
    mkdir "$outdir/Sample_$samplename/ariba" if ($st);
    mkdir "$outdir/Sample_$samplename/spades" if ($assemble == 3 or $assemble == 4);
    print STDERR "\nCreate sample directory and sub-directories for $samplename: $outdir/Sample_$samplename\n";
    print LOG "\nCreate sample directory and sub-directories for $samplename: $outdir/Sample_$samplename\n";
  }
  else{
    print STDERR "\nSkip directory creation for $samplename: it already exists.\n\n";
    print LOG "\nSkip directory creation for $samplename: it already exists.\n\n";
    $run = 0;
  }
  
  #command variables
  $paired = 1 if (-e "$indir/$samplefile2");
  my ($assemble_cmd,$breseq_cmd,$ariba_st_cmd,$ariba_ar_cmd,$ariba_vf_cmd,$trim_cmd,$fastqc_cmd,$fastqscreen_cmd,$sort_cmd1,$sort_cmd2,$markedup_cmd1,$markedup_cmd2,$bwa_cmd,$addgroups_cmd1,$addgroups_cmd2,$qualimap_cmd,$rg_cmd,$bamtools_cmd,$index_cmd,$filter_cmd,$centrifuge_cmd);
  
  ##commands for paired-end
  if ($paired){
    ##java -Xmx64g -jar ##in case we go back to 64g
    $rg_cmd = "$getRG -i $indir/$samplefile1 -o $outdir/Sample_$samplename/RG.txt";
    $trim_cmd = "$trimmomatic -threads $procs -phred33 -trimlog $outdir/Sample_$samplename/$sample\_trim.log $indir/$samplefile1 $indir/$samplefile2 $outdir/Sample_$samplename/$samplename\_R1_trimmed.fastq $outdir/Sample_$samplename/$samplename\_R1_trimmed_unpaired.fastq $outdir/Sample_$samplename/$samplename\_R2_trimmed.fastq $outdir/Sample_$samplename/$samplename\_R2_trimmed_unpaired.fastq ILLUMINACLIP:$adapters:$trimparams 2>$outdir/Sample_$samplename/$sample\_run_trimmomatic.log && rm $outdir/Sample_$samplename/$sample\_trim.log";
    $fastqc_cmd = "$fastqc -o $outdir/Sample_$samplename/fastqc/ --noextract -k 5 -t $procs -f fastq $indir/$samplefile1 $indir/$samplefile2 $outdir/Sample_$samplename/$samplename\_R1_trimmed.fastq $outdir/Sample_$samplename/$samplename\_R2_trimmed.fastq 2>$outdir/Sample_$samplename/fastqc/$sample\_run_fastqc.err";
    $fastqscreen_cmd = "$fastqscreen --outdir $outdir/Sample_$samplename/fastqscreen/ --subset 1000000 --threads $procs --nohits --conf $screenconfile --aligner bwa  $outdir/Sample_$samplename/$samplename\_R1_trimmed.fastq 2>$outdir/Sample_$samplename/fastqscreen/$sample\_fastqscreen.log ";
    $centrifuge_cmd = "$centrifuge -p $procs -x $centrifuge_ref --met-file $outdir/Sample_$samplename/centrifuge/$sample\_metrics.txt -1 $outdir/Sample_$samplename/$samplename\_R1_trimmed.fastq -2 $outdir/Sample_$samplename/$samplename\_R2_trimmed.fastq -S $outdir/Sample_$samplename/centrifuge/$samplename\_trimmed.centrifuge.out --report-file $outdir/Sample_$samplename/centrifuge/$sample\_centrifuge_report.tsv 2>$outdir/Sample_$samplename/centrifuge/centrifuge.err && $centrifuge-kreport -x $centrifuge_ref $outdir/Sample_$samplename/centrifuge/$samplename\_trimmed.centrifuge.out >$outdir/Sample_$samplename/centrifuge/$samplename\_trimmed.centrifuge.kraken";
    $ariba_st_cmd = "$ariba --threads $procs $ariba_st_ref $indir/$samplefile1 $indir/$samplefile2 $outdir/Sample_$samplename/ariba/mlst/ 2>$outdir/Sample_$samplename/ariba/ariba_st.err";
    $ariba_ar_cmd = "$ariba --threads $procs $ariba_ar_ref $indir/$samplefile1 $indir/$samplefile2 $outdir/Sample_$samplename/ariba/megares/ 2>$outdir/Sample_$samplename/ariba/ariba_ar.err";
    $ariba_vf_cmd = "$ariba --threads $procs $ariba_vf_ref $indir/$samplefile1 $indir/$samplefile2 $outdir/Sample_$samplename/ariba/vfdb/ 2>$outdir/Sample_$samplename/ariba/ariba_vf.err";
    $assemble_cmd = "$unicycler -t $procs -1 $indir/$samplefile1 -2 $indir/$samplefile2 -o $outdir/Sample_$samplename/unicycler/ --mode conservative --spades_path $spades" if ($assemble == 1);
    $assemble_cmd = "$unicycler -t $procs -1 $outdir/Sample_$samplename/$samplename\_R1_trimmed.fastq -2 $outdir/Sample_$samplename/$samplename\_R2_trimmed.fastq -o $outdir/Sample_$samplename/unicycler/ --mode conservative --spades_path $spades --no_correct" if ($assemble == 2);
    $assemble_cmd = "$spades --careful -t $procs -1 $indir/$samplefile1 -2 $indir/$samplefile2 -o $outdir/Sample_$samplename/spades/ 2>$outdir/Sample_$samplename/spades/spades.err" if ($assemble == 3);
    $assemble_cmd = "$spades --careful --only-assembler -t $procs -1 $outdir/Sample_$samplename/$samplename\_R1_trimmed.fastq -2 $outdir/Sample_$samplename/$samplename\_R2_trimmed.fastq -o $outdir/Sample_$samplename/spades/ 2>$outdir/Sample_$samplename/spades/spades.err" if ($assemble == 4);
    $bwa_cmd = "$bwa -t $procs $bwa_ref $outdir/Sample_$samplename/$samplename\_R1_trimmed.fastq $outdir/Sample_$samplename/$samplename\_R2_trimmed.fastq 2>$outdir/Sample_$samplename/$sample\_bwa-mem.err | $sambamba view -S -f bam /dev/stdin >$outdir/Sample_$samplename/$sample\_BWA-mem.bam";
    $addgroups_cmd1 = "$addreadgroup INPUT=$outdir/Sample_$samplename/$sample\_BWA-mem.bam OUTPUT=$outdir/Sample_$samplename/$sample\_BWA-mem.RG.bam VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=LENIENT COMPRESSION_LEVEL=5";
    $sort_cmd1 = "$sambamba sort -t $procs -m 32G --tmpdir $outdir/Sample_$samplename/temp -o $outdir/Sample_$samplename/$sample\_BWA-mem.RG.sorted.bam $outdir/Sample_$samplename/$sample\_BWA-mem.RG.bam 2>$outdir/Sample_$samplename/$sample\_sort_bwa.err && rm $outdir/Sample_$samplename/$sample\_BWA-mem.bam $outdir/Sample_$samplename/$sample\_BWA-mem.RG.bam";
    $markedup_cmd1 = "$sambamba markdup --sort-buffer-size=8192 --overflow-list-size=1000000 -t $procs --tmpdir=$outdir/Sample_$samplename/temp $outdir/Sample_$samplename/$sample\_BWA-mem.RG.sorted.bam $outdir/Sample_$samplename/$sample\_BWA-mem.RG.sorted.markedup.bam 2>$outdir/Sample_$samplename/$sample\_markdup_bwa.err && rm $outdir/Sample_$samplename/$sample\_BWA-mem.RG.sorted.bam $outdir/Sample_$samplename/$sample\_BWA-mem.RG.sorted.bam.bai";
    $filter_cmd = "$sambamba view -F \"not secondary_alignment and not supplementary\" -f bam $outdir/Sample_$samplename/$sample\_BWA-mem.RG.sorted.markedup.bam >$outdir/Sample_$samplename/$sample\_BWA-mem.RG.sorted.markedup.filtered.bam";
    $bamtools_cmd = "$bamtools stats -in $outdir/Sample_$samplename/$sample\_BWA-mem.RG.sorted.markedup.filtered.bam 1>$outdir/Sample_$samplename/$sample\_BWA-mem.RG.sorted.markedup.filtered.bam_stats 2>$outdir/Sample_$samplename/$sample\_BWA-mem.RG.sorted.markedup.filtered.bam_stats.err && $bamtools stats -in $outdir/Sample_$samplename/$sample\_BWA-mem.RG.sorted.markedup.bam 1>$outdir/Sample_$samplename/$sample\_BWA-mem.RG.sorted.markedup.bam_stats 2>$outdir/Sample_$samplename/$sample\_BWA-mem.RG.sorted.markedup.bam_stats.err";
    $qualimap_cmd = "$qualimap bamqc -bam $outdir/Sample_$samplename/$sample\_BWA-mem.RG.sorted.markedup.filtered.bam -ip -c -outdir $outdir/Sample_$samplename/qualimap -outformat PDF -nt $procs --java-mem-size=64G -nw 500 -p NON-STRAND-SPECIFIC 2>$outdir/Sample_$samplename/qualimap/$sample\_qualimap.err";
  }
  ##commands for single-end
  else{
    die "\nno single-end reads supported!\n\n";
  }
  ##run commands
  if ($run){
    print STDERR "\n\nRunning jobs for sample $i: $sample\n";
    print LOG "\n\nRunning jobs for sample $i: $sample\n";
    print STDERR "\n\nGetting flowcell and lane information from FASTQ file read headers for sample $sample...";
    print LOG "\n\nGetting flowcell and lane information from FASTQ file read headers for sample $sample...";
    print LOG "\nGet RG Command: $rg_cmd -k $library -s $sample 2>$outdir/Sample_$samplename/RG.err\n\n";
    my $run_rg = `$rg_cmd -k $library -s $sample 2>$outdir/Sample_$samplename/RG.err`;
    print STDERR " Done!\n\nRG: $run_rg\n\nProceeding with pipeline...";
    print LOG "Trim Command: $trim_cmd\n\n";
    my $run_trim = `$trim_cmd`;
    ##run optional assembly commands
    if ($assemble and (($assemble == 1 or $assemble == 2) and !-e "$outdir/Sample_$samplename/unicycler/assembly.fasta") or (($assemble == 3 or $assemble == 4) and !-e "$outdir/Sample_$samplename/spades/scaffolds.fasta")){
      print STDERR "\nRunning optional Unicycler assemblies:\nCommand: $assemble_cmd\n" if ($assemble == 1 or $assemble == 2);
      print LOG "\nRunning optional Unicycler assemblies:\nCommand: $assemble_cmd\n" if ($assemble == 1 or $assemble == 2);
      print STDERR "\nRunning optional Spades assemblies:\nCommand: $assemble_cmd\n" if ($assemble == 3 or $assemble == 4);
      print LOG "\nRunning optional Spades assemblies:\nCommand: $assemble_cmd\n" if ($assemble == 3 or $assemble == 4);
      my $run_assemble = `$assemble_cmd`;
    }
    elsif ($assemble){
      print STDERR "\nAssembly FASTA file already exists! Skipping optional assembly for sample: $sample\n";
      print LOG "\nAssembly FASTA file already exists! Skipping optional assembly for sample: $sample\n";
      #$run = 0 if (-e "$outdir/Sample_$samplename/RG.txt");
    }
    ##run optional Ariba commands
    if ($st and !-e "$outdir/Sample_$samplename/ariba/mlst/"){
      print STDERR "\nRunning optional Ariba ST (MLST):\nCommand: $ariba_st_cmd\n";
      print LOG "\nRunning optional Ariba ST (MLST):\nCommand: $ariba_st_cmd\n";
      my $run_ariba_st = `$ariba_st_cmd`;
      if ($st == 2){
        print STDERR "\nRunning optional Ariba AR (Megares):\nCommand: $ariba_ar_cmd\n";
        print LOG "\nRunning optional Ariba AR (Megares):\nCommand: $ariba_ar_cmd\n";
        my $run_ariba_ar = `$ariba_ar_cmd`;
        print STDERR "\nRunning optional Ariba VF (VFDB):\nCommand: $ariba_vf_cmd\n";
        print LOG "\nRunning optional Ariba VF (VFDB):\nCommand: $ariba_vf_cmd\n";
        my $run_ariba_vf = `$ariba_vf_cmd`;
      }
    }
    elsif ($st){
      print STDERR "\nAriba directories already exists! Skipping optional Ariba for sample: $sample\n";
      print LOG "\nAriba directories already exists! Skipping optional Ariba for sample: $sample\n";
      $run = 0 if (-e "$outdir/Sample_$samplename/RG.txt");
    }
    ##$readgroup_info = "RGID=$flowcell\_$sample RGSM=$samplename RGLB=$library RGPL=ILLUMINA RGPU=$flowcell.$lane RGCN=NCI-CCRSF";
    print LOG "Fastqc Command: $fastqc_cmd\n\n";
    my $run_fastqc = `$fastqc_cmd`;
    print LOG "Fastqscreen Command: $fastqscreen_cmd\n\n";
    my $run_fastqscreen = `$fastqscreen_cmd`;
    print LOG "Centrifuge Command: $centrifuge_cmd\n\n";
    my $run_centrifuge = `$centrifuge_cmd`;
    print LOG "BWA-mem Command: $bwa_cmd\n\n";
    my $run_bwa = `$bwa_cmd`;
    
    #run addgroups, from FASTQ file scan info
    print LOG "Picard (AddReadGroups) BWA Command: $addgroups_cmd1 $run_rg\n\n";
    my $run_addgroups = `$addgroups_cmd1 $run_rg 2>$outdir/Sample_$samplename/picard_addgroups_bwa.err`;
    print LOG "Sambamba (sort) BWA Command: $sort_cmd1\n\n";
    my $run_sort = `$sort_cmd1`;
    print LOG "Sambamba (markedup) BWA Command: $markedup_cmd1\n\n";
    my $run_markedup = `$markedup_cmd1`;
    print LOG "Sambamba View/Filter BWA Command: $filter_cmd\n\n";
    my $run_filter = `$filter_cmd`;
    print LOG "Bamtools stats BWA Command: $bamtools_cmd\n\n";
    my $run_bamtools = `$bamtools_cmd`;
    print LOG "Qualimap BWA Command: $qualimap_cmd\n\n";
    my $run_qualimap = `$qualimap_cmd`;
    print STDERR " Done!\n\n";
    $x += 1;
  }
  else{
    print STDERR "\n\nSkipping jobs for sample $i: $sample\n";
    print LOG "\n\nSkipping jobs for sample $i: $sample\n";
  }
  
  #report lists for Ariba
  push @ar_reports,"$outdir/Sample_$samplename/ariba/megares/report.tsv\t$sample";
  push @vf_reports,"$outdir/Sample_$samplename/ariba/vfdb/report.tsv\t$sample";
}
##print some stats to STDERR/LOG
print STDERR "\nTotal input samples found: $i\n";
print STDERR "Samples processed: $x\n";
print LOG "\nTotal samples found: $i\n";
print LOG "Samples processed: $x\n";

##write report file lists for Ariba summary
my $ar_file = "$outdir/$project\_Ariba_AR_FileList.txt";
my $vf_file = "$outdir/$project\_Ariba_VF_FileList.txt";
if ($st){
  open (AR, ">$ar_file") or die "Cannot create $ar_file: $!\n";
  foreach my $file (@ar_reports){
    print AR "$file\n";
  }
  close(AR);
  open (VF, ">$vf_file") or die "Cannot create $vf_file: $!\n";
  foreach my $file (@vf_reports){
    print VF "$file\n";
  }
  close(VF);
}
print STDERR "\n\nAll individual sample level jobs have finished! Running project level jobs...\n\n";
print LOG "\n\nAll individual sample level jobs have finished! Running project level jobs...\n\n";

##run Ariba Summary stats
if ($st){
  my $st_summary_cmd = "$gatherSTs -d $outdir -o $outdir/$project\_MLST_SummaryReport.txt 2>$outdir/$project\_MLST_SummaryReport.err";
  print STDERR "Running Ariba ST summary (GatherSTinfo.pl): $st_summary_cmd\n\n";
  print LOG "Ariba ST summary (GatherSTinfo.pl): $st_summary_cmd\n\n";
  my $summary_run = `$st_summary_cmd`;
  my $ar_summary_cmd = "$summary --col_filter n --known_variants --novel_variants -f $ar_file $outdir/$project\_AR_SummaryReport 2>$outdir/$project\_AR_SummaryReport.err";
  my $vf_summary_cmd = "$summary --col_filter n --known_variants --novel_variants -f $vf_file $outdir/$project\_VF_SummaryReport 2>$outdir/$project\_VF_SummaryReport.err";
  print STDERR "Running Ariba AR Summary: $ar_summary_cmd\n\n";
  print LOG "Ariba AR Summary command: $ar_summary_cmd\n\n";
  $summary_run = `$ar_summary_cmd`;
  print STDERR "Running Ariba VF Summary: $vf_summary_cmd\n\n";
  print LOG "Ariba VF Summary command: $vf_summary_cmd\n\n";
  $summary_run = `$vf_summary_cmd`;
}
print STDERR "\nFinished project level jobs!\nPipeline finished!\n\n";
print LOG "\nFinished project level jobs!\nPipeline finished!\n\n";

##run MultiQC
my $multiqc_cmd = "$multiqc -n $outdir/$project\_MultiQC_SummaryReport -o $outdir $outdir 2>$outdir/MultiQC.err";
print STDERR "Running MultiQC: $multiqc_cmd\n\n";
print LOG "MultiQC command: $multiqc_cmd\n\n";
my $run_multiqc = `$multiqc_cmd`;

##get all config values sorted into integers for selection
sub GetConfigs{
  my ($reftype,$configfile) = @_;
  my %refs;
  my %refs2;
  open (CON, "<$configfile") or die "Cannot open $configfile: $!\n";
  while (my $line = <CON>){
    chomp $line;
    if ($line =~ m/^REF\_([^\t]+)\t([^\t]+)\t(.+)/){
      $refs{$1}{$2} = $3;
    }
    elsif ($line =~ m/^HREF\_([^\t]+)\t([^\t]+)\t(.+)/){
      $refs{$1}{$2} = $3;
    }
    elsif ($line =~ m/^$reftype\_([^\t]+)\t(.+)/){
      $refs{$1} = $2;
    }
  }
  close(CON);
  my $n = 0;
  if ($reftype eq 'REF'){
    foreach my $ref (sort keys %refs){
      $n++;
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
    if ($line =~ m/^REF_([^\t]+)\t([^\t]+)\t([^\t]+)/){
      $refs{$1} += 1;
    }
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

