#!/usr/bin/perl

# ************
# Author: J. Cristobal Vera
# email: jvera8888@gmail.com

#PrepGenome: prepare a virus or host genome for use in the VirusPipeline for counting virus in samples. Can 'hollow out' the genome
#by removing interior portion of the virus sequence (to prevent 'noise' read counts and/or avoid most CRISPR sites).  Can also
#remove target virus sequence from Host Reference genome, since this pipeline expects a virus free reference in order to work correctly.
#Can also insert a virus genome into the host genome at the specified position (see notes below).
##notes:
##-to insert a virus sequence, first reverse compliment the virus genome if needed (i.e. if it is in RC orientation in the host genome)
##-add an 'N' followed by the 5bp duplicated sequence to the end of the virus genome if its a mu phage

# defaults, initializations, and constants
my $help = "\n\nPrepGenome.\n".
            "\t-i  Option: Input file (tab-delimited). Default=STDIN.\n".
            "\t-o  Option: Output file (tab-delimited). Default=STDOUT.\n".
            "\t-b  Option: Border lengths on Virus Genome. Optional, Default=1000\n".
            "\t-r  Option: Virus region coordinates to remove from or add to Host Genome (Chr:Start:Stop). Optional.\n\t\t    E.G. -r 1:2500:7000\n".
            "\t-v  Option: Virus genome to add to host genome. Optional.\n\t    Use option -r for insert coordinates.\n".
            "\t-f  Option: Fill emptied virus region in Host/Virus genome with Ns. Default=0.\n".
            "\n************\nAuthor: J. Cristobal Vera, email: jvera8888\@gmail.com\n"; 
my $usage = "\nPrepGenome.pl -i [Input File] -o [Output File] -b [Border Lengths] -r [Remove/Add Coordinates] -v [Virus Genome] -f [Fill]\n";
my $infh = 'STDIN';
my $outfh = 'STDOUT';
my $i = 0;
my $border = 1000;
my $virusregion = 0;
my $virus = 0;
my $Ns = 1;
my $fill = 0;
my @virus;

#process command line custom script tags
my %cmds = ReturnCmds(@ARGV);
die "\n$help\n$usage\n" if ($cmds{'h'});
if ($cmds{'i'}) {
  $infh = 'IN';
  open ($infh, "<$cmds{'i'}") or die "Cannot open $cmds{'i'}: $!\n";
}
if ($cmds{'o'}) {
  $outfh = 'OUT';
  open ($outfh, ">$cmds{'o'}") or die "Cannot create $cmds{'o'}: $!\n";
}
if ($cmds{'v'}) {
  open (VIR, "<$cmds{'v'}") or die "Cannot open $cmds{'v'}: $!\n";
  chomp(my @viruses = <VIR>);
  close(VIR);
  @virus = split /\t/,$viruses[0];
  $virus = 1;
}
$border = $cmds{'b'} if ($cmds{'b'});
$virusregion = $cmds{'r'} if ($cmds{'r'});
$fill = $cmds{'f'} if ($cmds{'f'});

# Main
while (my $line = <$infh>) {
  chomp $line;
  next if ($line eq '');
  $i += 1;
  my @line = split /\t/,$line;
  my @seq = split //,$line[1];
  my $header = $line[0];
  $header =~ s/^([^ ]+) .+/$1/;
  my $seqlen = length $line[1];
  my @borders = split /:/,$virusregion;
  if (!$virusregion){
    die "\nError: border regions are too long (i.e. longer than sequence length):\n\tsequence length: $seqlen\n\ttotal border length: ".($border*2)."\n\n" if ($border*2 > $seqlen);
    $Ns = 100;  ##edit to add fill functionality here
    splice @seq,($border),-($border),('N' x $Ns); #removes internal portion of sequence; adds Ns inbetween
  }
  elsif ($virusregion and !$virus and $header eq $borders[0]){
    my $chunklen = ($borders[2]-$borders[1])+1;
    $Ns = $chunklen if ($fill);
    splice @seq,($borders[1]-1),$chunklen,('N' x $Ns);
  }
  elsif ($virusregion and $virus and $header eq $borders[0]){
    splice @seq,$borders[1],0,'N'.$virus[1];
  }
  $line[1] = join "",@seq;
  print $outfh "$line[0]\t$line[1]\n";
}
print STDERR "Total lines parsed: $i\n\n";

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
