#!/usr/bin/perl

# ************
# Author: J. Cristobal Vera
# email: jvera8888@gmail.com

##Description: combine two SAM files and keep only unique alignment rows. Keeps headers from first SAM (-i).

# defaults, initializations, and constants
my $help = "\n\nCombineSAMs.\nDescription: combine two SAM files and keep only unique alignment rows. Keeps headers from first SAM (-i).\n".
            "\t-i  Option: Input file. Required.\n".
            "\t-j  Option: Input file. Required.\n".
            "\t-o  Option: Output file. Default=STDOUT.\n".
            "\n************\nAuthor: J. Cristobal Vera, email: jvera8888\@gmail.com\n"; 
my $usage = "\nCombineSAMs.pl -i [Input File] -o [Output File]\n";
my $outfh = 'STDOUT';
my $i = my $n = 0;
my %sam;

#process command line custom script tags
my %cmds = ReturnCmds(@ARGV);
die "\n$help\n$usage\n" if ($cmds{'h'});
if ($cmds{'i'}) {
  open (IN1, "<$cmds{'i'}") or die "Cannot open $cmds{'i'}: $!\n";
  chomp(my @sam = <IN1>);
  close(IN1);
  GetSam(1,\%sam,@sam);
}
if ($cmds{'j'}) {
  open (IN2, "<$cmds{'j'}") or die "Cannot open $cmds{'j'}: $!\n";
  chomp(my @sam = <IN2>);
  close(IN2);
  GetSam(0,\%sam,@sam);
}
if ($cmds{'o'}) {
  $outfh = 'OUT';
  open ($outfh, ">$cmds{'o'}") or die "Cannot open $cmds{'o'}: $!\n";
}

# print out
my @headers = split /\=\+\=/,$sam{'header'};
foreach my $header (@headers){
  print $outfh "$header\n";
}
delete $sam{'header'};
foreach my $line (sort keys %sam){
  print $outfh "$line\n";
  $i += 1;
  $n += $sam{$line} - 1;
}
print STDERR "\nTotal unique alignments kept: $i\n";
print STDERR "Total redundant alignments dropped: $n\n\n";

sub GetSam{
  my ($x,$sam,@sam) = @_;
  foreach my $line (@sam){
    if ($line =~ m/^\@/){
      $sam->{'header'} .= "=+=$line" if ($x and exists $sam->{'header'});
      $sam->{'header'} = $line if ($x and !exists $sam->{'header'});
    }
    else{
      $sam->{$line} += 1;
    }
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
