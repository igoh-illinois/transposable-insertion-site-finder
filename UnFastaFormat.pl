#!/usr/bin/perl

# UnfastaFormat 1.0
# Accepts FASTA formatted input and returns tab-delimited output.  Output fields are 
# title and sequence.  Optionally, removes single quotes (') characters from title strings.
# Adds extra spaces at new lines for quality scores if option -q is on.
# Optionally, UnFastaFormat replaces spaces in titles with another character.
# ************
# Author: J. Cristobal Vera
# email: jvera8888@gmail.com


# defaults, initializations, and constants
my $help = "\nUnFastaFormat v1.0\nAccepts FASTA formatted input and returns tab-delimited output.  Outputs two columns with fields for title and sequence.\n".
          "\nOptions-Switches:\n".
          "\t-i  Option: allows input file to be specified.  Optional,\n\t    default is STDIN.\n".
          "\t-o  Option: allows output file to be specified.  Optional,\n\t    default is STDOUT.\n".
          "\t-s  Switch: causes single quote (') characters to be\n\t    removed from titles.\n".
          "\t-r  Option: allows spaces in titles to be replaced with\n\t    another character.  Optional.\n".
          "\t-q  Switch: convert quality scores to fasta format.\n".
          "\t-p  Switch: platform safe mode.  Removes ALL newlines and\n\t    carriage returns as well as tabs and spaces from end of strings.\n".
          "\n************\nAuthor: J. Cristobal Vera, email: cris.vera\@nih.gov\n";
my $usage = "\nUsage:\nUnFastaFormat -i [Input File] -o [Output File] -r [Replacement Character] -q {Quality Scores} -s {Remove Single Quotes} -p {Platform Safe}\n";
my $infh = 'STDIN';
my $outfh = 'STDOUT';
my $rep = 0;
my $qual = 0;
my $p = 0;
my $seq;
my $s;
my $i = 0;

#process command line custom script options
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
$rep = $cmds{'r'} if ($cmds{'r'});
$qual = 1 if ($cmds{'q'});
$s = 1 if ($cmds{'s'});
$p = 1 if ($cmds{'p'});


while ($line = <$infh>) {
  chomp $line;
  $line = SuperChomp($line) if ($p);
  if ($line =~ m/^>/){
    $i += 1;
    $line =~ s/^>//;
    $line =~ s/'//g if ($s);
    $line =~ s/ /$rep/g if ($rep);
    if ($i > 1){
      $seq =~ s/^ +//;
      $seq =~ s/ +$//;
      $seq =~ s/ {2,}/ /g;
      print $outfh "$seq\n";
      $seq = '';
    }
    print $outfh "$line\t"    
    }
  else {
    if ($qual){
      $seq .= $line.' ';
    }
    else{
      $seq .= $line;
    }
  }
}
$seq =~ s/ {2,}/ /g;
$seq =~ s/^ +//;
$seq =~ s/ +$//;
print $outfh "$seq\n";
print STDERR "\nSequences converted to tab-delimited: $i\n";

# Removes ALL newlines and carriage returns from a string.  Also removes tabs and spaces from end of string.
sub SuperChomp{
  my ($string) = @_;
  my $gone = 0;
  $string =~ s/\n//g;
  $string =~ s/\r//g;
  until ($gone){
    $gone = 1;
    if ($string =~ s/ +$//){
      $gone = 0;
    }
    if ($string =~ s/\t+$//){
      $gone = 0;
    }
  }
  return $string;
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
