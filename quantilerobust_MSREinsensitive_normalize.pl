#takes 3 arguments
#1 - name of file with the data
#2 - name of file with the UNCUTSNP sets (probe_ids)
#3 - name of file with the CUTSNP sets (probe_ids)
#4 - name of file with the model data (order of probe_ids doesn't matter)

use strict;
use warnings;

use Getopt::Long;

my $help;
my $datafile='';
my $mnrfile='';
my $mprfile='';
my $modelfile='';

GetOptions(
  'datafile|i=s' => \$datafile,
  'mnrfile|u=s' => \$mnrfile,
  'mprfile|c=s' => \$mprfile,
  'modelfile|m=s' => \$modelfile,
  'help|h' => \&usage);


# usage() if ( @ARGV < 1 or $help eq 1);

if (defined $datafile && $datafile eq '') {
print("raw intensity datafile not defined\nUse --datafile/-i\n");
die;
}
if (defined $mnrfile && $mnrfile eq '') {
print("MNR probe id file not defined\nUse --mnrfile/-u\n");
die;
}
if (defined $mprfile && $mprfile eq '') {
print("MPR probe id file not defined\nUse --mprfile/-c/\n");
die;
}
if (defined $modelfile && $modelfile eq '') {
print("model sample raw intensity file not defined\nUse --modelfile/-m\n");
die;
}

##filenames
my $rawuncutdata = "rawUNCUT.txt";
my $rawcutdata = "rawCUT.txt";
my $quantuncutdata = "quantUNCUT.txt";
my $adjustdata = "quant-rawUNCUTadjustments.txt";
my $temp = "tempfile.txt";


#open model file and combine with datafile
#open model file, read into hash with probeid as key
open (MODEL, $modelfile);
my ($probeid, $data, $modelcolheader, $line);
my $modelheader =<MODEL>;
my %model;
while (<MODEL>) {
chomp;
 ($probeid, $data) = split /\t/;
  $model{$probeid}=$data;
}
close MODEL;

#open data file, read into hash with probeid as key
open (ALLRAW, $datafile);
my $dataheader=<ALLRAW>;
chomp($dataheader);
my %allraw;
my @allrawline;
while (<ALLRAW>) {
chomp;
  @allrawline = split /\t/, $_;
  $probeid = $allrawline[0];
  $line = join "\t", @allrawline;
  $allraw{$probeid}=$line;
}
close ALLRAW;

#combine this data in a new file (temp)
open (TEMP, ">$temp");
print TEMP "$dataheader\tMODEL\n";

my $value;
foreach $value (keys %allraw) {
#print $model{$value};
print TEMP $allraw{$value}."\t".$model{$value}."\n";

}
close TEMP;


 #extract UNCUT SNP data from temp datafile
 open (TEMP, "<$temp");
 open (UNCUTSNPS, $mnrfile);
 open (RAWUNCUT, ">$rawuncutdata");
 pullsnps (*TEMP, *UNCUTSNPS, *RAWUNCUT);
 close UNCUTSNPS;
 close RAWUNCUT;
 close TEMP;
print "UNCUTSNP data extracted\n";
#extract CUT SNP data from datafile)
open (TEMP, "<$temp");
open (CUTSNPS, $mprfile);
open (RAWCUT, ">$rawcutdata");
pullsnps (*TEMP, *CUTSNPS, *RAWCUT);
close CUTSNPS;
close RAWCUT;
close TEMP;
print "CUTSNP data extracted\n";

#run R quantile normalization in R on rawUNCUTdata and ouput both the 
# quantile normalized UNCUTdata results as well as the adjustments
open (QUANTUNCUT, ">$quantuncutdata");
open (ADJUSTS, ">$adjustdata");
system ("R --no-save < ./bin/Rquantilenorm.robust.R $rawuncutdata $quantuncutdata $adjustdata");
close QUANTUNCUT;
close ADJUSTS;

#open a data file to grab the header ($header) and count the number of cel files you are analysing ($numdatafields)
open (RAWCUT, "<$rawcutdata");
my $header=<RAWCUT>;
chomp $header;
my @headerfields=split (/\t/, $header);
my $numdatafields=@headerfields;
#print "number of cel files analysed: $numdatafields\n";
close RAWCUT;


#load the data files into a bunch of hashes 
#this is done one data column at a time using a while loop and the $fieldcount and $numdatafields variables
#(with common keys for ADJ and UNCUT)
my $fieldcount=1;
while ($fieldcount<$numdatafields) {
  #get the header for the fileoutput
  my $fileheader=$headerfields[$fieldcount];
  my $fileoutput=$fileheader."_quant.txt";
  open (FILEOUTPUT, ">$fileoutput");
  print FILEOUTPUT "probeID\tIntensity\n";
  
  #load in the raw UNCUT data, after chopping off the header
  #split the data and stuff into hash, 
  #key=first column=probeid ($uncutfields[0]) 
  #value= intensity data from whichever column is referred to by $fieldcount ($uncutfields[$fieldcount])
  #secondvalue=cutstatus of this probe
  open (RAWUNCUT, "<$rawuncutdata");
  <RAWUNCUT>;
  my %all;
  while (<RAWUNCUT> ) {
    my @fields =  split /\t/;
    my $data=$fields[$fieldcount];
    my $probeid=$fields[0];
    $all{$probeid}[0]=$data;
    $all{$probeid}[1]="UN";
  }
  close RAWUNCUT;
  print "$fieldcount: RAWUNCUT loaded\n";
  #load in the raw CUT data, after chopping off the header
  #same hash loading as above
  open (RAWCUT, "<$rawcutdata");
  <RAWCUT>;
  while (<RAWCUT> ) {
    my @fields =  split /\t/;
    my $data=$fields[$fieldcount];
    my $probeid=$fields[0];
    $all{$probeid}[0]=$data;
    $all{$probeid}[1]="CUT";
  }
  close RAWCUT;
  print "$fieldcount: RAWCUT loaded\n";
  #load in the adjustment values for each of the raw UNCUT data points, after chopping off the header
  #same hash loading as above
  open (ADJUSTS, "<$adjustdata"); 
   <ADJUSTS>;  
   my %ADJHASH;
   while (<ADJUSTS> ) {
     my @fields =  split /\t/;
     my $data=$fields[$fieldcount];
     my $probeid=$fields[0];
     $ADJHASH{$probeid}=$data;
   }
  close ADJUSTS;
  print "$fieldcount: ADJUSTS loaded\n";
  #load in the data you have already calculated for the uncut probes (you will print this out later after you have the cut values)
  open (QUANTUNCUT, "<$quantuncutdata"); 
  <QUANTUNCUT>;  
  my %QUANT;
  while (<QUANTUNCUT> ) {
    chomp;
    my @fields =  split /\t/;
    my $data=$fields[$fieldcount];
    my $probeid=$fields[0];
    $QUANT{$probeid}=$data;
  }
  close QUANTUNCUT;
  print "$fieldcount: QUANTILE normalized UNCUT data loaded\n";



  #sort the array of data  (%all) by the intensities and put in same order onto two arrays, one with the SNP names (@snparray), and the other with the cut/uncut statuses (@cutarray)
  my @snparray;
  my @cutarray;
  my $count=0;
  my $probeid;
  foreach $probeid (sort { $all{$a}[0] <=> $all{$b}[0] } keys %all) { 
    $all{$probeid}[2]=$count;
    $count=$count+1;
    push (@snparray, $probeid);
    push (@cutarray, $all{$probeid}[1]);
  }
  print "$fieldcount: HASH sorted, and pushed onto arrays\n";
  


  #this is the algorithm for finding the next highest/lowest uncut probe for each cut probe  
  my $totalprobes=@cutarray;
  my $probeindex;
  my $oversnp;
  my $undersnp;
  my $overtest;
  my $i;
  foreach $probeid (sort { $all{$a}[2] <=> $all{$b}[2] } keys %all) { 
    $probeindex=$all{$probeid}[2];
    if ($cutarray[$probeindex]=~ /CUT/) 
      { 
	#start the counter at the current index in the sort
	$i=$probeindex;
	#temporarily set the $oversnp value to be the same as the current probe, in case there is not a higher uncut snp
	$oversnp=$snparray[$i];
	#step through the array until you get to the top, and exit loop when you find an uncut probe, dump that value into the $oversnp variable
	while ($i<$totalprobes) {
	  if ($cutarray[$i]=~/UN/) {
	    $oversnp=$snparray[$i];
	    last;
	  }
	  $i=$i+1;
	}
	#reset the counter
	$i=$probeindex;
	#set the $undersnp to be the same as the $oversnp, in case there is not a lower uncut snp
	$undersnp=$oversnp;
	while (0<=$i) {
	  if ($cutarray[$i]=~/UN/) {
	    $undersnp=$snparray[$i];
	    last;
	  }
	  $i=$i-1;  
	}
	#in case there was no $oversnp found and is still set the same as the current probe, set it to be the same as the $undersnp!
	if ($oversnp=$snparray[$probeindex]) {
	  $oversnp=$undersnp;
	}
	my $overvalue=$all{$oversnp}[0];
	my $undervalue=$all{$undersnp}[0];
	my $currentvalue=$all{$probeid}[0];
	my $over=abs($overvalue-$currentvalue);
	my $under=abs($currentvalue-$undervalue);
	my $adjcutid;
	if ($over>=$under) {
	  $adjcutid=$undersnp;
	} else {
	  $adjcutid=$oversnp;
	}
	# using the probe id of the UNCUT data point which has the closest value to the CUT data point
	# grab the adjustment factor for that UNCUT probe and apply it to the raw CUT data point
	my $adjcut = $ADJHASH{$adjcutid}*$currentvalue;
	#round up result to one decimal place
	$adjcut =  sprintf("%.1f", $adjcut);
	#push these values onto the results hash (with the uncut data)
	$QUANT{$probeid}=$adjcut;
      }
  }
  #sort the quantile normalized data by the snpid (uncut and cut) and print out
  foreach $probeid (sort keys %QUANT) {
    print FILEOUTPUT "$probeid\t$QUANT{$probeid}\n";
  }
  print $fileheader." printed to file\n";
  print "$fieldcount: all adjusted values printed to file\n";
  $fieldcount=$fieldcount+1;
  close FILEOUTPUT;
}

exit;


#SUBSCRIPTS
sub pullsnps {
  my $data = shift;
  my $snps = shift;
  my $outputfile = shift;
  my $header=<$data>;
  my $key;
  print $outputfile $header;	 
  my %DATAHASH;
  while (<$data>) {
      my @elements = split /\t/;
      my $line = join "\t", @elements;
      $DATAHASH{$elements[0]}=$line;
    }
  close $data;
  <$snps>;
  while (<$snps>) {
    chomp $_;
    $key = $_;
    if (exists $DATAHASH{$key}) {
      print  $outputfile $DATAHASH{$key};
    }
    }
    close $snps;
  }

sub usage {
  print "Unknown option: @_\n" if ( @_ );
  print "usage: program [--datafile/-i absolute path to datafile] [--mnrfile/u absolute path to MNR probeid file] [--mprfile/c absolute path to MPR probeid file] [--modelfile/m absolute path to HapMap-derived model raw data file]  [--help|-h]\n";
  exit;
}
