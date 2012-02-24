#!/usr/bin/perl

### Compare pa_ggg_benchmark.mav, pa_ggg_benchmark.ray, pa_ggg_benchmark.rav, pa_ggg_benchmark.vsw

use strict;
use warnings;

use Scalar::Util qw(looks_like_number);
use File::Basename;
#use Term::ANSIColor qw(:constants);
#$Term::ANSIColor::AUTORESET = 1;

#file1=$ARGV[0]
#file2=$ARGV[1]
#Example files:
#my $file1 = "$INSTALLDIR/install/benchmark_results/pa_ggg_benchmark.ray";
#my $file2 = "$INSTALLDIR/install/current_results/pa_ggg_benchmark.ray";

##check syntax
die "Usage: gggdiff.pl file1 file2\n" unless (($#ARGV+1) >= 2);

my @previous;
my @install;

##make sure both files are not empty
die "Warning: file $ARGV[0] is empty\n" if (-z $ARGV[0]);
die "Warning: file $ARGV[1] is empty\n" if (-z $ARGV[1]);

##check if input files are NOT .mav files
my(undef, undef, $ftype) = fileparse($ARGV[0], qr/\.[^.]*/);
die "This script only compares non .mav files?\n" if ($ftype eq ".mav");

##read in previous data
open (my $infile, $ARGV[0]) or die ("Can't open file $ARGV[0]");
my $firstrow = 0;
my $numhdr1 = 0;   #number of lines in header in the previous data file
while (<$infile>) {
   chomp;
   my @aline = split(' ');
   next if ($#aline == -1);   #skip empty lines

   if ($firstrow == 0) {
      $numhdr1 = $aline[0];   #number of lines in header, only works for *ray,*tav,*tsw,*vav,*vsw, *vav.* except *oof* and *vav.cew files
      $firstrow = 1;
   }

   push @previous, [@aline];
}
close ($infile);

##read in the new data from the second file
open ($infile, $ARGV[1]) or die ("Can't open file $ARGV[1]");
my $numhdr2 = 0;   #number of lines in header in the second file
$firstrow = 0;
while (<$infile>) {
   chomp;
   my @aline = split(' ');
   next if ($#aline == -1);   #skip empty lines

   if ($firstrow == 0) {
      $numhdr2 = $aline[0];   #number of lines in header
      $firstrow = 1;
   }

   push @install, [@aline];
}
close ($infile);

print "Comparing files $ARGV[0] and $ARGV[1]\n";

##compare number of rows in the headers of the 2 files 
if ($numhdr1 != $numhdr2) {
   print "Two files have differnent number of rows of headers: ";
   print "$numhdr1 vs. $numhdr2\n";
}
my $numbrow = @previous;       #total number of rows in the first file
my $numirow = $#install + 1;   #total number of rows in the second file
if ($numbrow != $numirow) {
   print "Two files have differnent total number of rows: ";
   print "$numbrow vs. $numirow\n";
}

##spectrum name row
my @specname= @{$previous[$numhdr1-1]};

##compare entire files row by row
my $count = 0;
my $maxdd = 0.;      #maximum difference
my $maxddfrac = 0.;  #maximum difference in percentage
my $maxrow;          #row where maximum difference occured
my $maxcol;          #col where maximum difference occured
my $ddfrac;
##first, compare the headers
my $equals;
my $hdrequal = 1;
for my $i (0 .. $numhdr1-1) {
    my @rowb = @{$previous[$i]};
    my @rowi = @{$install[$i]};
    my $rowbss = $previous[$i];
    my $rowbsscol = @$rowbss;

    if (@rowb != @rowi) {
       $equals = 0;
    } else {
       $equals = 1;
        foreach (my $j = 0; $j < @rowb; $j++) {
          if ($rowb[$j] ne $rowi[$j]) {
             $equals = 0;
             last;
          }
       }
    }
    if ($equals == 0) {
       print "Header line differs:\n";
       print "< @rowb\n";
       print "---\n";
       print "> @rowi\n\n";
       $hdrequal = 0;
    }
}
if ($numhdr2 > $numhdr1) {
   print "The second file has extra header lines:\n";
   for my $i ($numhdr1 .. $numhdr2-1) {
       print "> @{$install[$i]}\n\n";
   } 
}

##second, compare the numbers, against previous data format such as number of headers
##and number of rows of spectrum
for my $i ($numhdr1 .. $numbrow-1) {
    my $cur_brow = $previous[$i];
    my $bcol = @$cur_brow;
    if ($numirow >= $i) {
       my $cur_irow = $install[$i];
       my $icol = @$cur_irow;
       if ($bcol != $icol) {
          $count += 1;
          print "Number of columns is different in row $i\n";    
       } else {
          for my $j (0 .. $bcol-1) {
              if ($previous[$i][$j] ne $install[$i][$j]) {
                 $count += 1;
                 print "Difference found at row=".($i+1)." col=".($j+1)." specname=$specname[$j]:\n";
                 print "  benchmark=$previous[$i][$j] installation=$install[$i][$j]\n";
                 #make sure these 2 are numbers first!!
                 if ((looks_like_number($previous[$i][$j])) && (looks_like_number($install[$i][$j]))) {
                    my $dd = abs($previous[$i][$j] - $install[$i][$j]);
                    if ($previous[$i][$j] !=0) {
                       $ddfrac = ($dd/abs($previous[$i][$j]))*100.;
                    }
                    print "  absolute   difference=$dd\n";
                    print "  fractional difference=$ddfrac%\n";
                    if ($dd > $maxdd) {
                       $maxdd = $dd;
                       $maxddfrac = $ddfrac;
                       $maxrow = $i;
                       $maxcol = $j;
                    }
                 }
              }
          }
       } 
    }
}

##print out summary
if (($count == 0) && ($hdrequal== 1)) {
   print "\nTwo files are the ****SAME****\n";
} else {
   print "\nTwo files are ****DIFFERENT****\n";
   if ($count > 0) {
      print "The maximum difference is $maxdd ($maxddfrac%)\n";
      print "\tat [row=".($maxrow+1)." col=".($maxcol+1)."] or [$specname[0]=$previous[$maxrow][0],specname=$specname[$maxcol]]\n\n";
   }
}
print "\n";

