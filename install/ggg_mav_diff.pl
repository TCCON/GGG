#!/usr/bin/perl

### Compare pa_ggg_benchmark.mav files

use strict;
use warnings;

use Scalar::Util qw(looks_like_number);
use File::Basename;
#use Term::ANSIColor qw(:constants);
#$Term::ANSIColor::AUTORESET = 1;
use List::Util qw[min max];

#file1=$ARGV[0]
#file2=$ARGV[1]
#Example files:
#my $file1 = "$INSTALLDIR/install/benchmark_results/pa_ggg_benchmark.mav";
#my $file2 = "$INSTALLDIR/install/current_results/pa_ggg_benchmark.mav";

##check syntax
die "Usage: ggg_mav_diff.pl file1 file2\n" unless (($#ARGV+1) >= 2);

my @previous;
my @install;

##make sure both files are not empty
die "Warning: file $ARGV[0] is empty\n" if (-z $ARGV[0]);
die "Warning: file $ARGV[1] is empty\n" if (-z $ARGV[1]);

##check if input files are .mav files
my(undef, undef, $ftype) = fileparse($ARGV[0], qr/\.[^.]*/);
die "This script only compares .mav files?\n" unless ($ftype eq ".mav");

##read in previous data
open (my $infile, $ARGV[0]) or die ("Can't open file $ARGV[0]");
my $firstrow = 0;
my $numhdr1 = 0;   #number of lines in the header in the 1st file
my $numcol1 = 0;   #number of columns in the data in the 1st file
my $numrow1 = 0;   #number of rows in the data in the 1st file

while (<$infile>) {
   chomp;
   my @aline = split(' ');
   next if ($#aline == -1);   #skip empty lines

   if ($firstrow == 2) {
      $numhdr1 = $aline[0];   #number of lines in 1st header in the 1st file
      $numcol1 = $aline[1];   #number of columns in the 1st part of the data in the 1st file
      $numrow1 = $aline[2];   #number of rows in the 1st part of the data in the 1st file
   }

   push @previous, [@aline];
   $firstrow += 1;
}
close ($infile);

##read in the new data from the second file
open ($infile, $ARGV[1]) or die ("Can't open file $ARGV[1]");
my $numhdr2 = 0;   #number of lines in the header in the second file
my $numcol2 = 0;   #number of columns in the data in the second file
my $numrow2 = 0;   #number of rows in the data in the second file
$firstrow = 0;
while (<$infile>) {
   chomp;
   my @aline = split(' ');
   next if ($#aline == -1);   #skip empty lines

   if ($firstrow == 2) {
      $numhdr2 = $aline[0];   #number of lines in 1st header in the second file
      $numcol2 = $aline[1];   #number of columns in the 1st part of the data in the second file
      $numrow2 = $aline[2];   #number of rows in the 1st part of the data in the second file
   }
   push @install, [@aline];
   $firstrow += 1;
}
close ($infile);

print "Comparing files $ARGV[0] and $ARGV[1]\n";

##compare the 1st part of the data row by row
##compare number of rows in the headers of the 2 files 
if ($numhdr1 != $numhdr2) {
   print "Two files have differnent number of rows of headers in the 1st spectrum: ";
   print "$numhdr1 vs. $numhdr2\n";
}
if ($numcol1 != $numcol2) {
   print "Two files have differnent number of columns of data in the 1st spectrum: ";
   print "$numcol1 vs. $numcol2\n";
   print "Only the first ".min($numcol1, $numcol2)." columns will be compared\n";
}
if ($numrow1 != $numrow2) {
   print "Two files have differnent number of rows of data in the 1st spectrum: ";
   print "$numrow1 vs. $numrow2\n";
}
my $numbrow = @previous;       #total number of rows in the first file
my $numirow = $#install + 1;   #total number of rows in the second file
if ($numbrow != $numirow) {
   print "Two files have differnent total number of rows: ";
   print "$numbrow vs. $numirow\n";
}

##let's figure out how many spectrum there are in the file
my $numbrow1part = $numhdr1 + $numrow1 + 1;
my $numspectrum = ($numbrow-1)/$numbrow1part;

##init
my @specname;
my $count = 0;
my $maxdd = 0.;      #maximum difference
my $maxddfrac = 0.;  #maximum difference in percentage
my $maxrow;          #row where maximum difference occured
my $maxcol;          #col where maximum difference occured
my $ddfrac;
my $equals = 1;
my $hdrequal = 1;
my $minnumcol = min($numcol1, $numcol2);
my $thespectrum = 0;

##compare the first line of the two files
my @line1file1 = @{$previous[0]};
my @line1file2 = @{$install[0]};
if (@line1file1 != @line1file2) {
   $equals = 0;
} else {
   foreach (my $j = 0; $j < @line1file1; $j++) {
      if ($line1file1[$j] ne $line1file2[$j]) {
         $equals = 0;
         last;
      }
   }
}
if ($equals == 0) {
   print "Header line differs at row 1:\n";
   print "< @line1file1\n";
   print "---\n";
   print "> @line1file2\n\n";
   $hdrequal = 0;
}

##compare the rest
for my $ispec (0 .. $numspectrum-1) {
    ##first, compare the headers
    my $hdrs = 1 + $ispec * $numbrow1part;
    my $hdre = $hdrs + $numhdr1;
    for my $i ($hdrs .. $hdre) {
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
          print "Header line differs at row ".($i+1).":\n";
          print "< @rowb\n";
          print "---\n";
          print "> @rowi\n\n";
          $hdrequal = 0;
       }
    }
    ##second, compare the numbers
    @specname = @{$previous[$hdre]};
    my $dats = $hdre + 1;
    my $date = $dats + ($numrow1-1);
    for my $i ($dats .. $date) {
        #my $cur_brow = $previous[$i];
        #my $bcol = @$cur_brow;
        if ($numirow >= $i) {
           #my $cur_irow = $install[$i];
           #my $icol = @$cur_irow;
           for my $j (0 .. $minnumcol-1) {
               if ($previous[$i][$j] ne $install[$i][$j]) {
                  $count += 1;
                  print "Difference found at row=".($i+1)." col=".($j+1)." specname=$specname[$j]:\n";
                  print "  benchmark=$previous[$i][$j] installation=$install[$i][$j]\n";
                  #make sure these 2 are numbers first!!
                  if ((looks_like_number($previous[$i][$j])) && (looks_like_number($install[$i][$j]))) {
                     my $dd = abs($previous[$i][$j] - $install[$i][$j]);
                     if ($previous[$i][$j] !=0) {
                        $ddfrac = ($dd/abs($previous[$i][$j]))*100.;
                     } else {
                        $ddfrac = ($dd/abs($install[$i][$j]))*100.;
                     }
                     print "  absolute   difference=$dd\n";
                     print "  fractional difference=$ddfrac%\n";
                     if ($dd > $maxdd) {
                        $maxdd = $dd;
                        $maxddfrac = $ddfrac;
                        $maxrow = $i;
                        $maxcol = $j;
                        $thespectrum = $ispec;
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
      my $hdrs = 1 + $thespectrum* $numbrow1part;
      my $hdre = $hdrs + $numhdr1;
      @specname = @{$previous[$hdre]};
      print "The maximum difference is $maxdd ($maxddfrac%)\n";
      print "\tat [row=".($maxrow+1)." col=".($maxcol+1)."] or [$specname[0]=$previous[$maxrow][0], specname=$specname[$maxcol]]\n\n";
   }
}
print "\n";

