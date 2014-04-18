#!/usr/bin/perl

### Compare pa_ggg_benchmark.mav files

use strict;
use warnings;
 
use Scalar::Util qw(looks_like_number);
use File::Basename;
use List::Util qw[min max];
use List::Util 'first';

#subroutine declaration
sub getIndex;
sub findNewCols;

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
   print "Two files have different number of rows of headers in the 1st spectrum: ";
   print "$numhdr1 vs. $numhdr2\n";
}
if ($numcol1 != $numcol2) {
   print "Two files have different number of columns of data in the 1st spectrum: ";
   print "$numcol1 vs. $numcol2\n";
}
if ($numrow1 != $numrow2) {
   print "Two files have different number of rows of data in the 1st spectrum: ";
   print "$numrow1 vs. $numrow2\n";
}
my $numbrow = @previous;       #total number of rows in the first file
my $numirow = $#install + 1;   #total number of rows in the second file
if ($numbrow != $numirow) {
   print "Two files have different total number of rows: ";
   print "$numbrow vs. $numirow\n";
}

##let's figure out how many spectrum there are in the file
my $numbrow1part = $numhdr1 + $numrow1 + 1;
my $numspectrum1 = ($numbrow-1)/$numbrow1part;
my $numbrow2part = $numhdr2 + $numrow2 + 1;

##init
my @specname1;
my @specname2;
my $count = 0;
my $maxdd = 0.;      #maximum difference
my $maxddfrac = 0.;  #maximum difference in percentage
my $maxrow;          #row where maximum difference occured
my $maxcol;          #col where maximum difference occured
my $ddfrac;
my $equals = 1;
my $hdrequal = 1;
my $minhdrnum = min($numhdr1, $numhdr2);
my $thespectrum = 0;

##compare the first line of the two files
my @line1file1 = @{$previous[0]};
my @line1file2 = @{$install[0]};
my $join1 = join('  ', @line1file1);
my $join2 = join('  ', @line1file2);

if ($join1 ne $join2) { 
   print "Header line differs at row 1:\n";
   print "< @line1file1\n";
   print "---\n";
   print "> @line1file2\n\n";
   $equals = 0;
   $hdrequal = 0;
}

##compare the rest
for my $ispec (0 .. $numspectrum1-1) {
    ##first, compare the headers
    my $hdrs1 = 1 + $ispec * $numbrow1part;
    my $hdre1 = $hdrs1 + $numhdr1;
    my $hdrs2 = 1 + $ispec * $numbrow2part;
    my $hdre2 = $hdrs2 + $numhdr2;
    my $ii = 0;
    for my $i ($hdrs1 .. ($hdrs1+$minhdrnum-1)) { #common num. rows of 2 headers except the last one (gas names)
        my $m = $ii + $hdrs2;
        my @rowb = @{$previous[$i]};
        my @rowi = @{$install[$m]};
        
        my($filenameb, $directoriesb) = fileparse(@rowb);
        my($filenamei, $directoriesi) = fileparse(@rowi);

        if (@rowb != @rowi) {
           $equals = 0;
        } elsif ($filenameb eq $filenamei) {
           $equals = 1;
           foreach (my $j = 0; $j < @rowb; $j++) {
              if ($rowb[$j] ne $rowi[$j]) {
                 $hdrequal = 0;
                 last;
              }
           }
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
       $ii += 1;
    }
    if ($numhdr1 > $minhdrnum) {        #1st file has more header rows
       print "The first file has extra header lines:\n";
       for my $i (($hdrs1+$minhdrnum) .. ($hdre1-1)) {
           my @rowb = @{$previous[$i]};
           print "< @rowb\n\n";
       }
    } elsif ($numhdr2 > $minhdrnum) {   #2nd file has more header rows
       print "The second file has extra header lines:\n";
       for my $i (($hdrs2+$minhdrnum) .. ($hdre2-1)) {
           my @rowi = @{$install[$i]};
           print "> @rowi\n\n";
       }
    }
    ##compare the last row of both headers (gas names)
    my @rowb = @{$previous[$hdre1]};
    my @rowi = @{$install[$hdre2]};
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
       print "Header line differs at gas name row.\n";
       #find new columns in both files, if any
       findNewCols(\@rowb, \@rowi, $ARGV[0], $ARGV[1]);
       $hdrequal = 0;
    }
    ##second, compare the numbers
    @specname1 = @{$previous[$hdre1]};
    @specname2 = @{$install[$hdre2]};

    ## find common columns that are in both files by comparing window names
    my ($ind1ref, $ind2ref) = getIndex(\@specname1, \@specname2);
    my $num_comm_col = @$ind1ref;
    #
    my $dats1 = $hdre1 + 1;
    my $date1 = $dats1 + ($numrow1-1);
    my $dats2 = $hdre2 + 1;
    my $date2 = $dats2 + ($numrow2-1);
    $ii = 0;
    for my $i ($dats1 .. $date1) {
        my $m = $ii + $dats2;
        if ($numirow >= $i) {
           for my $j (0 .. $num_comm_col-1) {
               my $col1 = $ind1ref->[$j];
               my $col2 = $ind2ref->[$j];
               if ($previous[$i][$col1] ne $install[$m][$col2]) {
                  $count += 1;
                  print "Difference found at row=".($i+1)." col=".($col1+1)." specname=$specname1[$col1]:\n";
                  print "  benchmark=$previous[$i][$col1] installation=$install[$m][$col2]\n";
                  #make sure these 2 are numbers first!!
                  if ((looks_like_number($previous[$i][$col1])) && (looks_like_number($install[$m][$col2]))) {
                     my $dd = abs($previous[$i][$col1] - $install[$m][$col2]);
                     if ($previous[$i][$col1] !=0) {
                        $ddfrac = ($dd/abs($previous[$i][$col1]))*100.;
                     } else {
                        $ddfrac = ($dd/abs($install[$m][$col2]))*100.;
                     }
                     print "  absolute   difference=$dd\n";
                     print "  fractional difference=$ddfrac%\n";
                     if ($dd > $maxdd) {
                        $maxdd = $dd;
                        $maxddfrac = $ddfrac;
                        $maxrow = $i;
                        $maxcol = $col1;
                        $thespectrum = $ispec;
                     }
                  }
                }
           }
        } 
        $ii += 1;
    }
}

##print out summary
if ($count == 0) {
   if ($hdrequal == 1) {
      print "\nTwo files are the ****SAME****\n\n";
   } else {
      print "\nTwo files are ****DIFFERENT**** only in headers\n";
   }
} else {
   print "\nTwo files are ****DIFFERENT****\n";
   my $hdrs = 1 + $thespectrum* $numbrow1part;
   my $hdre = $hdrs + $numhdr1;
   @specname1 = @{$previous[$hdre]};
   print "The maximum difference is $maxdd ($maxddfrac%)\n";
   print "\tat [row=".($maxrow+1)." col=".($maxcol+1)."] or [$specname1[0]=$previous[$maxrow][0], specname=$specname1[$maxcol]]\n\n";
}
print "\n";

##compare two string arrays and return their index where they are the same
sub getIndex() {
    my @ind1;
    my @ind2;
    my $index;
    my $n = 0;

    my ($array1ref, $array2ref) = @_;
    my @array2 = @$array2ref;
    my $array1len = @$array1ref;

    for my $i (0 .. ($array1len - 1)) {
        my $search = $array1ref->[$i];
        my $index = first { $array2[$_] eq $search } 0 .. $#array2;
        if (defined $index) {
           $ind1[$n] = $i;
           $ind2[$n] = $index;
           $n++;
        }
    }
    return(\@ind1, \@ind2);
}

##find new columns in both files, if any
sub findNewCols() {
    my @ind1;
    my @ind2;
    my $n = 0;

    my ($array1ref, $array2ref, $file1, $file2) = @_;
    my @array1 = @$array1ref;
    my @array2 = @$array2ref;
    my $array1len = @$array1ref;
    my $array2len = @$array2ref;

    #search for new columns in the first file
    for my $i (0 .. ($array1len - 1)) {
        my $search = $array1ref->[$i];
        my $match = grep {/$search/} @array2;
        if (!$match) {
           $ind1[$n] = $i;
           $n++;
        }
    }
    if ($n != 0) {
       print "File $file1 has new columns:\n";
       for my $i (0 .. ($n-1)) {
           print "New column: ".($ind1[$i]+1).", $array1ref->[$ind1[$i]]\n";
       }
       print "\n";
    }

    #search for new columns in the second file
    $n = 0;
    for my $i (0 .. ($array2len - 1)) {
        my $search = $array2ref->[$i];
        my $match = grep {/$search/} @array1;
        if (!$match) {
           $ind2[$n] = $i;
           $n++;
        }
    }
    if ($n != 0) {
       print "File $file2 has new columns:\n";
       for my $i (0 .. ($n-1)) {
           print "New column: ".($ind2[$i]+1).", $array2ref->[$ind2[$i]]\n";
       }
       print "\n";
    }
}
