#!/usr/bin/perl

### Compare pa_ggg_benchmark.mav, pa_ggg_benchmark.ray, pa_ggg_benchmark.vav, pa_ggg_benchmark.vsw

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
die "Usage: gggdiff.pl file1 file2\n" unless (($#ARGV+1) >= 2);

my @previous;
my @install;

##make sure both files are not empty
die "Warning: file $ARGV[0] is empty\n" if (-z $ARGV[0]);
die "Warning: file $ARGV[1] is empty\n" if (-z $ARGV[1]);

##check if input files are NOT .mav files
my(undef, undef, $ftype) = fileparse($ARGV[0], qr/\.[^.]*/);
die "This script only compares non .mav files.\n" if ($ftype eq ".mav");

##read in previous data
open (my $infile, $ARGV[0]) or die ("Can't open file $ARGV[0]");
my $firstrow = 0;
my $numhdr1 = 0;   #number of lines in header in the previous data file
my $numcol1 = 0;   #number of columns in the data in the 1st file

while (<$infile>) {
   chomp;
   my @aline = split(' ');
#   next if ($#aline == -1);   #skip empty lines

   if ($firstrow == 0) {
      $numhdr1 = $aline[0];   #number of lines in header, only works for *ray,*tav,*tsw,*vav,*vsw, *vav.* except *oof* and *vav.cew files
      $numcol1 = $aline[1];   #number of columns in the data
      $firstrow = 1;
   }

   push @previous, [@aline];
}
close ($infile);

##read in the new data from the second file
open ($infile, $ARGV[1]) or die ("Can't open file $ARGV[1]");
my $numhdr2 = 0;   #number of lines in header in the second file
my $numcol2 = 0;   #number of columns in the data in the second file
$firstrow = 0;
while (<$infile>) {
   chomp;
   my @aline = split(' ');
#   next if ($#aline == -1);   #skip empty lines

   if ($firstrow == 0) {
      $numhdr2 = $aline[0];   #number of lines in header
      $numcol2 = $aline[1];   #number of columns in the data
      $firstrow = 1;
   }

   push @install, [@aline];
}
close ($infile);

print "Comparing files $ARGV[0] and $ARGV[1]\n";

##compare header information of the 2 files 
if ($numhdr1 != $numhdr2) {
   print "Two files have different number of rows of headers: ";
   print "$numhdr1 vs. $numhdr2\n";
}
if ($numcol1 != $numcol2) {
   print "Two files have different number of columns of data: ";
   print "$numcol1 vs. $numcol2\n";
}
my $numbrow = @previous;       #total number of rows in the first file
my $numirow = $#install + 1;   #total number of rows in the second file
if ($numbrow != $numirow) {
   print "Two files have different total number of rows: ";
   print "$numbrow vs. $numirow\n";
}

##spectrum name row
my @specname1 = @{$previous[$numhdr1-1]};
my @specname2 = @{$install[$numhdr2-1]};

## find common columns that are in both files by comparing spectrum names
my ($ind1ref, $ind2ref) = getIndex(\@specname1, \@specname2);
my $num_comm_col = @$ind1ref;

##compare entire files row by row
my $count = 0;
my $maxdd = 0.;      #maximum difference
my $maxddfrac = 0.;  #maximum difference in percentage
my $maxrow;          #row where maximum difference occured
my $maxcol;          #col where maximum difference occured
my $ddfrac;
my $minhdrnum = min($numhdr1, $numhdr2);

##first, compare the headers row by row
my $equals;
my $hdrequal = 1;
for my $i (0 .. ($minhdrnum-2)) {   #common number rows of 2 headers except the spectrum names
    my @rowb = @{$previous[$i]};
    my @rowi = @{$install[$i]};

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
if ($numhdr1 > $numhdr2) {
   print "The first file has extra header lines:\n";
   for my $i (($numhdr2-1) .. $numhdr1-2) {
       print "> @{$previous[$i]}\n\n";
   } 
}
if ($numhdr2 > $numhdr1) {
   print "The second file has extra header lines:\n";
   for my $i (($numhdr1-1) .. $numhdr2-2) {
       print "> @{$install[$i]}\n\n";
   } 
}
##compare the last row of both headers (spectrum names)
my @rowb = @{$previous[$numhdr1-1]};
my @rowi = @{$install[$numhdr2-1]};
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
   print "Header line differs at spectrum name row.\n";
   findNewCols(\@rowb, \@rowi, $ARGV[0], $ARGV[1]);
   $hdrequal = 0;
}

##second, compare the numbers.  It checks up to the number of rows in the first
##file (i.e. if the 2nd file has more rows than the 1st, then only the beginning
##of $numbrow is checked)
my $ii = 0;
for my $i ($numhdr1 .. $numbrow-1) {
    if ($numirow >= $i) {
       my $m = $numhdr2 + $ii;
       my $cur_irow = $install[$m];
       my $icol = @$cur_irow;
       for my $j (0 .. $num_comm_col-1) {
           my $col1 = $ind1ref->[$j];
           my $col2 = $ind2ref->[$j];
#           if ($previous[$i][$col1] ne $install[$m][$col2]) {
           if (defined($previous[$i][$col1]) and defined($install[$m][$col2]) and $previous[$i][$col1] ne $install[$m][$col2]) {
              $count += 1;
              print "Difference found at row=".($i+1)." col=".($col1+1)." colname=$specname1[$col1]:\n";
              print "  benchmark=$previous[$i][$col1] installation=$install[$m][$col2]\n";
              #make sure these 2 are numbers first!!
              if ((looks_like_number($previous[$i][$col1])) && (looks_like_number($install[$m][$col2]))) {
                 my $dd = abs($previous[$i][$col1] - $install[$m][$col2]);
                 if ($previous[$i][$col1] !=0) {
                    $ddfrac = ($dd/abs($previous[$i][$col1]))*100.;
                 }
                 print "  absolute   difference=$dd\n";
                 print "  fractional difference=$ddfrac%\n";
                 if ($dd > $maxdd) {
                    $maxdd = $dd;
                    $maxddfrac = $ddfrac;
                    $maxrow = $i;
                    $maxcol = $col1;
                 }
              }
           }
       }
    } 
    $ii += 1;
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
   print "The maximum difference is $maxdd ($maxddfrac%)\n";
   print "\tat [row=".($maxrow+1)." col=".($maxcol+1)."] or [$specname1[1]=$previous[$maxrow][1],colname=$specname1[$maxcol]]\n";
   print "\n";
}

#find new columns in both files, if any
findNewCols(\@specname1, \@specname2, $ARGV[0], $ARGV[1]);

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
        my $match = grep {$_ eq $search} @array2;
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
        my $match = grep {$_ eq $search} @array1;
        if (!$match) {
           $ind2[$n] = $i;
           $n++;
        }
    }
    if ($n != 0) {
       print "File $file2 has new columns:\n";
       for my $i (0 .. ($n-1)) {
           print "-->New column: ".($ind2[$i]+1).", $array2ref->[$ind2[$i]]\n";
       }
       print "\n";
    }
}
