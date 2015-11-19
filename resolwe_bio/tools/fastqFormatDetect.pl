#!/usr/bin/env perl
# Author: Martin Dahlo / modified Jacques Dainat
#
# Usage:  perl scriptname.pl <infile> [-a -t <max seconds to search>]
# ex.
# perl scriptname.pl reads.fq
# perl scriptname.pl reads.fq -a
# perl scriptname.pl reads.fq -a -t 90

use warnings;
use strict;
use Getopt::Std;


=pod

Used to detect the format of a fastq file. It has 2 different modes, normal and advanced.

In the normal mode, it can only differentiate between Sanger/Illumina1.8+ and Solexa/Illumina1.3+/Illumina1.5+.

In the advanced mode, it will try to pinpoint exactly which scoring system is used. It will look at all quality scores until either

* There is only one scoring system that matches
* It has been running for more than a specified time
* It reaches the end of the file without either finding a single matching scoring syste, or running out of time before reaching the end of the file.


To run the program in normal mode, give it only the name of the fastq file:

perl scriptname.pl <infile>
Ex.
perl scriptname.pl myReads.fq


To run it in advanced mode with default timout of 60 seconds, specify -a:

perl scriptname.pl <infile> -a
Ex.
perl scriptname.pl myReads.fq -a


To run it in advanced mode with a custom timeout, specify -a and -t:

perl scriptname.pl <infile> -a -t <max seconds to search>
Ex.
perl scriptname.pl myReads.fq -a -t 600

 

The output from the program reports the interval in which qualities were observed (raw ascii numbers, not adjusting for any offsets or phred etc), and the scoring systems matcing these values.

Ex.

Time limit reached, observed qualities in range [37,67].
Possible matches:
Illumina 1.8+;Sanger 



Can easily be copy/pasted into any other script and altered to do other things than die when it has determined the format.

Pseudo code

* Open the fastq file
* Look at each quality ASCII char and convert it to a number
* Depending on if that number is above or below certain thresholds,
  determine the format. 

=cut

# REMINDER
my $remain="\n/!\\ We remain that ASCII code used by the differents score system overlap each other. We differentiate them only looking the part non-overlaping. So we assume that the reads analyzed should statistically contains at least one value of the non-overlapping part. 
Indeed, fastq file are enough large to contains each value possible of the quality score system.\n(fastqc do the same assumption and never reported any error)
";

# get filename
my $usage = <<EOF;
Usage:  perl scriptname.pl <infile> [-a -t <max seconds to search>]

-a		Advanced mode. Can be used to find exactly which scoring system it is.
-t		Set the max search time in seconds to be used when using -a. Default is 60.
EOF

my $fq = shift or die $usage;

# get flags, if any
getopts('at:');
our($opt_a, $opt_t);
my $adv = $opt_a;
my $time = $opt_t || 999999999;

# open the files
open FQ, "<", $fq or die $!;




if(!$adv){

	# initiate
	my @line;
	my $l;
	my $number;
	 
	 
	# go thorugh the file
	while(<FQ>){

		# if it is the line before the quality line
		if($_ =~ /^\+/){

			$l = <FQ>; # get the quality line
			chomp($l); # remove newline and whitespaces
			@line = split(//,$l); # divide in chars

			for(my $i = 0; $i <= $#line; $i++){ # for each char

				$number = ord($line[$i]); # get the number represented by the ascii char

				# check if it is sanger or illumina/solexa, based on the ASCII image at http://en.wikipedia.org/wiki/FASTQ_format#Encoding
				if($number > 74){ # if solexa/illumina
						die "This file looks like Solexa/Illumina1.3+/Illumina1.5+ format.\n"; # print result to terminal and die
				}elsif($number < 59){ # if sanger
					die "This file looks like Sanger/Illumina 1.8+ format.\n"; # print result to terminal and die
				}
			}
		}
	}

	die "Unconclusive, could be either Sanger/Illumina1.8+ or Solexa/Illumina1.3+/Illumina1.5+\n";
}






# if the user wants the advanced mode
if($adv){

	# initiate
	my @line;
	my $l;
	my $number;
	my $max = -100;
	my $min = 99999999;
	my $start = time();

	# scoring system definitions, according to http://en.wikipedia.org/wiki/FASTQ_format#Encoding 
	# Feel free to add more on your own, following the system of the ones already in here.
	my %systems = (	'Sanger', [33,126], 
			'Solexa', [59,126], 
			'Illumina 1.3+', [64,126], 
			'Illumina 1.5+', [66,126], 
			'Illumina 1.8+', [35,126]);

	my %infoDisplay = ( 'Sanger' => 'Phred+33 - It could be Sanger or Illumina 1.8+. Sorry this is the only case impossible to really differentiate !\nAnyway, you will be happy to learn that one or the other have exactly the same quality score system (Phred+33)'.
				'\nAs the character <I> is not present we could assume that is Sanger !',
                        	'Solexa' => 'Solexa+64 - It could be Solexa',
                        	'Illumina 1.3+' => 'Phred+64 - It could be Illumina 1.3+',
                        	'Illumina 1.5+'=> 'Phred+64 - It could be Illumina 1.5+',
				'Illumina 1.8+' => 'Phred+33 - It could be Sanger or Illumina 1.8+. Sorry this is the only case impossible to really differentiate !\nAnyway, you will be happy to learn that one or the other have exactly the same quality score system (Phred+33)'. 
                                '\nWe know that you really want to know exactly which Quality score it is... So, as the character <I> is present we could assume that is Illumina 1.8+ !',
				'last' => 'Phred+33 - It could be Sanger or Illumina 1.8+. Sorry this is the only case impossible to really differentiate !\nAnyway, you will be happy to learn that one or the other have exactly the same quality score system (Phred+33)'.
                                '\nWe know that you really want to know exactly which Quality score it is... but there is ASCII value over 41 we absolutly cannot differentiate them.');

 	my $nb_line =	`awk 'END {print NR}' $fq`;
	my $nb_read = $nb_line/4;
	my $startP=time;
	print "Your file contains $nb_read reads. The analysis could take a while.\n";
	my $nb_read_checked=0;	

	# go thorugh the file
	while(<FQ>){

		#Display progression
		if ((30 - (time - $startP)) < 0) {
        		my $done = ($nb_read_checked*100)/$nb_read;
			$done = sprintf ('%.0f', $done);
    			print "Progression : $done % processed.\n";
			$startP= time;

		}

		# if it is the line before the quality line
		if($_ =~ /^\+/){
			$nb_read_checked++;
			$l = <FQ>; # get the quality line
			chomp($l); # remove newline and whitespaces
			@line = split(//,$l); # divide in chars

			for(my $i = 0; $i <= $#line; $i++){ # for each char

				$number = ord($line[$i]); # get the number represented by the ascii char

				# check if the new number is larger or smaller than the previous records  
				if($number < $min){

					# update min and check how many systems are matching
					$min = $number;
					check($min, $max, \%systems, \%infoDisplay);
				}
				if($number > $max){

					# update max and check how many systems are matching
					$max = $number;
					check($min, $max, \%systems, \%infoDisplay);
				}

				# terminate if time is up
				if((time() - $start) >= $time){

					# print message to screen
					die "Time limit reached, observed qualities in range [$min,$max].\nPossible matches:\n".join("\n", check($min, $max, \%systems, \%infoDisplay))."\n".$remain."\n";
				}
			}
		}
	}

	# reached the end of the file without finding a definite answer, without running out of time
	die "Reached end of file, observed qualities in range [$min,$max].\nPossible matches:\n".join("\n", check($min, $max, \%systems, \%infoDisplay))."\n".$remain."\n";

}




###subroutines

# check how many scoring systems are matching the current max min values
sub check{

	# get arguments
	my ($min, $max, $systems, $infoDisplay) = @_;

	# init
	my @matching;

	# check available systems
	foreach my $key (keys %{$systems}){

		# is it a match?
		if( ($min >= $systems->{$key}[0]) && ($max <= $systems->{$key}[1]) ){

			# save matching systems
			my $messageToDisplay = $infoDisplay->{$key};
			push(@matching, $messageToDisplay);

		}

	}


	# check if only one system matched
	if($#matching == 0){

		# print message to screen
		die "Only one possible match, observed qualities in range [$min,$max]:\n$matching[0]\n";

	}

	# If still not dtermined
        if($#matching >= 1){
		@matching=();
		

			if($min >= $systems->{'Illumina 1.5+'}[0]){
				my $messageToDisplay = $infoDisplay->{'Illumina 1.5+'};
                        	push(@matching, $messageToDisplay);
			}
			elsif($min >= $systems->{'Illumina 1.3+'}[0]){ 
                                my $messageToDisplay = $infoDisplay->{'Illumina 1.3+'};
                                push(@matching, $messageToDisplay);
                        }
			elsif($min >= $systems->{'Solexa'}[0]){ 
                                my $messageToDisplay = $infoDisplay->{'Solexa'};
                                push(@matching, $messageToDisplay);
                        }
			else{ #could be Illumina 1.8+ or Sanger
				if($max == $systems->{'Sanger'}[1]){
					my $messageToDisplay = $infoDisplay->{'Sanger'};
					push(@matching, $messageToDisplay);
				}
				elsif($max == $systems->{'Illumina 1.8+'}[1]){
					my $messageToDisplay = $infoDisplay->{'Illumina 1.8+'};
                                        push(@matching, $messageToDisplay);
				}
				else{
					my $messageToDisplay = $infoDisplay->{'last'};
                                        push(@matching, $messageToDisplay);
				}
			}
			
        }

	# return all matching systems
	return @matching;
}