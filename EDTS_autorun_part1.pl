#! /usr/bin/perl 
$|++;

# Original Author(s):
#     C.Y. Lin, E.I. Izgorodina, D. Brittain, and K. Zhang
#
# Original script can be found at: 
#     https://rsc.anu.edu.au/~mcoote/scripts.php
# 
# When using this codes please using the following citation: 
#     E. I. Izgorodina , C.Y. Lin, M. L. Coote. Phys. Chem. Chem. Phys. 9, 2507 (2007)
#
######################################################################################
#                                                                          
# EDTS_autorun.              PART 1/3
#
# Author ofchanges: Samuel James Higgs: s.higgs@surrey.ac.uk: 30/09/2021
#
# Based on ConfSearch.pl but modified to use SLURM batch job control system.                                                         
#                                                                          
# This code is for performing a conformer search, implementing the EDTS method. This method is designed to identify the lowest-energy           
# conformation of a molecule (or at least a structure within 1 kcal mol−1 of it) with a very high degree of reliability, without exhaustively     
# searching all conformations                                             
# 
# NOTE
# Use Confmaker program first to generate the Gaussian input files for the full conformational space of a molecule 
#                                                                          
######################################################################################
#                                                                          
# USAGE
#
#     $ module load perl/5.24.0-GCC-5.4.0-2.26
#
# # to execute the Perlcode, first make the script executable from the directory containing EDTS perl scripts.                       
#     $ chmod +x *.pl                                               
#
# EXE
#
#     $ ConfSearch $mol [$atom1 $atom2] [$gooddist] [$toldist]
#
# NOTE [] are for transition state optimisation, if r12 is fixed, this can be used to help determine the optimised geometry is valid.                                        #
#                                                                          
#  REQUIREMENTS 
#
#     1] Input file format "$mol.$dih.com" -> ensure the use of the period (.) as opposed to underscores (_) or hyphens (-) 
#     2] $mol.list in the input directory.  
#
#  DEFAULT VALUES
#     
#     $EC1 = 3 
#     $EC2 = 4
#     $Nmax = 5
#                                                                          
######################################################################################
#
# COMMENTS FROM AUTHORS
#                                                                          
#  Not yet deal with Large system where the round2/3 maybe 1000            
#  All the job.list is printout in CF-*                                    
#  CF-$mol.round1 (if 2nd lowest < 4kJ) -> CF-$mol.round2 -> CF-$mol.round3-$i.tosub
#                (if 2nd lowest > 4kJ)                  -> CF-$mol.round3-$i.tosub
#  A few other useful files:
#  1. CF-$mol.round$x.uh/lh (upper/lower half with the energy difference)
#  2. CF-$mol.*.window (those energies within 4kJ and the first Nmax with their energy difference)
#  3. CF-$mol.done (At every stage, the done file will tell you what has been done)
#  4. At the end of the program, you can look into CF-$mol.done.*, which tell you the all Energy difference.
#
######################################################################################

#use Term::ANSIColor;
#print color 'blue';
printoutput "###################################################################################\n"; 
printoutput "#                                        EDTS                                     #\n";
printoutput "###################################################################################\n\n"; 
#print color 'reset';

# Defaults
$Nmax = 5;  # Default
$EC1 = 3;  # Default 3kJ/mol
$EC2 = 4;  # Default 4kJ/mol

# Input arguments 
$mol=$ARGV[0];  

# Optional arguments (all must be given if used)
$atom1=$ARGV[1];
$atom2=$ARGV[2];

if ($ARGV[1]){
    if (!$ARGV[3] ){
        printoutput "No appropriate bond length specified, 1.33 A will be used.\n";
        printoutput "No appropriate bond length tolerance specified, 0.00 A will be used.\n";
        $gooddist=1.33;
        $toldist=0.00;
    }
    else{
        $gooddist=$ARGV[3];
        printoutput "Bond length specified at $gooddist A.\n";
        if (!$ARGV[4]){
            printoutput "No appropriate bond length tolerance specified, 0.00 A has been used.\n";
            $toldist=0.00;
        }
        else{
            $toldist=$ARGV[4];
            printoutput "Bond length tolerance specified at $toldist A.\n";
        }
    }

}
chomp($mol);

# remove any previous log file called $mol.log
if (-e "$DataDir/$mol.log") {
    system("rm $DataDir/$mol.log");
}

printoutput "\nConformer search for $mol.\n";

use Cwd 'abs_path';
use File::Basename;

# path to directory containing all edts scripts
my $CmdFile = abs_path(__FILE__);
my $CmdDir = dirname($CmdFile);
#$CmdDir = getcwd; 

# path to gaussian file directory
use Cwd;
$DataDir = cwd;
$DataDir = "$DataDir/g16";

printoutput "$CmdDir\n";
printoutput "$DataDir\n";

if (-e "$DataDir/log.txt") {
    system("rm $DataDir/log.txt");
}

###################################################################################
#                                       Round1                                    #
###################################################################################
#print color 'green';
printoutput "###################################################################################\n"; 
printoutput "#                                       Round1                                    #\n";
printoutput "###################################################################################\n"; 
#print color 'reset';

# create a list of individual rotations for first round of optimisation
$file = "$DataDir/CF-$mol.round1";
#print "$file\n";

if (-e $file) {
    printoutput "Old file found, removing $file\n";
    system("rm $file");
}

open(comlist,'<', "$DataDir/$mol.list") or die $!;
while ($dih = <comlist>){
    #$dih = $_;  # make code clearer in above assignment
    open(fr, '>>', "$file") or die $!;
    if (length($dih) <4 || (length($dih)<6 && $dih=~/a1/) ){
        print fr "$mol.$dih";
    }
    close fr;
}
close comlist;

# submits array of jobs based on input file list $mol.round1 and dependent job which 
# performs round 2 processing
my @args = ("$CmdDir/subarrayjob", "$file", "$CmdDir", "$DataDir", "$CmdDir/EDTS_autorun_part2.pl $mol $CmdDir $DataDir $atom1 $atom2 $gooddist $toldist");
exec ("/bin/bash", @args) == 0 or die "system @args failed: $?"

# subroutine to echo STDOUT to mol.log file
# change print to printoutput above

sub printoutput {
    my $text = @_;
    
    #print to STDOUT
    print $text;

    # append text to logfile called mol.log
    # open/print/close approach ensures log text flush to disk/ viewable as script runs through
    open(mollog,'>>', "$DataDir/$mol.log") or die $!;
    print mollog $text;
    close mollog;
}


__END__

