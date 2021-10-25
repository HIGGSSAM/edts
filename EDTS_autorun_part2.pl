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
# EDTS_autorun.              PART 2/3
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

# Defaults
$Nmax = 5;  # Default
$EC1 = 3;  # Default 3kJ/mol
$EC2 = 4;  # Default 4kJ/mol

# Input arguments
$mol = $ARGV[0];  
$CmdDir = $ARGV[1];
$DataDir = $ARGV[2];

# Optional arguments (all must be given if used)
$atom1 = $ARGV[3];
$atom2 = $ARGV[4];
$gooddist = $ARGV[5];
$toldist = $ARGV[6];

chomp($mol);

$log = "$DataDir/$mol.log";

printoutput($log, "\n###################################################################################\n"); 
printoutput($log, "#                                        EDTS                                     #\n");
printoutput($log, "###################################################################################\n\n"); 

printoutput($log, "Conformer search for $mol.\n");
printoutput($log, "Data directory: $DataDir\n\n");

printoutput($log, "###################################################################################\n");
printoutput($log, "#                                       Round2                                    #\n");
printoutput($log, "###################################################################################\n");

# initialise TSchecks result file
if (-e "$DataDir/CF-$mol.TScheck.results"){
    system("/bin/rm $DataDir/CF-$mol.TScheck.results");
}

system("cat $DataDir/CF-$mol.round1 | sort | uniq > $DataDir/CF-$mol.done");
open done, "$DataDir/CF-$mol.done" or die $!;
while ($done = <done>){
    chomp($done);
    foreach $y ($done){
        my @args = ("$y", "$mol", "$DataDir", "$CmdDir","$log", "$atom1", "$atom2", "$goodist", "$toldist");
        #print "@args";
        system($^X, "$CmdDir/lib/edts_TScheck.pl", @args) == 0 or die "system @args failed: $?";
        }
    }
close done;

$errcount=0;
$lnecount=0;
# count how many TScheck errors detected
open(errcount, '<', "$DataDir/CF-$mol.TScheck.results") or die $!;
while ($line = <errcount>){
    $linecnt++;
    if ($line=~ /Error/){
            $errcount++;
    }
}
close errcound;
print ("\nErrcount=$errcount of $linecnt\n");
# terminate if more than 20% of files errored
$errpercent=$errcount/$linecnt*100;
#print ("$errpercent");
if ($errpercent > 20) {
    printf "\n*** More than %3.2f % ***\n",$errpercent;
    print "*** Too many Gaussian errors to continue ***\n\n";
    exit;
}


# Rank the energies of optimised conformations
my @args = ("CF-$mol.round1", "$EC1","$DataDir", "$CmdDir", "$log", "$Nmax", "$atom1", "$atom2", "$goodist", "$toldist");
system($^X, "$CmdDir/lib/edts_opteng.pl", @args) == 0 or die "system @args failed: $?";

# count how many energy like conformations found in first optimisation run
open(lengthx, '<', "$DataDir/CF-$mol.round1.lh") or die $!;
$engcount++ while (<lengthx>);
close lengthx;

my @args = ("CF-$mol.round1.uh", "$DataDir");
system($^X, "$CmdDir/lib/edts_squeeze.pl", @args) == 0 or die "system @args failed: $?";

# round.ll contains list of energy like conformations from first optimisation run

# Yes there are more than one energy like conformations (case $1 ne 1) so for
# half of the space perform full conformational search
#if ($l ne "1"){
if ($engcount ne "1"){

    printoutput($log, ".....Round 2 performing full conformational search on half the space...\n\n");

    # for half the space (INT(Nrot/2) perform full conformational search
    my @args = ("CF-$mol.round1.lh", "$DataDir");
    system($^X, "$CmdDir/lib/edts_squeeze.pl", @args) == 0 or die "system @args failed: $?";

    # *** next module now redundant using a sort instead ***
    #my @args = ("CF-$mol.round1.lh.sq", "$DataDir");
    #system($^X, "$CmdDir/edts_group.pl", @args) == 0 or die "system @args failed: $?";
    system("sort $DataDir/CF-$mol.round1.lh.sq > $DataDir/CF-$mol.round1.lh.sq.gp");

    my @args = ("CF-$mol.round1.lh.sq.gp", "$DataDir", "$mol");
    system($^X, "$CmdDir/lib/edts_combine.pl", @args) == 0 or die "system @args failed: $?";

    #system("sort $DataDir/CF-$mol.i | uniq | sed 's/\\.//g' | sed 's/\^/$mol./' > $DataDir/CF-$mol.round1.todo");
    system("sort $DataDir/CF-$mol.i | uniq | sed 's/\^/$mol./' > $DataDir/CF-$mol.round1.todo");
    #system("rm $DataDir/CF-$mol.i");

    my @args = ("CF-$mol.round1.todo", "$DataDir", "$mol");
    system($^X, "$CmdDir/lib/edts_recombine.pl", @args) == 0 or die "system @args failed: $?";

    system("mv $DataDir/CF-$mol.round1.todo.uniq $DataDir/CF-$mol.round2");

    # initialise round 3 file
    if (-e "$DataDir/CF-$mol.round3"){
        system("/bin/rm $DataDir/CF-$mol.round3");
    }
    system("touch $DataDir/CF-$mol.round2 $DataDir/CF-$mol.round3");

    # submits array of jobs based on input file list $mol.round1 and dependent job which 
    # performs round 2 processing
    printoutput($log, "Submitting round 2 jobs ...\n");
    open(roundlist, '<', "$DataDir/CF-$mol.round2") or die $!;
    while ($round=readline(roundlist)){
        printoutput($log,"$round");
    }
    close roundlist;

    my @args = ("$CmdDir/subarrayjob", "$DataDir/CF-$mol.round2", "$CmdDir", "$DataDir", "$CmdDir/EDTS_autorun_part3.pl $mol $CmdDir $DataDir 0 $atom1 $atom2 $gooddist $toldist");
    exec("/bin/bash", @args) == 0 or die "system @args failed: $?";
}
else {
    printoutput($log, "Moving to round 3...\n\n");

    # initialise round 3 file
    if (-e "$DataDir/CF-$mol.round3"){
        system("/bin/rm $DataDir/CF-$mol.round3");
    }
    system("touch $DataDir/CF-$mol.round2 $DataDir/CF-$mol.round3");

    my @args = ("$CmdDir/EDTS_autorun_part3.pl", "$mol", "$CmdDir", "$DataDir", "0", "$atom1", "$atom2", "$gooddist", "$toldist");
    #print("@args\n\n");
    exec("/bin/perl",@args) == 0 or die "system @args failed: $?";
}

# subroutine to echo STDOUT to mol.log file
# change print to printoutput above

sub printoutput {
    #my $text = @_;
    my ($q,$j) = @_;
    #print to STDOUT
    #print "$text\n";
    print "$j";
    
    # append text to logfile called mol.log
    # open/print/close approach ensures log text flush to disk/ viewable as script runs through
    
    #print "$log\n";
    open(mollog, '>>', "$q") or die $!;
    print mollog $j;
    close mollog;
}

__END__

