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
# EDTS_autorun.              PART 3/3
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
# NOTE [] are for transition state optimisation, if r12 is fixed, this can be used to help determine the optimised geometry is valid.                    #
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
$mol=$ARGV[0];  
$CmdDir=$ARGV[1];
$DataDir=$ARGV[2];
$uhix=$ARGV[3];

# optional arguments
$atom1=$ARGV[4];
$atom2=$ARGV[5];
$gooddist=$ARGV[5];
$toldist=$ARGV[6];
 
chomp($mol);

$log = "$DataDir/$mol.log";
#printoutput($log, "$mol\n");
#printoutput($log, "$CmdDir\n");
#printoutput($log, "$DataDir\n");
#printoutput($log, "$uhix\n");

$curround=$uhix + 1;
printoutput($log, "\n###################################################################################\n"); 
printoutput($log, "#                                        EDTS                                     #\n");
printoutput($log, "###################################################################################\n\n"); 

printoutput($log, "Conformer search for $mol.\n");
printoutput($log, "Data directory: $DataDir\n\n");
printoutput($log, "###################################################################################\n");
printoutput($log, "#                                  Round3 - stage $curround                               #\n");
printoutput($log, "###################################################################################\n");

# initialise round 3 file
if (-e "$DataDir/CF-$mol.TScheck.results"){
   system("/bin/rm $DataDir/CF-$mol.TScheck.results");
}

if ($uhix ne "0") {
    $uhixx=$uhix - 1;  #look at previously run jobs
    open(R3tosub, "<", "$DataDir/CF-$mol.round3-$uhixx.tosub") or die $!;
    while ($R3tosub = <R3tosub>){
    chomp($R3tosub);
        foreach $y ($R3tosub){
            my @args = ("$y", "$mol", "$DataDir", "$CmdDir", "$log", "$atom1", "$atom2", "$goodist", "$toldist");
            system($^X, "$CmdDir/lib/edts_TScheck.pl", @args) == 0 or die "system @args failed: $?";
        }
    }
    close R3tosub;
}
else{
    open (round2, "<", "$DataDir/CF-$mol.round2")or die $!;
    while ($round2 = <round2>){
    chomp($round2);
        foreach $y ($round2){
            my @args = ("$y", "$mol", "$DataDir", "$CmdDir", "$log", "$atom1", "$atom2", "$goodist", "$toldist");
            system($^X, "$CmdDir/lib/edts_TScheck.pl", @args) == 0 or die "system @args failed: $?";
        }
    }
    close round2;
}

$errcount=0;
$lnecount=0;
# count how many TScheck errors detected
if (-e "$DataDir/CF-$mol.TScheck.results"){
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
}

open (r1uh, "<", "$DataDir/CF-$mol.round1.uh.sq") or die $!;

$uhi = 0;
while ($r1uh = <r1uh>){
    #print ("r1uh=$r1uh\n");
    if ($uhix eq $uhi) {
        foreach $x ($r1uh){
            open(r1uhlist, ">", "$DataDir/CF-$mol.round1.uh.sq.$uhi") or die $!;    
            print r1uhlist $x;
            close r1uhlist;
        }
        system("cat $DataDir/CF-$mol.round1 $DataDir/CF-$mol.round2 $DataDir/CF-$mol.round3 | sort | uniq > $DataDir/CF-$mol.done");
       
        my @args = ("CF-$mol.done", "$EC2","$DataDir", "$CmdDir", "$log", "$Nmax", "$atom1", "$atom2", "$goodist", "toldist");
        system($^X, "$CmdDir/lib/edts_opteng.pl", @args) == 0 or die "system @args failed: $?";

        my @args = ("CF-$mol.done.window", "CF-$mol.round1.uh.sq.$uhix", "$DataDir", "$mol");
        system($^X, "$CmdDir/lib/edts_window.pl", @args) == 0 or die "system @args failed: $?";

        # ensure dihs in mol.window sorted alphabetically
        system("sort $DataDir/CF-$mol.window | uniq > $DataDir/CF-$mol.window.sq");

        my @args = ("CF-$mol.window.sq", "$DataDir", "$mol");
        system($^X, "$CmdDir/lib/edts_combine.pl", @args) == 0 or die "system @args failed: $?";
        
        system("sort $DataDir/CF-$mol.i | sed 's/\^/CF-$mol./' > $DataDir/CF-$mol.order");
        # system("rm $DataDir/CF-$mol.i");

        # *** next module redundant use sort instead ***
        # my @args = ("CF-$mol.needreorder", "$DataDir", "$mol");
        # system($^X, "$CmdDir/edts_reorder.pl", @args) == 0 or die "system @args failed: $?";
        
        # Before submit round3-1, check if there is anything in DONE!
        my @args = ("CF-$mol.order", "$DataDir", "$mol");
        system($^X, "$CmdDir/lib/edts_recombine.pl",  @args) == 0 or die "system @args failed: $?";

        system("sort $DataDir/CF-$mol.order.uniq | uniq > $DataDir/CF-$mol.round3-$uhi.tosub");
        system("cat $DataDir/CF-$mol.round3-$uhi.tosub >> $DataDir/CF-$mol.round3");     

        $uhix++; 
        #printoutput($log, "Beginning round 3-$uhix\n");

        # insert fail save if .tosub file is empty.
        $curround=$uhix + 1;
        if (-s "$DataDir/CF-$mol.round3-$uhi.tosub"){
            # submits array of jobs based on input file list $mol.round1 and dependent job which 
            # performs round 2 processing
            
            printoutput($log, "Submitting round 3-$curround jobs ...\n");
            open(roundlist, '<', "$DataDir/CF-$mol.round3-$uhi.tosub") or die $!;
            while ($round=readline(roundlist)){
                printoutput($log,"$round");
            }
            close roundlist;

            my @args = ("$CmdDir/subarrayjob", "$DataDir/CF-$mol.round3-$uhi.tosub", "$CmdDir", "$DataDir", "$CmdDir/EDTS_autorun_part3.pl $mol $CmdDir $DataDir $uhix $atom1 $atom2 $gooddist $toldist");
            exec("/bin/bash", @args) == 0 or die "system @args failed: $?";
        }
        else{
            printoutput($log, "No new jobs to submit during round 3-$curround ...\n");
        }
    }
    $uhi++;
}
close r1uh; 
system("cat $DataDir/CF-$mol.round1 $DataDir/CF-$mol.round2 $DataDir/CF-$mol.round3 | sort | uniq > $DataDir/CF-$mol.done");

# pick lowest conformation 
my @args = ("CF-$mol.done", "$EC2", "$DataDir", "$CmdDir", "$log", "Nmax", "$atom1", "$atom2", "$goodist", "toldist");
system($^X, "$CmdDir/lib/edts_opteng.pl", @args) == 0 or die "system @args failed: $?";

# remove temp files  NOTE do "*.gp from group.pl" and "*.i  from combine.pl" files also need clearing?
system("rm $DataDir/CF-$mol.*order* $DataDir/CF-$mol.*sq* $DataDir/CF-$mol.*.ll");


printoutput($log, "\n************************************************\n"); 
printoutput($log, "     Conformer search completed for $mol\n");
printoutput($log, "\************************************************\n\n"); 

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
