#! /usr/bin/perl 
$|++;

############################################################################
#                                  OptEng                                  #
############################################################################    

#sub opteng{

# print("YES\n");

use Data::Dumper;

$roundlist = $ARGV[0];
#print("$roundlist\n");
$window = $ARGV[1];
$DataDir = $ARGV[2];
$CmdDir = $ARGV[3];
$atom1 = $ARGV[5];
$atom2 = $ARGV[6];
$gooddist = $ARGV[7];
$toldist = $ARGV[8];
$Nmax = $ARGV[4];

    $count = 0;
#   my ($roundlist,$window) = @_;
    open(roundlist, '<', "$DataDir/$roundlist") or die $!;
    while ($round = <roundlist>) {
        foreach $i ($round){
            chomp($i);
            @roundlist=split(/\./,$i);
            $dih[$count] = $roundlist[1];
            open(output, '<', "$DataDir/$i.out") or die $!;
            while (<output>) {
                if (/SCF Done/) {
                    @crap=split(/\s+/,$_);
                    #print Dumper @crap;
                    $eng[$count]=$crap[5];
                    #print Dumper $eng[$count];
                }
            }
            close output;
            
	    #print Dumper $eng[$count];
            
            if ($atom1){
                system($^X, "$CmdDir/edts_zmattoxyz.pl", "xyz", "$DataDir/$i.out");
                open(xyzfile, '<', "$DataDir/$i.xyz") or die $!;
                readline(xyzfile);
                readline(xyzfile);
                $natom=1;
                while ($atomline = readline(xyzfile)){
                    @atom = split(/\s+/,$atomline);
                    $atomx[$natom]= $atom[1];
                    $atomy[$natom]= $atom[2];
                    $atomz[$natom]= $atom[3];
                    $natom ++;
                }
                close xyzfile;
                $bonddist[$count] = 
                sqrt(($atomx[$atom1] - $atomx[$atom2])**2 + ($atomy[$atom1] - $atomy[$atom2])**2 + ($atomz[$atom1] - $atomz[$atom2])**2);

                if ( $bonddist[$count] >= $gooddist+$toldist or $bonddist[$count] <= $gooddist-$toldist){
                    print("AHHHH\n");
                    $eng[$count] = 0.0;
                }
                #print "$dih[$count] , $bonddist[$count]\n";
            }
            #print ("$dih[$count], $eng[$count]\n");
            #print Dumper $eng[$count];
            $count ++;
        }
    }
    #print("count is $count\n");
    $count1 = $count/2;
    #print("$count1\n");
    #print Dumper $eng[$count];
    
    close roundlist;

    sub numerically { $a <=> $b; }
    @sorteng = sort numerically @eng;


    $sortk=0;
    while ($sortk <$count){
        for ($k = 0; $k < $count ; $k++){
            if ($sorteng[$sortk] == $eng[$k] ){
                $sortdih[$sortk] = $dih[$k];
                $sortk++;
            }
        }
    }


    open(fla, '>', "$DataDir/$roundlist.all") or die $!;
    open(flh, '>', "$DataDir/$roundlist.lh") or die $!;
    open(flw, '>', "$DataDir/$roundlist.window") or die $!;
    open(fuh, '>', "$DataDir/$roundlist.uh") or die $!;

    $Lowest = ($sorteng[1]-$sorteng[0]) *2625.5;
    if ($Lowest < $window){
        for ($sortk = 0; $sortk < $count ; $sortk++){  
            $Deng[$sortk]= ($sorteng[$sortk]-$sorteng[0]) *2625.5;  
            if (($Deng[$sortk] < $window && $Deng[$sortk]> -$window)&& $sortk < ($Nmax - 1) ){
                printf flw "Less than %.2fkJ/mol= %s \t%f\n", $window, $sortdih[$sortk], $Deng[$sortk];
            }
            if ($sortk < ($count - 1)/2) {
                printf "Lower Half = %s \t%f\n", $sortdih[$sortk], $Deng[$sortk];
		printf flh "Lower Half = %s \t%f\n", $sortdih[$sortk], $Deng[$sortk];
            }
            else {
                printf "Upper Half = %s \t%f\n", $sortdih[$sortk], $Deng[$sortk];
                printf fuh "Upper Half = %s \t%f\n", $sortdih[$sortk], $Deng[$sortk];
            }   
            if ($sortk < ($count - 1)/2) {
                #printf "Lower Half = %s \t%f\n", $sortdih[$sortk], $Deng[$sortk];
		#printf flh "Lower Half = %s \t%f\n", $sortdih[$sortk], $Deng[$sortk];
                printf fla "Lower Half = %s \t%f kJ/mol or %f Hartree\n", $sortdih[$sortk], $Deng[$sortk], $sorteng[$sortk];
            }
            else {
                #printf "Upper Half = %s \t%f\n", $sortdih[$sortk], $Deng[$sortk];
                #printf fuh "Upper Half = %s \t%f\n", $sortdih[$sortk], $Deng[$sortk];
                printf fla "Upper Half = %s \t%f kJ/mol or %f Hartree\n", $sortdih[$sortk], $Deng[$sortk], $sorteng[$sortk];
            } 
        }

    print "\n-----------------------------------------------------------\n";
    printf "Lowest energy conformer is %s \t%f Hartree\n", $sortdih[0], $sorteng[0];
    printf "The next lowest is %s \t+%f kJ/mol\n", $sortdih[1] ,$Lowest ;
    print "-----------------------------------------------------------\n";
    }
    else {
        print "\n-----------------------------------------------------------\n";
        printf "Lowest 1 = %s \t%f Hartree\n", $sortdih[0], $sorteng[0];
        printf "No other conformer has energy less than %.2fkJ/mol\n", $window;
        print "-----------------------------------------------------------\n";
        printf flh "Lowest 1 = %s \t%f\n", $sortdih[0], $Deng[0];
        printf fla "Lowest 1   = %s \t%f kJ/mol or %f Hartree\n", $sortdih[0], $Deng[0], $sorteng[0];
        printf flw "Less than %.2fkJ/mol= %s \t%f\n", $window, $sortdih[0], $Deng[0];
        for ($sortk = 1; $sortk < $count ; $sortk++){  
            $Deng[$sortk]= ($sorteng[$sortk]-$sorteng[0]) *2625.5;
            printf fuh "Upper Rest = %s \t%f\n", $sortdih[$sortk], $Deng[$sortk];
            printf fla "Upper Rest = %s \t%f kJ/mol or %f Hartree\n", $sortdih[$sortk], $Deng[$sortk], $sorteng[$sortk];
        }
    }

    close(fla);
    close(flh);
    close(flw);
    close(fuh);




__END__
