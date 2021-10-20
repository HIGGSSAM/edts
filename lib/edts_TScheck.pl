#! /usr/bin/perl 
$|++;

############################################################################
#       Check Output file is Valid (Good Distance) for TS Optimisation     # 
############################################################################    

#sub TScheck {
#    my $file = shift(@_);

$file =$ARGV[0];
$DataDir = $ARGV[1];
$CmdDir = $ARGV[2];
$LogFile = $ARGV[3];
$atom1=$ARGV[4];
$atom2=$ARGV[5];
$gooddist = $ARGV[6];
$toldist = $ARGV[7];


    open(logfile, '>>', "$LogFile") or die $!;

    $jobdone =-2 ;

    open(ifh,'<',"$DataDir/$file.out") or die $!;
    if ($atom1){ #doing ts
        system("$CmdDir/edts_zmattoxyz.pl xyz $DataDir/$file.out $LogFile");
        open(xyzfile,'<', "$DataDir/$file.xyz") or die $!;
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
        close(xyzfile);
        if ($atom[1] =~/Xx/){
            $bondlength =0.0;
        }
        else{
            $bondlength = 
            sqrt(($atomx[$atom1] - $atomx[$atom2])**2 + ($atomy[$atom1] - $atomy[$atom2])**2 + ($atomz[$atom1] - $atomz[$atom2])**2);
        }
    }
    while ($line = readline(ifh)) {
        if ($line=~ /Normal/){
            $jobdone=1;
        }
        elsif ($line=~ /Address/){
            $jobdone=-1;
        }
        elsif ($line=~ /kill/){
            $jobdone=-2;
        }
    }
    close (ifh);

    if (($jobdone <=-2) ) { 
        if ($atom1){
            if ( $bondlength<= $gooddist+$toldist and $bondlength >= $gooddist-$toldist){
                printf "%s killed, no resubmit\t %5.2f\t(within %5.2f-%5.2f )\n",$DataDir/$file,$bondlength,$gooddist-$toldist,$gooddist+$toldist;
                printf logfile "%s killed, no resubmit\t %5.2f\t(within %5.2f-%5.2f )\n",$DataDir/$file,$bondlength,$gooddist-$toldist,$gooddist+$toldist;
                $resub = "N";  # maybe yes
                #system("mailx -s 'Resubmit $file (killed), bond length is $bondlength within $gooddist-$toldist and $gooddist-$toldist' $email < $file.xyz ") ;
                print("Resubmit $file (killed), bond length is $bondlength within $gooddist-$toldist and $gooddist-$toldist");
                print logfile ("Resubmit $file (killed), bond length is $bondlength within $gooddist-$toldist and $gooddist-$toldist");
            }
            else{
                printf "%s killed, no resubmit\t %5.2f\t(shorter than %5.2f or longer than %5.2f)\n",$DataDir/$file,$bondlength,$gooddist-$toldist,$gooddist+$toldist;
                printf logfile "%s killed, no resubmit\t %5.2f\t(shorter than %5.2f or longer than %5.2f)\n",$DataDir/$file,$bondlength,$gooddist-$toldist,$gooddist+$toldist;
                $resub = "N";
                #system("mailx -s 'NO Resubmit $file (killed), bond length is $bondlength shorter than $gooddist-$toldist or longer than $gooddist+$toldist' $email < $file.xyz ") ;
                print("NO Resubmit $file (killed), bond length is $bondlength shorter than $gooddist-$toldist or longer than $gooddist+$toldist");
                print logfile ("NO Resubmit $file (killed), bond length is $bondlength shorter than $gooddist-$toldist or longer than $gooddist+$toldist");
            }
        }
        else{
            #printf "%s killed, resubmit the job\n",$DataDir/$file; 
            #$resub = "Y";
            printf "%s killed, No resubmit the job\n",$DataDir/$file; 
            printf logfile "%s killed, No resubmit the job\n",$DataDir/$file; 
            $resub = "N";

        }
        if ($resub eq "Y"){
            system("cat $DataDir/$file.out >> $DataDir/$file.out.2; /bin/rm $DataDir/$file.out"); 
            #return 0; 
        }
    }
    elsif (($jobdone ==-1) ) { 
        if ($atom1){
            if ( $bondlength<= $gooddist+$toldist and $bondlength >= $gooddist-$toldist){
                printf "%s Address Error, resubmit\t %5.2f\t(within %5.2f- %5.2f)\n",$DataDir/$file,$bondlength,$gooddist-$toldist,$gooddist+$toldist;
                printf logfile "%s Address Error, resubmit\t %5.2f\t(within %5.2f- %5.2f)\n",$DataDir/$file,$bondlength,$gooddist-$toldist,$gooddist+$toldist;
                $resub = "N"; # maybe yes
                #system("mailx -s 'Resubmit $file (Address Error), bond length is $bondlength (within $gooddist-$toldist and $gooddist +$toldist), but we adivice you to check what causes this Address Error Termination' $email < $file.xyz ") ;
                print("Resubmit $file (Address Error), bond length is $bondlength (within $gooddist-$toldist and $gooddist +$toldist), but we adivice you to check what causes this Address Error Termination");
                print logfile ("Resubmit $file (Address Error), bond length is $bondlength (within $gooddist-$toldist and $gooddist +$toldist), but we adivice you to check what causes this Address Error Termination");
            }
            else{
                printf "%s Address Error, no resubmit\t %5.2f\t(shorter than %5.2f or longer than %5.2f)\n",$DataDir/$file,$bondlength,$gooddist-$toldist,$gooddist+$toldist;
                printf logfile "%s Address Error, no resubmit\t %5.2f\t(shorter than %5.2f or longer than %5.2f)\n",$DataDir/$file,$bondlength,$gooddist-$toldist,$gooddist+$toldist;
                $resub = "N";
                #system("mailx -s 'NO Resubmit $file (Address Error), bond length is $bondlength (shorter than $gooddist -$toldist or longer than $gooddist+$toldist), ' $email < $file.xyz ") ;
                print("NO Resubmit $file (Address Error), bond length is $bondlength (shorter than $gooddist -$toldist or longer than $gooddist+$toldist)");
                print logfile ("NO Resubmit $file (Address Error), bond length is $bondlength (shorter than $gooddist -$toldist or longer than $gooddist+$toldist)");
            }
        }
        else{
            #printf "%s  Address Error, resubmit the job\n",$DataDir/$file; 
            #$resub = "Y";
            printf "%s  Address Error, No resubmit the job\n",$DataDir/$file; 
            printf logfile "%s  Address Error, No resubmit the job\n",$DataDir/$file; 
            $resub = "N";
        }
        if ($resub eq "Y"){
            system("cat $DataDir/$file.out >> $DataDir/$file.out.2; /bin/rm $DataDir/$file.out"); 
            #return 0; 
        }
    }    
    elsif ($jobdone ==1){
        if ($atom1){
            if ( $bondlength<= $gooddist+$toldist and $bondlength >= $gooddist-$toldist){
                #print color 'green';
                printf "%s terminated normally\t %5.2f\t(within %5.2f- %5.2f)\n",$DataDir/$file,$bondlength,$gooddist-$toldist,$gooddist+$toldist;
                printf logfile "%s terminated normally\t %5.2f\t(within %5.2f- %5.2f)\n",$DataDir/$file,$bondlength,$gooddist-$toldist,$gooddist+$toldist;
                #print color 'reset';
            }
            else{
                #print color 'green';
                printf "%s terminated normally.\t %5.2f\t(shorter than %5.2f or longer than %5.2f)\n",$DataDir/$file,$bondlength,$gooddist-$toldist,$gooddist+$toldist;
                printf logfile "%s terminated normally.\t %5.2f\t(shorter than %5.2f or longer than %5.2f)\n",$DataDir/$file,$bondlength,$gooddist-$toldist,$gooddist+$toldist;
                #print color 'reset';
            }
        }
        else{
            #print color 'green';
            printf "%s terminated normally.\n\n",$file;
            printf logfile "%s terminated normally.\n\n",$file;
            #print color 'reset';
        }
        #return 0;
    }

    close(logfile);

__END__
