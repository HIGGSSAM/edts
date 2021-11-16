#! /usr/bin/perl 
$|++;


############################################################################
#                                  ReCombine                               #
############################################################################    

#sub recomb{
 #   my $torc = shift(@_);

$torc = $ARGV[0];
$DataDir = $ARGV[1];
$mol = $ARGV[2];

    chomp($torc);
    open(torc, '<', "$DataDir/$torc") or die $!;
    open(moldone, '<', "$DataDir/CF-$mol.done") or die $!;
    $count1=0;
    $count2=0;
    while ($crap1[$count1] = readline(torc)){
        chomp($crap1[$count1]);
        @crap2 = split(/\./,$crap1[$count1]);
        $todo[$count1] = $crap2[1];
        $count1++;
    }
    while ($crap2[$count2] = readline(moldone)){
        chomp($crap2[$count2]);
        @crap3 = split(/\./,$crap2[$count2]);
        $done[$count2] = $crap3[1];
        $count2++;
    }
    open(todo, '>', "$DataDir/$torc.uniq") or die $!; 
    for ($kk=0;$kk< $count1 ; $kk++){
        for ($jj=0;$jj<$count2;$jj++){
            if ("$done[$jj]" eq "$todo[$kk]"){
                $iii=1;
                last;   #break out of loop when match                
            }
            else{
                $iii=0;
            }
        }
        if ($iii ==0){
            print todo "$mol.$todo[$kk]\n";
        }
    }
    close moldone;
    close todo;
    close torc;

__END__