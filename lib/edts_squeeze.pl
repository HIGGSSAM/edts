#! /usr/bin/perl 
$|++;



############################################################################
#                                  Squeeze                                 #
############################################################################    

#sub squeeze{
#    my $tosq=shift(@_);

$tosq = $ARGV[0];
$DataDir = $ARGV[1];


open(filesq, '>', "$DataDir/$tosq.sq") or die $!;
open(tosq, '<', "$DataDir/$tosq") or die $!;
    
$count=0;
while ($dihhh= readline(tosq)){
    @crap=split(/\s/,$dihhh);
    $dih = $crap[3];
    if (length($dih) <4 ){
        print filesq "$dih\n";
    }
    else {
        @IndDih= split(/a1/,$dih);
        print filesq "$IndDih[1]\n";
    }
    }
close tosq;
close filesq;


__END__