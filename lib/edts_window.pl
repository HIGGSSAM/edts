#! /usr/bin/perl 
$|++;


############################################################################
#                          Window                                          #
#     < CF-$mol.done.window CF-$mol.round1.uh.sq.$i > CF-$mol.window       #
############################################################################    

#sub window {
#    my ($donewd,$uhi)=@_;

$donewd =$ARGV[0];
$uhi = $ARGV[1];
$DataDir = $ARGV[2];
$mol = $ARGV[3];

    open (wd,">", "$DataDir/CF-$mol.window") or die $!;

    open (uhi,"<", "$DataDir/$uhi") or die $!;
    while ($dihh= readline(uhi)){
        @crap=split(/\s/,$dihh);    
        $dih = $crap[0];
        print wd "$dih ";
    }
    print wd "\n";
    close uhi;

    open (donewd, "<", "$DataDir/$donewd") or die $!;
    $count=0;
    while ($dihh= readline(donewd)){
        @crap=split(/\s/,$dihh);
        $dih = $crap[3];
        # *** change to print out dih as separate records and no full stops
        #$dih =~ s/(\d)/$1\./g;
        #print wd "$dih ";
        print wd "$dih\n";
    }

    close donewd;
    close wd;



__END__