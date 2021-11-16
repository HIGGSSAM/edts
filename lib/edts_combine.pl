#! /usr/bin/perl
$|++;
use Data::Dumper;


############################################################################
#                                  Combine without.                        #
############################################################################    

#sub combine{

@data = ();

#    my $tocomb = shift(@_);
$tocomb =$ARGV[0];
#print("$tocomb\n");
$DataDir = $ARGV[1];
#print("$DataDir\n");
$mol = $ARGV[2];
#print("$mol\n");

open('tocomb', '<', "$DataDir/$tocomb") or die $!;
while ($nnn=readline('tocomb')){
    push(@data,(split(/\s+/,$nnn)));;
}
#print "data\n";
#print Dumper @data;

close 'tocomb';

open(moli, '>', "$DataDir/CF-$mol.i") or die $!;
# start at 1 to skip the empty selection
for ($i = 1; $i < 2**@data; $i++) {
    #print("i is $i\n");
    @bin = dec2bin($i);
    #print "bin\n";
    #print Dumper @bin;
    my @tmp = ();
    for ($j = 0; $j < @bin ;$j++) {
        push(@tmp,$data[$j]) if ($bin[$j]);;
    }
    #print "STARTING perm with tmp\n";
    #print "...\n";
    #print "tmp\n";
    #print Dumper @tmp;
    perm("",@tmp);
}
close moli;

    
sub dec2bin {
    my $str = unpack("B32", pack("N", shift));
    $str =~ s/^0+(?=\d)//;
    return split(//,reverse($str));
}

sub perm {
    my ($s,@a) = @_;
    my $i;
    #print "STARTING perm\n";
    #print "@ a is:\n";
    #print Dumper @_; 
    #print Dumper @a;
    #print " s is:\n";
    #print Dumper $s;
    #print " i is:\n";
    #print $i;
    
    
    if (@a) {
        #print("something is in a! a is printed.\n");
        #print Dumper @a;
        #print "$a[0]\n";
        #print "@{$a[0]}\n";
        
	#my @b = @{$a}[0];
	#print("b is printed.\n");
        #print Dumper @b;
        
	my @b = shift (@a);
        #print $_;
        #print("shifted a is printed.\n");
        #print Dumper @a;
	
	for (@b) {
	    #print("something is in b! b is printed.\n");
            #print Dumper @b;
            #change this line later
            #print Dumper "$_.";
            
            # *** don't need full stops ***
	    #perm($s."$_.",@a);
	    perm($s."$_",@a);
        }
    }else {
            
        #print length($s),"\n";
        #print "$s\n";
        
        # multiple rotations
        # *** test less than 2 now full stops removed
        #if (length($s)>3){
        if (length($s)>2){
            #print("length greater than 2\n");
            if ($s =~ /a/) {
                # insert test for multi rotations on same bond
                my $rot = $s;
                $rot =~ s/\d//g;
                #$rot =~ s/\.//g;
                my $str = $rot;
                #print("$str\n");
                #$rot =~ tr/a-z/a-z/s;
                $rot =~ s/(.)(?=.*?\1)//g;
                #print("$rot\n"); 
                if ($str eq $rot) { 
                    # $s doesnt contain multiple letters then print 
                    print moli "$s\n";
                }
	    }else{
	        # insert test for multi rotations on same bond
	        my $rot = $s;
                $rot =~ s/\d//g;
                #$rot =~ s/\.//g;
                my $str = $rot;
                #print("$str\n");
                #$rot =~ tr/a-z/a-z/s;
                $rot =~ s/(.)(?=.*?\1)//g;
                # $s doesnt contain multiple letters then print
                if ($str eq $rot) {
                    #print moli "a1.$s\n";
                    print moli "a1$s\n";
                }
            }
        }# single rotation
        else{
            #print("length => 2\n");
            # if $s doesnt contain 'a' then print
            if (!$s eq /[a]/) {
                #print "$s add a1 to front.\n";
                #print moli "a1.$s\n";
                print moli "a1$s\n";
                }
	    }
        
        }
        return;
    }

    

__END__