#! /usr/bin/perl 
use Text::Wrap;

################################################################
#                         edts_zmattoxyz 1.0
#			            
#  Based on CRIN.PL by David Brittain
#                         11/10/06
#
#  Modified Leaf 04Jan08
#
# Program to take file containing geometries, extract
# the geometries, and create Gaussian Input
#
# Syntax: crin [output format] $mol.ext $diffuse $temp 
#		
#       $mol.[out]
#	$diffuse blank or 0 = no diffuse function, 1= add diffuse function
#	$temp = temperature. Blank defaults to 298.15K
#
# Which is which?
# .out : Gaussian output
#
# 
#
# Features:
#	Converts Gaussian zmat output file format into xyz format file
#
#	Works with both zmat and cartesian coordinates in all input
#	formats. Even supports old g98 archive zmat format.
#
#	Gaussian:
#	Supports multiple --Link1-- jobs. So long as the header and
#	footer are matched in that the %chk lines are the same, the 	
#	output will run correctly. The charge/multiplicity line is
#	corrected from whatever it is in the footer file to whatever
#	it is in the input geometry file.
#
#	All comment lines are changed to the current (stripped)
#	filename, even if there are multiple --Link1-- jobs.
#
#
#	Puts (stripped) filename in the title card.
#
#
#Changelog:
#
#
#
################################################################

#declarations:


$grid="INT(grid=ultrafine)"; 
#$grid="";

if(@ARGV <2){
    die "More arguments needed.\nSyntax: crin [output format] $mol.[ext] $diffuse $temp  \n";
}

$job = $ARGV[0];
if($job =~ /xyz/i  ) {$output_type="gaussian";}
else         {die "unsupported output type $job.\n";}

# ARGV[2] is the diffuse fun switch
if (!$ARGV[2] or $ARGV[2] =~ /0/){
    #print "no diffuse function\n";
    $basis = "6-31G*";
    $genbasis = "6-31Gd";
    $heavybasis = "lacvp + polarisation + ecp";
}
else{
    print "add diffuse function\n";
    $basis = "6-31+G*";
    $genbasis = "6-31pGd";
    $heavybasis = "lacvp + polarisation + diffuse + ecp";
}
if (!$ARGV[3]){
    $temperature = 298.15;
}
else{
    $temperature = $ARGV[3];
}

#check to see all of the files exist

#strip the geometries from whatever the input file is.
$inputfile = $ARGV[1];
if ($inputfile =~  /\.out/ ){
    if($inputfile =~  /\.298\.mecn\.out/ and $job =~ /mecn/i){
        $inputfile_stripped = stripname($inputfile,"\.298\.mecn\.out");
    }
    else{
        $inputfile_stripped = stripname($inputfile,"\.out");	# strip .out from $inputfilename
    }
    $geometry = striparc($inputfile,0);				# get geometry string
    if ($geometry eq ""){
        print "Caution! $inputfile This is not optimised geometry!\n";
        $stripped_geometry = striplastxyz($inputfile);
        print("$stripped_geometry");
    }
    else{
        $stripped_geometry = stripcom($geometry);
	print("$stripped_geometry");			# cull leading and trailing rubbish
    }
    @atom = split(/\n/,$stripped_geometry);
    #@atomline = split(/\n/,$stripped_geometry);
    for $i (1..$#atom){
        $atom[$i] =~ s/-*\d+.*//g;
        $atom[$i] =~ s/\s+//g;
        $atomn = atomicN($atom[$i]);
        if ($atomn > 18){
            $Rassolov=1;
            if ($atomn > 36){
                print "use $heavybasis for $atom[$i]\n";
                $heavy=1;
            }
        }
    }
}

else {
    print "Unknown file extension. Supported filetypes: .out (Gaussian)\n";
    die;
}


#Now, create the output.


if ($output_type eq "gaussian"){
    $cm = getcm("$stripped_geometry");
    $iszmat = check_geometry_type($stripped_geometry);
    #if it is, only print the first 7 fields of the zmat
    #if the no. of fields is 8.
    ##$finaldeck = ${stripped_geometry}."\n";
    #$finaldeck = addcomment($inputfile_stripped, $finaldeck);
    if($job =~  /xyz/i ){
        if ($inputfile =~ /\.out$/){
            @xyz = split(/\n/,$stripped_geometry);
            $header ="$#xyz\n";
            $stripped_geometry =~ s/,/ /g;
        }
        else{
            if ($iszmat){
                $header = "";
            }
            else{
                @xyz = split(/\n/,$stripped_geometry);
                $header ="$#xyz\n";
            }
        }
        $footer ="";
    }

    $footer = ${footer}."\n";
    $finaldeck = ${header}.${stripped_geometry}."\n".${footer};
    if($job =~  /xyz/i)   {writedeck($finaldeck,$inputfile_stripped.".xyz");}
}



sub stripcom{
#stripcom("geometry_string"), returns stripped geometry string.

#declarations:
    my($geometry, @file, $cm_flag, $bk, @stripped_geometry, $cm_index, $i);
    my($temp2, @temp, $iszmat, $stripped_geometry);

    $geometry = $_[0];
    @file = split(/\n/,$geometry);
    $cm_flag = 0;   # whether we've hit the (first) charge/multiplicity section or not.
    $bk = 0;        #blank line counter for zmat section
    @stripped_geometry = ();
    $cm_index = 9999999;

    for($i = 0; $i <= $#file; $i++){


        #match the charge/multiplicity section.
        if ($file[$i] =~/^[ \t]*-*\+*[ \t]*[0-9]{1,2}[ \t]*\+*[ \t]*[0-9]{1,2}[ \t]*$/){
            $cm_flag +=1;			#means you only get the top geometry.
            if ($cm_flag == 1){
                push(@stripped_geometry, $file[$i]);		#put cm into array
                $cm_index = $i;				#mark where cm is
                $temp2 = $file[$i+1];				#get next line
                $temp2 =~ s/^[ \t]*//;				#strip leading whitespace
                $temp2 =~ s/[ \t]*$//;				#strip trailing whitespace

                @temp = split(/[ \t]+/,$temp2);		#split


                if ($temp[0] =~ /\b[a-z]{1,2}\b/i){		#see if zmat
                    if (@temp == 1){
                        $iszmat = 1;
                    }else{
                        if (@temp == 4 or @temp == 5){
                            $iszmat = 0;
                        }else{
                            die "$temp2 : Not xyz or zmat!\n";
                        }
                    }
                }
                else{
                    die "$temp[0] : not an atom!!\n";
                }
            }
        }
        $file[$i] =~s/ 0 //;
        if ($i > $cm_index && $cm_flag == 1){	#discard useless header and footer
            if($iszmat == 0 && $file[$i] =~/^[ \t]*[a-zA-Z]{1,2}[ \t]{1,}-?\+? *\.?[0-9]{1,}\.?[0-9]*[ \t]{1,}-?\+? *\.?[0-9]{1,}\.?[0-9]*[ \t]{1,}-?\+? *\.?[0-9]{1,}\.?[0-9]*[ \t]*$/){		#if zmat and a coordinate
                push(@stripped_geometry, $file[$i]);
            }
            if($iszmat == 1){
                if($file[$i] =~ /^[ \t]*$/){
                    $bk+=1;				#count blank lines in zmat
                }
                if ($bk < 2){		# if zmat and before the end of it, then
                    push(@stripped_geometry, $file[$i]);
                }
            }
        }

    }





    $stripped_geometry = join("\n",@stripped_geometry);
    return "$stripped_geometry\n";
}

sub striplastxyz {
    my @stripped_geometry=("");
    @stripped_geometry = ();
    my $charge="";
    my $multi="";
    @symbo=("H" ,                                                                                "He",
        "Li","Be",                                                   "B", "C", "N", "O", "F","Ne",
        "Na","Mg",                                                  "Al","Si", "P", "S","Cl","Ar",
        "K" ,"Ca","Sc","Ti", "V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
        "Rb","Sr", "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te", "I","Xe",
        "Cs","Ba",
        "La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
        "Hf","Ta", "W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
        "Fr","Ra",
        "Ac","Th","Pa", "U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr");

    open(LOG, $_[0]) || die "Can't open $_[0] for reading.\n";
    while (<LOG>){
        if (/Multiplicity/){
            @mc = split(/\s+/,$_);
            $charge = $mc[3];
            $multi = $mc[6];
        }
        if (/orientation/){
            @last_geometry="";
            $line = readline(LOG);
            $line = readline(LOG);
            $line = readline(LOG);
            $line = readline(LOG);
            $stop = 0;
            while ($stop < 1){
                $line = readline(LOG);
                if ($line =~ /-------/){
                    $stop=1;
                }
                else{
                    $stop=0;
                    push @last_geometry, $line;
                }
            }
        }
    }
    close LOG;
    push @stripped_geometry, "$charge $multi";
    for ($i = 1 ;$i <= $#last_geometry ; $i++){
        @crap= split(/\s+/,$last_geometry[$i]);
        $atom[$i] = $crap[2];
        $x[$i] = $crap[4];
        $y[$i] = $crap[5];
        $z[$i] = $crap[6];
        push @stripped_geometry, "$symbo[$atom[$i]-1]\t$x[$i]\t$y[$i]\t$z[$i]";
    }
    $stripped_geometry = join("\n",@stripped_geometry);
    return "$stripped_geometry\n";
}

sub atomicN {

    my ($a1) = @_;

    @symbo=("h" ,                                                                                "he",
        "li","be",                                                   "b", "c", "n", "o", "f","ne",
        "na","mg",                                                  "al","si", "p", "s","cl","ar",
        "k" ,"ca","sc","ti", "v","cr","mn","fe","co","ni","cu","zn","ga","ge","as","se","br","kr",
        "rb","sr", "y","zr","nb","mo","tc","ru","rh","pd","ag","cd","in","sn","sb","te", "i","xe",
        "cs","ba",
        "la","ce","pr","nd","pm","sm","eu","gd","tb","dy","ho","er","tm","yb","lu",
        "hf","ta", "w","re","os","ir","pt","au","hg","tl","pb","bi","po","at","rn",
        "fr","ra",
        "ac","th","pa", "u","np","pu","am","cm","bk","cf","es","fm","md","no","lr");


    my %index ;
    @index{@symbo} = (0..$#symbo);
    my $index = $index{lc($a1)} + 1 ;
    return $index; 
}


sub striparc {
#striparc("inputfile",1|0) : returns geometry from .log, with 1 gives archive form, and 0 breaks into xyz lines

    my ( @file);
    my ( @archive);
    my ($i);
    local $/;


    $arc = $_[1];

    open(LOG, $_[0]) || die "Can't open $_[0] for reading.\n";
    undef $/;
    @file = split(/\n\n+/,<LOG>,-1);
    close(LOG)||die "Can't close $_[0].\n";
    for ($i= 0;$i <= $#file;$i++){
        if ($file[$i]=~/GINC(.+?)@/s){
            push(@archive, $1);
        }
    }
    if ($arc){
            $archive[0] =~ s/\s\n//g;
            $archive[0] =~ s/\n+\s//g;
            $archive[0] =~ s/DipoleDeriv.*PG/PG/g;
            $archive[0] =~ s/NImag.*//s;

        return $archive[0];
    }
    else{
        if ($archive[0]){
            $archive[0] =~ s/ //g;
            $archive[0] =~ s/\n//g;
            $archive[0] =~ s/\r//g;
            $archive[0] =~ s/\\/\n/g;
            $archive[0] =~ s/,/ /g;
            $archive[0] =~ s/ /  /g;
            $archive[0] =~ s/  -/ -/g;
            return $archive[0];
        }
        else{
            return "";
        }
    }
}


sub stripname {
    my ($strippedfile);
    $strippedfile = $_[0];
    $strippedfile =~ s/$_[1]\b//;
    return $strippedfile;
}



#sub strip {
# strip("filename") : gets a txt file, turns it into a string.
#    my ($i,$file,$temp);
#    local $/;

#    open(INPUTFILE, $_[0]) ||die "Can't open $_[0] for reading.\n";
#    undef $/;
#    $file = <INPUTFILE>;
#    close(INPUTFILE)|| die "Can't close $_[0].\n";
#    return $file;
#}

sub writedeck{
#writedeck("string", "filename") : writes "string" to "filename"
    my ($temp);

    open(OUTFILE, ">$_[1]") ||die "Can't open $_[0] for writing.\n";
    $temp = $\;
    print OUTFILE $_[0];

    close(OUTFILE) || die "Can't close $_[1].\n";
}

sub getcm {
#getcm("geometry") : gets charge/multiplicity line from "geometry" string.
    my ($geometry, @geometry, $cm);

    $geometry = $_[0];
    @geometry = split(/\n/, $geometry);
    $cm = $geometry[0];
    return $cm;
}



sub addcomment{
#addcomment("inputfile", "comment") : substitutes the comment
#line(s) in a (complete) gaussian .com string with "comment"

    my($comment, @inputfile, $inputfile, $i);

    ($comment, $inputfile) = ($_[0],$_[1]);
    @inputfile = split(/\n/, $inputfile,-1);
    for($i = 0; $i <= $#inputfile; $i++){
        if ($inputfile[$i] =~ /^[ \t]*-*\+*[ \t]*[0-9]{1,2}[ \t]*\+*[ \t]*[0-9]{1,2}[ \t]*$/){
            $inputfile[$i-2] = $comment;
        }
    }
    $inputfile = join("\n",@inputfile);
    return $inputfile;


}

sub check_geometry_type{
# check_geometry_type("stripped_geometry") : checks to see if a stripped
#geometry is xyz or zmat.

    my($stripped_geometry,@stripped_geometry, $i, $temp2, $temp, @temp, $iszmat);

    $stripped_geometry = $_[0];

    @stripped_geometry = split(/\n/,$stripped_geometry);


    for($i = 0; $i <= $#stripped_geometry; $i++){
        if ($stripped_geometry[$i] =~/^[ \t]*-*\+*[ \t]*[0-9]{1,2}[ \t]*\+*[ \t]*[0-9]{1,2}[ \t]*$/){
            $temp2 = $stripped_geometry[$i+1];                          #get next line
            $temp2 =~ s/^[ \t]*//;                         #strip leading whitespace
            $temp2 =~s/[ \t]*$//;                          #strip trailing whitespace
            @temp = split(/[ \t]+/,$temp2);                #split

            if ($temp[0] =~ /\b[a-z]{1,2}\b/i){            #see if zmat
                if (@temp == 1 ){
                    $iszmat = 1;
                }else{
                    if (@temp == 4){
                        $iszmat = 0;
                    }else{
                        die "$temp2 : Not xyz or zmat!\n";
                    }
                }
            }
            else{
                die "$temp[0] : not an atom!!\n";
            }
        }
    }
    return $iszmat;
}



__END__