#!usr/bin/perl -w

#Define filename
my $targetfile = $ARGV[0];
my $Randomizations = $ARGV[1];

open (INFILE, "$targetfile") || die "Can't open the infile!\n";
my @Fasta = <INFILE>;
close(INFILE) || die "can't close file";

my $N = @Fasta;


#Put FASTA seqs into a Hash

my %FastaHash = ();
my $CurrentSeq = "";

#Create array of names to keep track of order
my @SeqNames = ();

for ($i = 0; $i < $N; $i++) {

    $Line = $Fasta[$i];
    chomp $Line;

    #Place title and seq into their respective places in the FastaHash
	
	#Test if title
    if (substr($Line, 0, 1) eq ">") {

        my $SeqName = $Line;
        chomp $SeqName;
		$SeqName =~ s/\R//g;
        $SeqName =~ s/>//g;
        $SeqName =~ s/:/-/g;
        $CurrentSeq = $SeqName;
        push (@SeqNames, $CurrentSeq);
        }
    
	#Test if sequence data (must be IUPAC code)
    if ( $Line =~ m/(A|G|C|U|T|Y|R|K|M|B|D|H|V|N|S|W)/g && substr($Line, 0, 1) ne ">") {

    $FastaHash{$CurrentSeq} .= $Line;

    }
}

foreach my $SequenceName (@SeqNames) {

    chomp $SequenceName;    
    my $Sequence = $FastaHash{$SequenceName};
    chomp $Sequence;
    $Sequence =~ s/\s+//g;
    $Sequence =~ s/-//g;
    $Sequence = uc $Sequence;

	    #my $Command = "echo " . $Sequence . " | RNAfold";
	    #my @Out = `$Command`;


        #Put the native sequence in the first position always!
        my @seqarray = ();
        push (@seqarray, $Sequence);
        my @ScrambledSeqs = Scramble ($Sequence);
        push (@seqarray, @ScrambledSeqs);

    #my $NTFreqs = &NucFreqs($SEQ);     
	my @EnergyArray = &Energy(@seqarray);
	my $Zscore = &ZScore(@EnergyArray);
    my $NativeDG = $EnergyArray[0];
	chomp $NativeDG;
	my $PValue = &PValue(@EnergyArray);
	
   	print "$SequenceName\t$NativeDG\t$Zscore\t$PValue\n";
    #print "\n\n";    
	#foreach my $Element (@EnergyArray) { print "$Element";}    

}


`rm *.seq`;
`rm  rna.ps`;

######Sub-routine to scramble RNAs################################
sub Scramble {

        my $InSeq = $_[0];
        my @Out = ();

        for (my $i = 0; $i < $Randomizations; $i++) {
                my $OutSeq = "";
                my @InSeq = split ("", $InSeq);

                while (@InSeq > 0) {
                        my $Rand = rand(@InSeq);
                        my $RandBase = splice(@InSeq, $Rand, 1);
                        $OutSeq .= $RandBase;
                        }
                push(@Out, $OutSeq);
        }
        return @Out;
}

######Sub-routine to calculate MFEs using RNAfold#################
sub Energy {

    my @engarr = @_;
    my $k = 0;
    my @returnarray = ();

    foreach my $Sequence (@engarr) {

	my $Command = "echo " . $Sequence . " | RNAfold";
	my @Out = `$Command`;
    #print @Out;
    #print "\n";
    my $STR_EN =  $Out[1];
    my @STR_EN = split(/\s+\(/, $STR_EN);
    my $EN = $STR_EN[1];
    #print "$EN\n";
    $EN =~ s/\(//g;
    $EN =~ s/\)//g;
    
    push (@returnarray, $EN);

}
return @returnarray;
@returnarray = ();
}

######Sub-routine to calculate Nucleotide frequencies#######
sub NucFreqs {

        my $InSeq = $_[0];
        $InSeq =~ s/T/U/g;

        $Gs = 0;
        $Cs = 0;
        $As = 0;
        $Us = 0;

        while ($InSeq =~ /G/g) {$Gs += 1;}
        while ($InSeq =~ /C/g) {$Cs += 1;}
        while ($InSeq =~ /A/g) {$As += 1;}
        while ($InSeq =~ /U/g) {$Us += 1;}
        
	my $length = $As + $Cs + $Gs + $Us;
        my $OutFreqs = "$As\t$Gs\t$Cs\t$Us";
        return  $OutFreqs;
}

######Sub-routine to calculate Z-scores########
sub ZScore {

	my @arr = @_;
	#print "This is my arr @arr\n";
	my $sum = 0;
	my $count = @arr;
	foreach my $l(@arr){
	

		$sum += $l;
                #print "DG: $l\n";

}
	my $Average = $sum/$count;
        my $average = ($sum-$arr[0])/($count-1);
        #print "sum: $sum\n";
	$sum = "";
	#print "AVG: $Average\n";

	my $Sigma = 0;
	my $Sum = 0;
	
	foreach my $m(@arr){
		#print "My i is $i\n";
		$Sigma = ($m - $Average)**2;
		#print "My Sigma is $Sigma\n";
		$Sum += $Sigma;
}
	my $SD = sqrt($Sum/$count);
        #print "SUM: $Sum\n";
        #print "SD: $SD\n";
	$Sum = "";

	my $return = (($SD ne 0) ? (($arr[0]-$Average)/$SD):"Undefined");
        my $ZScore = "$return\t$average";
        chomp $ZScore;
        $ZScore =~ s/\\n//g;
        #print "My return is $return\n";
		my $Output = substr($ZScore, 0, 5);
	return $Output;
}

######Sub-routine to calculate P-value (fraction of scrambled dG < native dG########
sub PValue {

	my @arr = @_;
	my $BelowNative = 0;
	my $TotalCount = @arr;
	
	my $Native = $arr[0]; 
	foreach my $l (@arr) {
	    if ($l < $Native) {$BelowNative += 1;}	    
	}
	
	my $Fraction = ($BelowNative / $TotalCount);

	
	return $Fraction;
}
