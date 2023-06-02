##############################################
# Module for evolutionary analysis version 2 #
##############################################

use strict ;
use warnings ;

# check_stop_codon
# clean_up_cds_2seq
# codonSyn
# codonSub
# codon4aa
# kaks2


# get DNA alignment of CDS sequecnes
# This uses clustal-omega
sub get_align_cds_o{
    my(@cds) = @_;
    my$n_cds = scalar@cds;
    #print"# of cds = $n_cds\n";
    # check CDS
    my@stop = ();
    foreach my$cds(@cds){
        my$l = length$cds;
        # see if lengths of CDS is a multiple of 3
        if($l%3 != 0){
            print"Error in sub get_align_cds_o\nLength\n";exit;
        }
        # keep stop codon (if avaiable)
        my$tmp = substr($cds, $l-3, 3);
        if(&check_stop_codon(\$tmp)){
            $cds = substr ($cds, 0, $l-3) ;
        }
        else {
            $tmp = "---";
        }
        push(@stop, $tmp);
    }
    #print"@stop\n";
    # translate
    my@pep = ();
    foreach my$cds(@cds){
        my$tmp = &mycodon($cds);
        push(@pep, $tmp);
    }
    #print"@pep\n";
    # align pep
    open(TMP, ">pep.fas")||die"$!";
    for(my$i=0; $i < $n_cds; $i++){
        print TMP">x$i\n$pep[$i]\n";
    }
    close (TMP);
    # run clustalo
    system("clustalo -i pep.fas -o pep.fas --force");
    # algin CDS
    my%pep = ();
    &fopen_hash_fas3(\%pep, \ "pep.fas");
    my@res = ("") x $n_cds;
    for(my$k=0; $k < $n_cds; $k++){
        # Convert into array
        my@tmpp = split("", $pep{"x$k"});
        my@tmpc = split("", $cds[$k]);
        # insert gap -> cds
        my@outc = (); 
        my$j = -1;
        foreach(@tmpp){
            if($_ eq "-"){
                push(@outc, "-"); push (@outc, "-"); push (@outc, "-");
            }
            else{
                push(@outc, $tmpc[++$j]);
                push(@outc, $tmpc[++$j]);
                push(@outc, $tmpc[++$j]);
            }
        }
        $res[$k] = join("", @outc);
    }
    return@res;
}


# hash
my(%syn_change_len) = (
	'TCA' => 1.0,	# Serine
	'TCC' => 1.0,	# Serine
	'TCG' => 1.0,	# Serine
	'TCT' => 1.0,	# Serine
	'TTC' => 1/3,	# Phenylalanine
	'TTT' => 1/3,	# Phenylalanine
	'TTA' => 2/3,	# Leucine
	'TTG' => 2/3,	# Leucine
	'TAC' => 1.0,	# Tyrosine
	'TAT' => 1.0,	# Tyrosine
	#'TAA' => '*',    # Stop
	#'TAG' => '*',    # Stop
	'TGC' => 0.5,	# Cysteine
	'TGT' => 0.5,	# Cysteine
	#'TGA' => '*',    # Stop
	'TGG' => 0.0,	# Tryptophan
	'CTA' => 4/3,	# Leucine
	'CTC' => 1.0,	# Leucine
	'CTG' => 4/3,	# Leucine
	'CTT' => 1.0,	# Leucine
	'CCA' => 1.0,	# Proline
	'CCC' => 1.0,	# Proline
	'CCG' => 1.0,	# Proline
	'CCT' => 1.0,	# Proline
	'CAC' => 1/3,	# Histidine
	'CAT' => 1/3,	# Histidine
	'CAA' => 1/3,	# Glutamine
	'CAG' => 1/3,	# Glutamine
	'CGA' => 1.5,	# Arginine
	'CGC' => 1.0,	# Arginine
	'CGG' => 4/3,	# Arginine
	'CGT' => 1.0,	# Arginine
	'ATA' => 2/3,	# Isoleucine
	'ATC' => 2/3,	# Isoleucine
	'ATT' => 2/3,	# Isoleucine
	'ATG' => 0.0,	# Methionine
	'ACA' => 1.0,	# Threonine
	'ACC' => 1.0,	# Threonine
	'ACG' => 1.0,	# Threonine
	'ACT' => 1.0,	# Threonine
	'AAC' => 1/3,	# Asparagine
	'AAT' => 1/3,	# Asparagine
	'AAA' => 1/3,	# Lysine
	'AAG' => 1/3,	# Lysine
	'AGC' => 1/3,	# Serine
	'AGT' => 1/3,	# Serine
	'AGA' => 5/6,	# Arginine
	'AGG' => 2/3,	# Arginine
	'GTA' => 1.0,	# Valine
	'GTC' => 1.0,	# Valine
	'GTG' => 1.0,	# Valine
	'GTT' => 1.0,	# Valine
	'GCA' => 1.0,	# Alanine
	'GCC' => 1.0,	# Alanine
	'GCG' => 1.0,	# Alanine
	'GCT' => 1.0,	# Alanine
	'GAC' => 1/3,	# Aspartic Acid
	'GAT' => 1/3,	# Aspartic Acid
	'GAA' => 1/3,	# Glutamic Acid
	'GAG' => 1/3,	# Glutamic Acid
	'GGA' => 1.0,	# Glycine
	'GGC' => 1.0,	# Glycine
	'GGG' => 1.0,	# Glycine
	'GGT' => 1.0,	# Glycine
) ;

my(%genetic_code4aa) = (
    
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '*',    # Stop
    'TAG' => '*',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '*',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
) ;

my(%gc4) = (    
    'TCA' => 1,    # Serine
    'TCC' => 1,    # Serine
    'TCG' => 1,    # Serine
    'TCT' => 1,    # Serine
	'AGC' => 0,    # Serine
    'AGT' => 0,    # Serine
	
    'TTC' => 0,    # Phenylalanine
    'TTT' => 0,    # Phenylalanine
    
	'TTA' => 0,    # Leucine
    'TTG' => 0,    # Leucine
	'CTA' => 1,    # Leucine
    'CTC' => 1,    # Leucine
    'CTG' => 1,    # Leucine
    'CTT' => 1,    # Leucine
	
    'TAC' => 0,    # Tyrosine
    'TAT' => 0,    # Tyrosine
	
    'TAA' => 0,    # Stop
    'TAG' => 0,    # Stop
    'TGA' => 0,    # Stop
	
	'TGC' => 0,    # Cysteine
    'TGT' => 0,    # Cysteine
    
	'TGG' => 0,    # Tryptophan
    
    'CCA' => 1,    # Proline
    'CCC' => 1,    # Proline
    'CCG' => 1,    # Proline
    'CCT' => 1,    # Proline
	
    'CAC' => 0,    # Histidine
    'CAT' => 0,    # Histidine
	
    'CAA' => 0,    # Glutamine
    'CAG' => 0,    # Glutamine
	
    'CGA' => 1,    # Arginine
    'CGC' => 1,    # Arginine
    'CGG' => 1,    # Arginine
    'CGT' => 1,    # Arginine
	'AGA' => 0,    # Arginine
    'AGG' => 0,    # Arginine
    
	'ATA' => 0,    # Isoleucine
    'ATC' => 0,    # Isoleucine
    'ATT' => 0,    # Isoleucine
    
	'ATG' => 0,    # Methionine
    
	'ACA' => 1,    # Threonine
    'ACC' => 1,    # Threonine
    'ACG' => 1,    # Threonine
    'ACT' => 1,    # Threonine
    
	'AAC' => 0,    # Asparagine
    'AAT' => 0,    # Asparagine
    
	'AAA' => 0,    # Lysine
    'AAG' => 0,    # Lysine
	
	'GTA' => 1,    # Valine
    'GTC' => 1,    # Valine
    'GTG' => 1,    # Valine
    'GTT' => 1,    # Valine
    
	'GCA' => 1,    # Alanine
    'GCC' => 1,    # Alanine
    'GCG' => 1,    # Alanine
    'GCT' => 1,    # Alanine
    
	'GAC' => 0,    # Aspartic Acid
    'GAT' => 0,    # Aspartic Acid
    
	'GAA' => 0,    # Glutamic Acid
    'GAG' => 0,    # Glutamic Acid
    
	'GGA' => 1,    # Glycine
    'GGC' => 1,    # Glycine
    'GGG' => 1,    # Glycine
    'GGT' => 1,    # Glycine
);


# Check if the codon is four-degenerated
# $res = &check_gc4(\$str);
# if True, return 1
sub check_gc4{
	use strict;
    my($cdn4) = @_;
    $$cdn4 = uc $$cdn4;
    if (exists$gc4{$$cdn4}){
        return$gc4{$$cdn4};
    }
	else{
		print"Error in  check_gc4 in evolve2\n";
        print"Bad codon: $$cdn4\n";
        exit;
    }
}

sub get_align_DNA{
	use strict;

	my($foo1, $foo2) = @_;
	$foo1 = uc$foo1;
	$foo2 = uc$foo2;
	open (TMP, ">tmp2.fas")|| die "iijjnn\n" ;
	print TMP ">a\n$foo1\n" ;
	print TMP ">b\n$foo2\n" ;
	close (TMP) ;
	system ("clustalw -align -INFILE=tmp2.fas -TYPE=DNA -OUTFILE=tmp.fas -OUTPUT=GDE -OUTORDER=INPUT -CASE=UPPER > a.out") ;
	
	my (@tmp) = &fileopen ("tmp.fas") ;
	my ($n_tmp) = scalar @tmp ;
	my (@res) = ("", "") ;
	my ($f1) = 0 ;
	foreach (@tmp){
		if ( $_ eq "#a" ){
			$f1 = 1 ;
		}
		elsif ( $_ eq "#b" ){
			$f1 = 2 ;
		}
		elsif ( $f1 == 1 ){
			$res[0] = $res[0].$_ ;
		}
		elsif ( $f1 == 2 ){
			$res[1] = $res[1].$_ ;			
		}
		else {
			print "Error in align_dna\n" ; exit ;
		}
	}
	if ( length $res[0] != length $res[1] ){
		print "Error2 in align_dna\n" ; exit; 
	}
	
	system ("rm a.out") ;
	system ("rm tmp.fas") ;
	system ("rm tmp2.fas") ;
	system ("rm tmp2.dnd") ;
	
	
	return @res ;
}
	
sub get_align_cds{
	use strict ;
	my (@cdsAC) = @_ ;
	my ($n_cdsAC) = scalar @cdsAC ;
	#print "# of cds = $n_cdsAC\n" ;
	
	# check length %3
	foreach (@cdsAC){
		if ( (length $_)%3 != 0 ){print "Error in sub get_align_cds, %3\n" ; }
	}
	
	# keep stop codon (if avaiable)
	my (@stopAC) = () ;
	foreach (@cdsAC){
		my ($lenAC) = length $_ ;
		my ($tmpAC) = substr ($_, $lenAC-3, 3) ;
		if ( &check_stop_codon (\$tmpAC) ){
			$_ = substr ($_, 0, $lenAC-3) ;
		}
		else {
			$tmpAC = "---" ;
		}
		push (@stopAC, $tmpAC) ;
	}
	#print "@stopAC\n" ;
	
	# translate
	my (@pepAC) = () ;
	foreach (@cdsAC){
		my ($tmpAC) = &mycodon ($_) ;
		push (@pepAC, $tmpAC) ;
	}
	
	# align pep
	open (TMPAC, ">pep.fas")|| die "error tmp.fas\n" ;
	for ( my ($iAC)=0; $iAC < $n_cdsAC; $iAC++ ){
		print TMPAC ">$iAC\n$pepAC[$iAC]\n" ;
	}
	close (TMPAC) ;
	system ("clustalw -align -INFILE=pep.fas -TYPE=PROTEIN -OUTFILE=pep.fas -OUTPUT=GDE -OUTORDER=INPUT -CASE=UPPER > a.out") ;
	
	my (@tmpAC) = () ;
	&fileopen2 (\@tmpAC, \"pep.fas") ;
	my ($n_tmpAC) = scalar @tmpAC ;
	@pepAC = ("") x $n_cdsAC ;
	my ($fAC1) = 0 ;
	for ( my ($iAC)=1; $iAC < $n_tmpAC; $iAC++ ){
		if ( $tmpAC[$iAC] =~ /\%/ ){
			$fAC1++ ;
		}
		else {
			$pepAC[$fAC1] = $pepAC[$fAC1].$tmpAC[$iAC] ;
		}
	}
	
	# align cds
	my (@resAC) = ("") x $n_cdsAC ;
	for ( my ($kAC)=0; $kAC < $n_cdsAC; $kAC++){
		# Divided for array
		my (@tmppAC) = split ("", $pepAC[$kAC]) ;
		my (@tmpcAC) = split ("", $cdsAC[$kAC]) ;
		
		# insert gap -> cds
		my (@outcAC) = () ; my ($jAC) = -1 ;
		foreach (@tmppAC){
			if ( $_ eq "-" ){
				push (@outcAC, "-"); push (@outcAC, "-"); push (@outcAC, "-");
			}
			else {
				push (@outcAC, $tmpcAC[++$jAC]); push (@outcAC, $tmpcAC[++$jAC]); push (@outcAC, $tmpcAC[++$jAC]);
			}
		}
		$resAC[$kAC] = join ("", @outcAC) ;
	}
	
	unless ( scalar grep (/\-\-\-/, @stopAC) == $n_cdsAC ){
		for ( my ($iAC)=0; $iAC < $n_cdsAC; $iAC++ ){
			$resAC[$iAC] = $resAC[$iAC].$stopAC[$iAC] ;
		}
	}
	
	#open (TMP, ">tmp.fas")|| die "pp\n" ;
	#print TMP "\>0\n$resAC[0]\n\>1\n$resAC[1]\n" ;
	#close (TMP) ;
	
	system ("rm a.out") ;
	system ("rm pep.fas") ;
	system ("rm pep.dnd") ;
	
	return @resAC ;
}


# check stop codon if stop, 1 return
# usage: &check_stop_codon (\$codon) 
sub check_stop_codon {
	use strict ;
	my ($hogeCSC) = @_ ;
	if ( length $$hogeCSC != 3 ){print "Error in check_stop_codon, $$hogeCSC\n"; }
	if ( $$hogeCSC eq "TAA" || $$hogeCSC eq "TAG" || $$hogeCSC eq "TGA" ){
		return 1 ;
	}
	else {
		return 0 ;
	}
}

# remove codons with gap or not ATGC sites in multiple sequences
# assume first codon is ATG and last is terminate codon.
sub clean_up_cds_multiseq {
	use strict ;
    my (@seCUM) = @_ ;
	my ($n_seCUM) = scalar @seCUM ;
	
    my ($leCUM) = length $seCUM[0] ;

	@seCUM = map {uc $_} @seCUM ;

    my (@resCUM) = ("")x$n_seCUM ;
	my (@ntCUM1) = () ;
	my (@ntCUM2) = () ;
	my ($check) ;

    for ( my ($pCUM)=3; $pCUM < $leCUM-3; $pCUM+=3 ){
		@ntCUM1 = map {substr ($_, $pCUM, 3)} @seCUM ;
		@ntCUM2 = &union (@ntCUM1) ;
	
		if ( scalar grep (!/[ATGC][ATGC][ATGC]/, @ntCUM2) > 0 ){next ; }
		$check = 0 ;
		foreach my $ntCUM3 (@ntCUM2){
			if ( &check_stop_codon (\$ntCUM3) ){
				$check = 1 ; last ;
			}
		}
		if ( $check == 1 ){next ; }
		
		for ( my ($qCUM)=0; $qCUM < $n_seCUM; $qCUM++ ){
			$resCUM[$qCUM] = $resCUM[$qCUM].$ntCUM1[$qCUM] ;
		}
    }

	for ( my ($qCUM)=0; $qCUM < $n_seCUM; $qCUM++ ){
		if ( (length $resCUM[$qCUM])%3 != 0 ){print "Error in claen up cds multiseq\n" ; exit ; }
	}

    return @resCUM ;

}

# remove codons with gap or not ATGC sites in two sequences
# assume first codon is ATG and last is terminate codon.
sub clean_up_cds_2seq{
	use strict;
    my@se = @_;
    my$le = length$se[0];
    # Convert upper case
	@se = map{uc $_}@se;
    # remove unnecessary codon
    my@res = ("", "");
    my$nt1;
    my$nt2;
    for(my$p=0; $p < $le; $p+=3){
		$nt1 = substr($se[0], $p, 3);
		$nt2 = substr($se[1], $p, 3);
        # remove stop codon
		if(&check_stop_codon(\$nt1)){next;}
		if(&check_stop_codon(\$nt2)){next;}
        # other than ATGC
		unless($nt1 =~ /^[ATGC][ATGC][ATGC]$/){next;}
		unless($nt2 =~ /^[ATGC][ATGC][ATGC]$/){next;}
		$res[0] = $res[0] . $nt1;
		$res[1] = $res[1] . $nt2;
    }
    my$lk1 = length$res[0];
    my$lk2 = length$res[1];
    if($lk1%3 != 0){print "Error in claen up cds1, $lk1\n"; exit;}
    if($lk2%3 != 0){print "Error in claen up cds2, $lk2\n"; exit;}
    if($lk1 != $lk2){print"Error in clean up cds3 $lk1 $lk1\n"; exit;}
    return @res;

}

# the number of potentially synonymous changes for codon
sub codonSyn {
	use strict;
	my($codonS) = @_;
	$codonS = uc $codonS;

	if (exists $syn_change_len{$codonS}){
		return $syn_change_len{$codonS};
	}
	else {
		print "Error in codonSyn, unknown or stop codon, $codonS\n"; exit;
	}
	
}

# number of total substituion in codon
sub codonSub {
	use strict ;
	my (@hogeSU) = @_ ;
	my ($n_hogeSU) = scalar @hogeSU ;
	if ($n_hogeSU != 2){print "Error number of array in codonSub, @hogeSU\n" ; exit ; }
	
	my ($lenSU1) = length $hogeSU[0] ;
	my ($lenSU2) = length $hogeSU[1] ;
	if ($lenSU1 != $lenSU2 || $lenSU1 != 3){print "Error length of codon in codonSub, @hogeSU\n" ; exit ; }
	
	my ($difSU) = 0 ;
	for ( my ($qw)=0; $qw < $lenSU1; $qw++ ){
		if ( substr ($hogeSU[0], $qw, 1) ne substr ($hogeSU[1], $qw, 1) ){
			$difSU++;
		}
	}
	
	if ( $difSU > 3 ){print "Error # of substitution > 3 in codonSub, @hogeSU\n" ; exit ; }
	
	return $difSU ;
	
}

sub codon4aa {
	use strict ;
    my($codon4aa) = @_ ;

    $codon4aa = uc $codon4aa ;

	if ( exists $genetic_code4aa{$codon4aa}) {
        return $genetic_code4aa{$codon4aa};
    }
	else{ print "Bad codon $codon4aa\n" ; exit ; }
}

sub codonNG1 {
	use strict ;
	my (@hogeNG1) = @_;
	my ($n_hogeNG1) = scalar @hogeNG1 ;
	if ($n_hogeNG1 != 2){print "Error number of array in codonNG1, @hogeNG1\n"; exit; }
	
	my ($pepNG1) = &codon4aa ($hogeNG1[0]) ;
	my ($pepNG2) = &codon4aa ($hogeNG1[1]) ;
	
	if ( $pepNG1 eq "*" || $pepNG2 eq "*" ){print "Stop codon are present in codonNG1, @hogeNG1\n" ; exit ; }

	if ( $pepNG1 ne $pepNG2 ){
		return (0, 1) ;
	}
	elsif ( $pepNG1 eq $pepNG2 ){
		return (1, 0) ;
	}
	else {print "Error in codonNG1, $pepNG1, $pepNG2\n" ; exit ; }
}

sub codonNG2 {
	use strict ;
	my (@hogeNG2) = @_ ;
	my ($n_hogeNG2) = scalar @hogeNG2 ;
	if ( $n_hogeNG2 != 2 ){print "Error number of array in codonNG2, @hogeNG2\n" ; exit ; }
	
	# position of 2 mutaion
	my (@posiNG2) = (-1, -1) ;
	my ($jjNG2) = 0 ;
	if ( substr ($hogeNG2[0], 0, 1) ne substr ($hogeNG2[1], 0, 1) ){$posiNG2[$jjNG2] = 0 ; $jjNG2++ ;}
	if ( substr ($hogeNG2[0], 1, 1) ne substr ($hogeNG2[1], 1, 1) ){$posiNG2[$jjNG2] = 1 ; $jjNG2++ ;}
	if ( substr ($hogeNG2[0], 2, 1) ne substr ($hogeNG2[1], 2, 1) ){$posiNG2[$jjNG2] = 2 ; $jjNG2++ ;}
	if ( $jjNG2 != 2|| $posiNG2[0] == -1 || $posiNG2[1] == -1 ){print "No two mutations, @hogeNG2\n" ; exit ; }
	#print "@posiNG2\n";
	
	# create intermediate codon
	# 1
	my (@tmpNG21) = split ("", $hogeNG2[0]) ;
	my (@tmpNG22) = @tmpNG21 ;
	$tmpNG21[$posiNG2[0]] = substr ($hogeNG2[1], $posiNG2[0], 1) ; 
	my ($mmCod1) = join ("", @tmpNG21) ;
	# 2
	$tmpNG22[$posiNG2[1]] = substr ($hogeNG2[1], $posiNG2[1], 1) ;
	my ($mmCod2) = join ("", @tmpNG22) ;
	#print "$mmCod1, $mmCod2\n";
	
	# check stop codon
	my ($stopNG21) = 0 ;
	my ($stopNG22) = 0 ;
	if ( $mmCod1 eq "TAA" || $mmCod1 eq "TAG" || $mmCod1 eq "TGA" ){$stopNG21 = 1 ; }
	if ( $mmCod2 eq "TAA" || $mmCod2 eq "TAG" || $mmCod2 eq "TGA" ){$stopNG22 = 1 ; }
	#print "$stopNG21, $stopNG22\n";
	if ( $stopNG21 + $stopNG22 == 2 ){print "Both pathway include stop codon in codonNG2, @hogeNG2\n"; exit; }
	
	# count syn and rep
	my ($pepNG21) =  &codon4aa ($hogeNG2[0]) ;
	my ($pepNG23) =  &codon4aa ($hogeNG2[1]) ;
	my (@difNG2) = (0, 0) ;
	my ($numPathNG2) = 0 ;
	# path1
	if ($stopNG21  == 0){
		my ($pepNG22) =  &codon4aa ($mmCod1) ;
		#print "$pepNG21, $pepNG22, $pepNG23\n";
		   if ( $pepNG21 eq $pepNG22 ){$difNG2[0]++ ; }
		elsif ( $pepNG21 ne $pepNG22 ){$difNG2[1]++ ; }
		   if ( $pepNG22 eq $pepNG23 ){$difNG2[0]++ ; }
		elsif ( $pepNG22 ne $pepNG23 ){$difNG2[1]++ ; }
		#print "@difNG2\n";
		$numPathNG2++;
	}
	# path2
	if ( $stopNG22  == 0 ){
		my ($pepNG24) =  &codon4aa ($mmCod2) ;
		#print "$pepNG21, $pepNG24, $pepNG23\n";
		   if ( $pepNG21 eq $pepNG24 ){$difNG2[0]++ ; }
		elsif ( $pepNG21 ne $pepNG24 ){$difNG2[1]++ ; }
		   if ( $pepNG24 eq $pepNG23 ){$difNG2[0]++ ; }
		elsif ( $pepNG24 ne $pepNG23 ){$difNG2[1]++ ; }
		#print "@difNG2\n";
		$numPathNG2++ ;
	}
	@difNG2 = map {$_/$numPathNG2} @difNG2 ;
	#print "@difNG2\n";
	return @difNG2 ;

}

sub codonNG3 {
	use strict;
	my (@hogeNG3) = @_ ;
	my ($n_hogeNG3) = scalar @hogeNG3 ;
	if ( $n_hogeNG3 != 2 ){print "Error number of array in codonNG3, @hogeNG3\n" ; exit ; }
	
	# count syn and rep
	my (@tmpNG31) = split ("", $hogeNG3[0]) ;
	my (@tmpNG32) = split ("", $hogeNG3[1]) ;
	my ($pepNG31) =  &codon4aa ($hogeNG3[0]) ;
	my ($pepNG32) =  &codon4aa ($hogeNG3[1]) ;
	my (@difNG3) = (0, 0) ;
	my ($numPathNG3) = 0 ;
	my ($mmCod31) ;
	my ($mmCod32) ;
	
	# path1
	$mmCod31 = "$tmpNG32[0]$tmpNG31[1]$tmpNG31[2]" ;
	$mmCod32 = "$tmpNG32[0]$tmpNG32[1]$tmpNG31[2]" ;
	#print "$mmCod31, $mmCod32\n";
	if ( $mmCod31 ne "TAA" && $mmCod31 ne "TAG" && $mmCod31 ne "TGA" && $mmCod32 ne "TAA" && $mmCod32 ne "TAG" && $mmCod32 ne "TGA" ){
		my ($mmPep31) = &codon4aa ($mmCod31);
		my ($mmPep32) = &codon4aa ($mmCod32);
		#print "$mmPep31, $mmPep32\n";
		   if ($pepNG31 eq $mmPep31){$difNG3[0]++; }
		elsif ($pepNG31 ne $mmPep31){$difNG3[1]++; }
		   if ($mmPep31 eq $mmPep32){$difNG3[0]++; }
		elsif ($mmPep31 ne $mmPep32){$difNG3[1]++; }
		   if ($mmPep32 eq $pepNG32){$difNG3[0]++; }
		elsif ($mmPep32 ne $pepNG32){$difNG3[1]++; }
		$numPathNG3++;
	}
	# path2
	$mmCod31 = "$tmpNG32[0]$tmpNG31[1]$tmpNG31[2]";
	$mmCod32 = "$tmpNG32[0]$tmpNG31[1]$tmpNG32[2]";
	#print "$mmCod31, $mmCod32\n";
	if ($mmCod31 ne "TAA" && $mmCod31 ne "TAG" && $mmCod31 ne "TGA" && $mmCod32 ne "TAA" && $mmCod32 ne "TAG" && $mmCod32 ne "TGA"){
		my ($mmPep31) = &codon4aa ($mmCod31);
		my ($mmPep32) = &codon4aa ($mmCod32);
		#print "$mmPep31, $mmPep32\n";
		   if ($pepNG31 eq $mmPep31){$difNG3[0]++; }
		elsif ($pepNG31 ne $mmPep31){$difNG3[1]++; }
		   if ($mmPep31 eq $mmPep32){$difNG3[0]++; }
		elsif ($mmPep31 ne $mmPep32){$difNG3[1]++; }
		   if ($mmPep32 eq $pepNG32){$difNG3[0]++; }
		elsif ($mmPep32 ne $pepNG32){$difNG3[1]++; }
		$numPathNG3++;
	}
	# path3
	$mmCod31 = "$tmpNG31[0]$tmpNG32[1]$tmpNG31[2]";
	$mmCod32 = "$tmpNG32[0]$tmpNG32[1]$tmpNG31[2]";
	#print "$mmCod31, $mmCod32\n";
	if ($mmCod31 ne "TAA" && $mmCod31 ne "TAG" && $mmCod31 ne "TGA" && $mmCod32 ne "TAA" && $mmCod32 ne "TAG" && $mmCod32 ne "TGA"){
		my ($mmPep31) = &codon4aa ($mmCod31);
		my ($mmPep32) = &codon4aa ($mmCod32);
		#print "$mmPep31, $mmPep32\n";
		   if ($pepNG31 eq $mmPep31){$difNG3[0]++; }
		elsif ($pepNG31 ne $mmPep31){$difNG3[1]++; }
		   if ($mmPep31 eq $mmPep32){$difNG3[0]++; }
		elsif ($mmPep31 ne $mmPep32){$difNG3[1]++; }
		   if ($mmPep32 eq $pepNG32){$difNG3[0]++; }
		elsif ($mmPep32 ne $pepNG32){$difNG3[1]++; }
		$numPathNG3++;
	}
	# path4
	$mmCod31 = "$tmpNG31[0]$tmpNG32[1]$tmpNG31[2]";
	$mmCod32 = "$tmpNG31[0]$tmpNG32[1]$tmpNG32[2]";
	#print "$mmCod31, $mmCod32\n";
	if ($mmCod31 ne "TAA" && $mmCod31 ne "TAG" && $mmCod31 ne "TGA" && $mmCod32 ne "TAA" && $mmCod32 ne "TAG" && $mmCod32 ne "TGA"){
		my ($mmPep31) = &codon4aa ($mmCod31);
		my ($mmPep32) = &codon4aa ($mmCod32);
		#print "$mmPep31, $mmPep32\n";
		   if ($pepNG31 eq $mmPep31){$difNG3[0]++; }
		elsif ($pepNG31 ne $mmPep31){$difNG3[1]++; }
		   if ($mmPep31 eq $mmPep32){$difNG3[0]++; }
		elsif ($mmPep31 ne $mmPep32){$difNG3[1]++; }
		   if ($mmPep32 eq $pepNG32){$difNG3[0]++; }
		elsif ($mmPep32 ne $pepNG32){$difNG3[1]++; }
		$numPathNG3++;
	}
	# path5
	$mmCod31 = "$tmpNG31[0]$tmpNG31[1]$tmpNG32[2]";
	$mmCod32 = "$tmpNG32[0]$tmpNG31[1]$tmpNG32[2]";
	#print "$mmCod31, $mmCod32\n";
	if ($mmCod31 ne "TAA" && $mmCod31 ne "TAG" && $mmCod31 ne "TGA" && $mmCod32 ne "TAA" && $mmCod32 ne "TAG" && $mmCod32 ne "TGA"){
		my ($mmPep31) = &codon4aa ($mmCod31);
		my ($mmPep32) = &codon4aa ($mmCod32);
		#print "$mmPep31, $mmPep32\n";
		   if ($pepNG31 eq $mmPep31){$difNG3[0]++; }
		elsif ($pepNG31 ne $mmPep31){$difNG3[1]++; }
		   if ($mmPep31 eq $mmPep32){$difNG3[0]++; }
		elsif ($mmPep31 ne $mmPep32){$difNG3[1]++; }
		   if ($mmPep32 eq $pepNG32){$difNG3[0]++; }
		elsif ($mmPep32 ne $pepNG32){$difNG3[1]++; }
		$numPathNG3++;
	}
	# path5
	$mmCod31 = "$tmpNG31[0]$tmpNG31[1]$tmpNG32[2]";
	$mmCod32 = "$tmpNG31[0]$tmpNG32[1]$tmpNG32[2]";
	#print "$mmCod31, $mmCod32\n";
	if ($mmCod31 ne "TAA" && $mmCod31 ne "TAG" && $mmCod31 ne "TGA" && $mmCod32 ne "TAA" && $mmCod32 ne "TAG" && $mmCod32 ne "TGA"){
		my ($mmPep31) = &codon4aa ($mmCod31);
		my ($mmPep32) = &codon4aa ($mmCod32);
		#print "$mmPep31, $mmPep32\n";
		   if ($pepNG31 eq $mmPep31){$difNG3[0]++; }
		elsif ($pepNG31 ne $mmPep31){$difNG3[1]++; }
		   if ($mmPep31 eq $mmPep32){$difNG3[0]++; }
		elsif ($mmPep31 ne $mmPep32){$difNG3[1]++; }
		   if ($mmPep32 eq $pepNG32){$difNG3[0]++; }
		elsif ($mmPep32 ne $pepNG32){$difNG3[1]++; }
		$numPathNG3++;
	}
	#print "@difNG3, $numPathNG3\n";
	@difNG3 = map {$_/$numPathNG3} @difNG3;
	#print "@difNG3, $numPathNG3\n";
	if ($numPathNG3 == 0){print "No pathway for stop codon in codonNG3, @hogeNG3\n"; exit; }
	return @difNG3 ;
}



# Nei-Gojobiri method ver. 2
# before running, start, stop, gap and not-ATGC sites should be removed
# @res = &kaks2(@($str1, $str2));
sub kaks2{
	use strict ;
	my(@seqN) = @_;
	@seqN = map uc, @seqN;
	# cut gap after checking frameshift mutation
	my$ltmp1 = length$seqN[0];
	my$ltmp2 = length$seqN[1];
	if($ltmp1 != $ltmp2){
        print"Error in kaks2 in evolve2\n";
        print"Length of input seqs are different, $ltmp1, $ltmp2\n"; 
        exit;
    }
	if($ltmp1%3 != 0 || $ltmp2%3 != 0){
        print"Error in kaks2 in evolve2\n";
        print"Length of seq is not a multiple of 3, $ltmp1, $ltmp2\n";
        exit;
    }
	my$LL = $ltmp1 ;
	# the number of potentially synonymous changes for codon
	my$Lsyn1 = 0;
	my$Lsyn2 = 0;
	my$Lstmp;
	# count synonymous and nonsynonymous changes
	my ($dSS) = 0;
	my ($dNN) = 0;
	my ($ccNN) ;
	for (my$i=0; $i < $LL; $i+=3){
		# get codon
		my$tN1 = substr($seqN[0], $i, 3);
		my$tN2 = substr($seqN[1], $i, 3);
		# synonymous change sites
		$Lstmp = &codonSyn("$tN1");
		$Lsyn1 = $Lsyn1 + $Lstmp;
		$Lstmp = &codonSyn("$tN2");
		$Lsyn2 = $Lsyn2 + $Lstmp;
		if ( $tN1 eq $tN2 ){next;} 
		$ccNN = &codonSub (($tN1, $tN2)) ;
		my (@difSN) ;
		# number of substitution = 1
		if ( $ccNN == 1 ){
			@difSN = &codonNG1 (($tN1, $tN2));
		}
		# 2
		elsif ($ccNN == 2){
			@difSN = &codonNG2 (($tN1, $tN2));
		}
		# 3
		elsif ($ccNN == 3){
			@difSN = &codonNG3 (($tN1, $tN2));
		}
		else {print "Error in kaks2\n"; exit; }
		$dSS += $difSN[0] ;
		$dNN += $difSN[1] ;
	}
	my$Lrep1 = $LL - $Lsyn1;
	my$Lrep2 = $LL - $Lsyn2;
	my$SS = ($Lsyn1+$Lsyn2)/2;
	my$NN = ($Lrep1+$Lrep2)/2;
	my$pSS;
	my$pNN;
	   if ($SS != 0){$pSS = $dSS/$SS;}
	elsif ($SS == 0){$pSS = 0;}
	   if ($NN != 0){$pNN = $dNN/$NN;}
	elsif ($NN == 0){$pNN = 0;}
	#print "$SS, $dSS, $pSS\n$NN, $dNN, $pNN\n";
	# 0 -> # of syn site
	# 1 -> # of syn change
	# 2 -> # of syn dif per site
	# 3 -> # of rep site
	# 4 -> # of rep change
	# 5 -> # of rep dif per site
	my$resNG = "$SS\t$dSS\t$pSS\t$NN\t$dNN\t$pNN";
	return$resNG;
}


# Apply Jukes & Cantor (1969) correlation
# $res = &jc69($float);
sub jc69 {
	use strict;
	my($hoge) = @_;
	if($hoge == 0){
		return 0;
	}
	else{
		return -3/4 * log (1-4/3*$hoge);
	}
}

sub lensynsites22 {
	use strict;
	my ($seqLN) = @_;
	$seqLN = uc $seqLN;
	#print "$seqN[0]\n\n$seqN[1]\n";

	# cut gap after checking frameshift mutation
	my ($ltmp11) = length $seqLN;
	if ($ltmp11%3 != 0){print "Length of seq is not 3 times in lensynsites, $ltmp11\n"; exit; }

	# the number of potentially synonymous changes for codon
	my ($LsynL1) = 0;
	for (my ($iL)=0; $iL < $ltmp11; $iL+=3){
		my ($tN1) = substr ($seqLN, $iL, 3);
		my ($Lstmp) = &codonSyn ("$tN1");
		$LsynL1 = $LsynL1 + $Lstmp;
		#print "$Lsyn1, $Lsyn2\n";
	}
	#print "$Lsyn1, $Lrep1, $Lsyn2, $Lrep2\n";

	return $LsynL1;
}

1;
