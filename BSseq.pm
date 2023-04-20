####################
# Module for BSseq #
####################

# List of sub
	


# count C and T on plus strand, paired-end results
# args: result hash, chrom, short read 1, position 1 (zero-start), short read 2, position 2 (zero-start)  
sub count_CT_plus_bismark_pe {
	use strict ;
	my ($resCT1, $conCT1, $tmCT1, $posCT1, $tmCT2, $posCT2) = @_ ;
	#print "$$tmCT1\n" ;
	my ($ntCT1, $ntCT2, $posCTCT) ;
	my (%checkCT1) = () ;
	my (@w1CT1) ;
	
	# 5' mate
	$$tmCT1 = uc $$tmCT1 ;
	my ($lCT1) = length $$tmCT1 ;
	my ($subCT1) = uc substr ($$conCT1, $$posCT1, $lCT1) ;
	my ($lCT2) = length $subCT1 ;
	#print "$$tmCT1\n$subCT1\n" ;
	
	for ( my ($yCT1)=0; $yCT1 < $lCT2; $yCT1++ ){
		$ntCT1 = substr ($$tmCT1, $yCT1, 1) ; # BSseq	
		$ntCT2 = substr ($subCT1, $yCT1, 1) ; # con
		
		# this ver doesn't allow any non-bisulfite mismatches and indels
		#if ( $ntCT1 ne $ntCT2 && $ntCT1 ne "T" && $ntCT2 ne "C" ){print "  1 BS-seq $ntCT1 1con $ntCT2\n" ; }
		if ( $ntCT2 ne "C" && $ntCT1 ne $ntCT2 ){last ; }
		if ( $ntCT2 eq "C" && $ntCT1 ne "C" && $ntCT1 ne "T"){last ; }
		if ( $ntCT2 ne "C" ){next ; }
		
		$posCTCT = $$posCT1+$yCT1 ;
		$checkCT1{$posCTCT} = 1 ;
		
		if ( exists $$resCT1{$posCTCT} ){
			@w1CT1 = split ("\t", $$resCT1{$posCTCT}) ;
			$w1CT1[0]++ ;
			if ( $ntCT1 eq "C" ){$w1CT1[1]++ ; }
			elsif ( $ntCT1 eq "T" ){}
			else {print "bar3 $ntCT1 $ntCT2\n" ; exit ; }
			$$resCT1{$posCTCT} = join ("\t", @w1CT1) ;
			
		}
		else {
			if ( $ntCT1 eq "C" ){
				$$resCT1{$posCTCT} = "1\t1\t+" ;
			}
			elsif ( $ntCT1 eq "T" ){
				$$resCT1{$posCTCT} = "1\t0\t+" ;
			}
			else {print "bar1 $ntCT1 $ntCT2\n" ; exit ; }
		}
	} # $yCT1
	
	# 3' mate
	$$tmCT2 = uc $$tmCT2 ;
	$lCT1 = length $$tmCT2 ;
	$subCT1 = uc substr ($$conCT1, $$posCT2, $lCT1) ;
	$lCT2 = length $subCT1 ;
	#print "$$tmCT2\n$subCT1\n" ;
	
	for ( my ($yCT1)=0; $yCT1 < $lCT2; $yCT1++ ){
		$ntCT1 = substr ($$tmCT2, $yCT1, 1) ; # BSseq	
		$ntCT2 = substr ($subCT1, $yCT1, 1) ; # con
		
		# this ver doesn't allow any non-bisulfite mismatches and indels
		#if ( $ntCT1 ne $ntCT2 && $ntCT1 ne "T" && $ntCT2 ne "C" ){print "  2 BS-seq $ntCT1 1con $ntCT2\n" ; }
		if ( $ntCT2 ne "C" && $ntCT1 ne $ntCT2 ){last ; }
		if ( $ntCT2 eq "C" && $ntCT1 ne "C" && $ntCT1 ne "T"){last ; }
		if ( $ntCT2 ne "C" ){next ; }
		
		$posCTCT = $$posCT2+$yCT1 ;
		if ( exists $checkCT1{$posCTCT} ){next ; }
		
		if ( exists $$resCT1{$posCTCT} ){
			@w1CT1 = split ("\t", $$resCT1{$posCTCT}) ;
			$w1CT1[0]++ ;
			if ( $ntCT1 eq "C" ){$w1CT1[1]++ ; }
			elsif ( $ntCT1 eq "T" ){}
			else {print "bar3 $ntCT1 $ntCT2\n" ; exit ; }
			$$resCT1{$posCTCT} = join ("\t", @w1CT1) ;
			
		}
		else {
			if ( $ntCT1 eq "C" ){
				$$resCT1{$posCTCT} = "1\t1\t+" ;
			}
			elsif ( $ntCT1 eq "T" ){
				$$resCT1{$posCTCT} = "1\t0\t+" ;
			}
			else {print "bar1 $ntCT1 $ntCT2\n" ; exit ; }
		}
	} # $yCT1
	
}

# count C and T on plus strand, paired-end results
# args: result hash, chrom, short read 1, position 1 (zero-start), short read 2, position 2 (zero-start)  
sub count_CT_minus_bismark_pe {
	use strict ;
	my ($resCT1, $conCT1, $tmCT1, $posCT1, $tmCT2, $posCT2) = @_ ;
	#print "$$tmCT1\n" ;
	my ($ntCT1, $ntCT2, $posCTCT) ;
	my (%checkCT1) = () ;
	my (@w1CT1) ;
	
	# 5' mate
	$$tmCT1 = uc $$tmCT1 ;
	my ($lCT1) = length $$tmCT1 ;
	my ($subCT1) = uc substr ($$conCT1, $$posCT1, $lCT1) ;
	my ($lCT2) = length $subCT1 ;
	#print "$$tmCT1\n$subCT1\n" ;
	
	for ( my ($yCT1)=0; $yCT1 < $lCT2; $yCT1++ ){
		$ntCT1 = substr ($$tmCT1, $yCT1, 1) ; # BSseq	
		$ntCT2 = substr ($subCT1, $yCT1, 1) ; # con
		
		# this ver doesn't allow any non-bisulfite mismatches and indels
		#if ( $ntCT1 ne $ntCT2 && $ntCT1 ne "T" && $ntCT2 ne "C" ){print "  1 BS-seq $ntCT1 1con $ntCT2\n" ; }
		if ( $ntCT2 ne "G" && $ntCT1 ne $ntCT2 ){last ; }
		if ( $ntCT2 eq "G" && $ntCT1 ne "G" && $ntCT1 ne "A"){last ; }
		if ( $ntCT2 ne "G" ){next ; }
		
		$posCTCT = $$posCT1+$yCT1 ;
		$checkCT1{$posCTCT} = 1 ;
		
		if ( exists $$resCT1{$posCTCT} ){
			@w1CT1 = split ("\t", $$resCT1{$posCTCT}) ;
			$w1CT1[0]++ ;
			if ( $ntCT1 eq "G" ){$w1CT1[1]++ ; }
			elsif ( $ntCT1 eq "A" ){}
			else {print "bar3 $ntCT1 $ntCT2\n" ; exit ; }
			$$resCT1{$posCTCT} = join ("\t", @w1CT1) ;
			
		}
		else {
			if ( $ntCT1 eq "G" ){
				$$resCT1{$posCTCT} = "1\t1\t-" ;
			}
			elsif ( $ntCT1 eq "A" ){
				$$resCT1{$posCTCT} = "1\t0\t-" ;
			}
			else {print "bar1 $ntCT1 $ntCT2\n" ; exit ; }
		}
	} # $yCT1
	
	# 3' mate
	$$tmCT2 = uc $$tmCT2 ;
	$lCT1 = length $$tmCT2 ;
	$subCT1 = uc substr ($$conCT1, $$posCT2, $lCT1) ;
	$lCT2 = length $subCT1 ;
	#print "$$tmCT2\n$subCT1\n" ;
	
	for ( my ($yCT1)=0; $yCT1 < $lCT2; $yCT1++ ){
		$ntCT1 = substr ($$tmCT2, $yCT1, 1) ; # BSseq	
		$ntCT2 = substr ($subCT1, $yCT1, 1) ; # con
		
		# this ver doesn't allow any non-bisulfite mismatches and indels
		#if ( $ntCT1 ne $ntCT2 && $ntCT1 ne "T" && $ntCT2 ne "C" ){print "  2 BS-seq $ntCT1 1con $ntCT2\n" ; }
		if ( $ntCT2 ne "G" && $ntCT1 ne $ntCT2 ){last ; }
		if ( $ntCT2 eq "G" && $ntCT1 ne "G" && $ntCT1 ne "A"){last ; }
		if ( $ntCT2 ne "G" ){next ; }
		
		$posCTCT = $$posCT2+$yCT1 ;
		if ( exists $checkCT1{$posCTCT} ){next ; }
		
		if ( exists $$resCT1{$posCTCT} ){
			@w1CT1 = split ("\t", $$resCT1{$posCTCT}) ;
			$w1CT1[0]++ ;
			if ( $ntCT1 eq "G" ){$w1CT1[1]++ ; }
			elsif ( $ntCT1 eq "A" ){}
			else {print "bar3 $ntCT1 $ntCT2\n" ; exit ; }
			$$resCT1{$posCTCT} = join ("\t", @w1CT1) ;
			
		}
		else {
			if ( $ntCT1 eq "G" ){
				$$resCT1{$posCTCT} = "1\t1\t-" ;
			}
			elsif ( $ntCT1 eq "A" ){
				$$resCT1{$posCTCT} = "1\t0\t-" ;
			}
			else {print "bar1 $ntCT1 $ntCT2\n" ; exit ; }
		}
	} # $yCT1
	
	
}




1;