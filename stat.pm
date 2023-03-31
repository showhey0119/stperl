#########################
# Module for statistics #
#########################

# List of sub
# sum
# mean
# sd
# cov
# cc
# max
# min
# permutation
# columnTab
# union
# column
# columnSeq
# position
# log10
# log2
# factorial (k!)


# Calculate summation
sub sum{
    use strict;
    my(@tmp) = @_;	# argument of sample
    my$res = 0;
    foreach my$tmp(@tmp){
        $res += $tmp;
    }
    return $res;
}


# Calculate mean
sub mean{
    use strict;
    my(@foo) = @_;           # argument of sample
    my$num = scalar@foo;  # number of sample
    my$sigma = 0;          # sum of sample
    foreach my $tmp(@foo){
        $sigma += $tmp;
    }
    return $sigma / $num;
}


# Calculate unbiased standard deviation.
sub sd{
    use strict;
    my(@foo) = @_;            # argument of sample
    my$num = scalar@foo;   # number of sample
    my$mean = &mean(@foo); # mean of sample
    my$sd = 0;              #standard deviation of @foo
    for(my$i=0; $i < $num; $i++){
        $sd += ($foo[$i]-$mean) ** 2;
    }
    $sd = sqrt($sd/($num-1));
    return $sd;
}


# covariance
sub cov{
    use strict;	
    my($foo1, $foo2) = @_; # argument of two sample
    my$num1 = scalar@$foo1; # number of each sample
    my$num2 = scalar@$foo2; 
    # check if array sizes are the same
    unless($num1 == $num2){
        print"Error in cov in stat\n";
        print"Array sizes are different\n"; 
        exit;	
    }
    # calculate covariance
    my$cov = 0;
    my$mean1 = &mean(@$foo1);
    my$mean2 = &mean(@$foo2);
    for(my$i=0; $i < $num1; $i++){
        $cov += ($$foo1[$i]-$mean1) * ($$foo2[$i]-$mean2);
    }
    return $cov / ($num1-1);
}


# Calculate correlation coefficient
sub cc{
    use strict;
    my($foo1,$foo2) = @_;
    my$num1 = scalar@$foo1;
    my$num2 = scalar@$foo2; 
    # check if array sizes are the same
    unless($num1 == $num2){
        print"Error in cc in stat\n";
        print"Array sizes are different\n";
        exit;
    }
    my$cov = &cov(\@$foo1, \@$foo2);
    my$sd1 = &sd(@$foo1);
    my$sd2 = &sd(@$foo2);
    return $cov / ($sd1 * $sd2);
}


# Return the max value
sub max{
    use strict;
    my(@hoge) = @_;
    my$max = $hoge[0];
    foreach(@hoge){
        if($max < $_){
            $max = $_;
        }
    }
    return $max;
}


# Return the  min value
sub min{
	use strict;
	my(@hoge) = @_;
	my$min = $hoge[0];
	foreach(@hoge){
		if($min > $_){
			$min = $_;
		}
	}
	return $min;
}


# Permutate the array
sub permutation{
	use strict;
	my(@hoge) = @_;
	my$n_hoge = scalar @hoge;
	my@out = ();	
	for(my$i=0; $i < $n_hoge; $i++){
		my$r = int(rand($#hoge+1));
		push(@out, $hoge[$r]);
		splice(@hoge, $r, 1);
	}
	return @out ;
}


# Extract a column from a tab array.
# Two argues are required
# e.g., @res = &columnTab(\@tmp, \2);
sub columnTab{
	use strict;
	my($hoge, $foo) = @_;
	my$n_hoge = scalar @$hoge;
	my@out = ();
	foreach(@$hoge){
		my@tmp = split(/\s+/, $_);
		push (@out, $tmp[$$foo]);
	}
	return @out;
}


# Union the array
sub union{
	use strict;
	my(@hoge) = @_;
	my%seen;
	@hoge = grep (!$seen{$_}++, @hoge);
	return @hoge;
}


# Extract a column from the 2D array
# Two argues are required
# e.g., @res = &column(\@tmp, \2);
sub column{
	use strict;
	my($hoge, $foo) = @_;
	my$n_hoge = scalar@$hoge;
	my$c_hoge = scalar@{@$hoge[0]};
	if($$foo >= $c_hoge){
        print"Error in column in stat\n";
        print "The second arg is bigger than the number of columns\n";
        exit;
    }
	my@out = ();
	for(my$i=0; $i < $n_hoge; $i++){
		push(@out, $$hoge[$i][$$foo]);
	}
	return @out;
}


# Extract characters of the position of str
# Two argues are required
# e.g., @res = &columnSeq(\@tmp, \2);
sub columnSeq {
	use strict;
	my($hoge, $foo) = @_;
	my@out = map{substr($_, $$foo, 1)}@$hoge;
	return @out;
}


# Return the first index of a str  in the array.
# Two argues are required
# Note that this sub could not be used for numbers.
sub position{
	use strict;
	my($hoge, $foo) = @_;
	my$cc = 0;
	foreach my$tmp(@$hoge){
		if($tmp eq $$foo){
			return $cc;
		}
		$cc++;
	}
	# not fount
	return -1;
}


# Calculate log10
sub log10{
	use strict;
    my($x) = @_;
    return log($x) / log (10);
}


# Calculate log2
sub log2{
	use strict;
	my($x) = @_ ;
	return log($x) / log (2);
}


# Calculate median
sub median{
	use strict ;
	return unless @_;
	return $_[0] unless @_ > 1;
	@_= sort{$a<=>$b}@_;
	return$_[$#_/2] if @_&1;
	my$mid= @_/2;
	return($_[$mid-1]+$_[$mid])/2;
}


# Calculate the factorial of arg
sub factorial {
    use strict;
	my($fac) = @_;
	my$res_fac = 1;
	if($fac == 0){
		return 1;
	}
	else{
		for( my$i=1; $i <= $fac; $i++ ){
			$res_fac *= $i;
		}
		return $res_fac ;
	}
}


1;
