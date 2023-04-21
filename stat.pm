#########################
# Module for statistics #
#########################

# add "use stat;" to your perl script

# List of functions
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
# $res = &sum(@array);
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
# $res = &mean(@array);
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
# $res = &sd(@array);
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
# $res = &cov(\@array1, \@array2);
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
# $res = &cc(\@array1, \@array2);
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
# $res = &max(@array);
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
# $res = &min(@array);
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
# @res = &permutation(@array);
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


# Extract $int th column from a tab array.
# @res = &columnTab(\@array, \$int);
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
# @res = &union(@array);
sub union{
    use strict;
    my(@hoge) = @_;
    my%seen;
    @hoge = grep (!$seen{$_}++, @hoge);
    return @hoge;
}


# Extract $int th column from the 2D array
# @res = &column(\@array, \$int);
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


# Extract $int th characters of str array
# @res = &columnSeq(\@array, \$int);
sub columnSeq {
    use strict;
    my($hoge, $foo) = @_;
    my@out = map{substr($_, $$foo, 1)}@$hoge;
    return @out;
}


# Return the first index of a str in the array.
# Note that this sub could not be used for numbers.
# $res = &position(\@array, \$str);
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
# $res = &log10($float);
sub log10{
    use strict;
    my($x) = @_;
    return log($x) / log (10);
}


# Calculate log2
# $res = &log2($float);
sub log2{
    use strict;
    my($x) = @_ ;
    return log($x) / log (2);
}


# Calculate median
# $res = &median(@array);
sub median{
    use strict;
    return unless @_;
    return $_[0] unless @_ > 1;
    @_= sort{$a<=>$b}@_;
    return$_[$#_/2] if @_&1;
    my$mid= @_/2;
    return($_[$mid-1]+$_[$mid])/2;
}


# Calculate the factorial of arg
# $res = &factorial($int);
sub factorial{
    use strict;
    my($fac) = @_;
    my$res_fac = 1;
    if($fac == 0){
        return 1;
    }
    else{
        for(my$i=1; $i <= $fac; $i++){
            $res_fac *= $i;
        }
        return $res_fac;
    }
}


1;
