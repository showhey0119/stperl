###########################
# Module for vcf analysis #
###########################

# add "use vcf;" to your perl script

# List of functions

# Remove unnecessary info from the array
# &shift9(\@arr);
sub shift9{
    my($arr) = @_;
    for(my$i=0; $i < 9; $i++){
        shift@$arr;
    }
}


# get genotypes as the array
# @res = &get_genotype(\$line);
sub get_genotypes{
    my($line) = @_;
    my@tmp = split(/\t/, $$line);
    &shift9(\@tmp);
    return@tmp;
}


# Convert the array of alleles from vcf genotype
# @allele = &get_allele(\$genotype);
sub get_allele{
    my($gen) = @_;
    if($$gen =~ /^([01\.])[\/\|]([01\.])/){
        my$a1 = $1;
        my$a2 = $2;
        if($a1 eq "."){
            if($a2 ne "."){
                print"E1 in get_allele\n$$gen\n";
                exit;
            }
            return(-9, -9);
        }
        else{
            return($a1, $a2);
        }
    }
    else{
        print"E2 in get_allele\n$$gen\n";
        exit;
    }
}


# Calculate proportion of missing data
# $p = &prop_missing(\@genotype_array);
sub prop_missing{
    my($gen) = @_;
    my$n_gen = scalar@$gen;
    my$n_miss = 0;
    foreach my$g(@$gen){
        my@allele = &get_allele(\$g);
        if($allele[0] == -9){
            $n_miss++;
        }
    }
    $n_miss /= $n_gen;
    return$n_miss;
}


# Set hash of VCF info
# &hash_vcf_info(\%hash, \@header, \$vcf);
sub hash_vcf_info{
    my($hs, $header, $line) = @_;
    my@vcf = split(/\t/, $$line);
    my$n_vcf = scalar@vcf;
    if($n_vcf != scalar@$header){
        print"E1 in hash_vcf_info\n";
        print"@$header\n$$line\n";
        exit;
    }
    %$hs = ();
    for(my$i=0; $i < $n_vcf; $i++){
        $$hs{$$header[$i]} = $vcf[$i];
    }
}





1;
