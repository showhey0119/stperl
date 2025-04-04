###########################
# Module for vcf analysis #
###########################

# add "use vcf;" to your perl script

# List of functions
# get_format
# check_missing
# get_header_acc
# shift9
# get_genotypes
# get_allele
# prop_missing
# hash_vcf_info


# Fileter a genothpe
# &filter_genotype(\$gen);
sub filter_genotype{
    my($gen) = @_;
    if($$gen =~ /^[01][\/\|][01]/){
        $$gen =~ s/^[01][\/\|][01]/\.\/\./;
    }
    else{
        print"E1 in filter_genotype\n";
        print"Unknown format: $$gen\n";
        exit;
    }
}



# Get the FORMAT field of a genotype as a hash
# &get_format(\@format, \$gen, \%res);
sub get_format{
    my($format, $gen, $res) = @_;
    my@tmp = split(/\:/, $$gen);
    if(scalar@$format != scalar@tmp){
        print"E1 in get_format: # columns are different\n";
        print"@$format\n@tmp\n"; exit;
    }
    %$res = ();
    for(my$i=0; $i < scalar@$format; $i++){
        if($tmp[$i] =~ /^\.$/){
            $tmp[$i] = 0;
        }
        $$res{$$format[$i]} = $tmp[$i];
    }
 }


# Check if a genotype is missing
# &check_missing(\$gen);
sub check_missing{
    my($gen) = @_;
    if($$gen =~ /^\.\/\./){
        return 1;
    }
    elsif($$gen =~ /^[01][\/\|][01]/){
        return 0;
    }
    else{
        print"E1 in check_missing: unknown format\n";
        print"$$gen\n";
        exit;
    }
}


# Get header and sample names from #CHROM
# &get_header_acc(\@header, \@acc, \$n_acc, \$of);
sub get_header_acc{
    my($header, $acc, $n_acc, $line) = @_;
    @$header = split(/\t/, $$line);
    @$acc = @$header;
    &shift9(\@$acc);
    $$n_acc = scalar@$acc;
}


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
