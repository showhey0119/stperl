#########################
# Module for CSV format #
#########################

# add "use csv;" to your perl script

# List of functions
# add_hash_csv


# Add CSV info to the array
# &add_hash(\@header, \$line, \%df);
sub add_hash_csv{
    my($header, $line, $df) = @_;
    my@tmp = split(/\,/, $$line);
    my$n_t = scalar@tmp;
    if(scalar@$header != $n_t){
        print"Error in add_hash_csv\n";
        print"No. columns are different between header and a line\n";
        exit;
    }
    %$df = ();
    for(my$i=0; $i < $n_t; $i++){
        $$df{$$header[$i]} = $tmp[$i];
    }
}


1;
