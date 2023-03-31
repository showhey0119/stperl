###############################################
# Module for changing FORMAT to UNIX and Perl #
###############################################

# qsub
# add_hash_count
# linend
# fileopen
# fileopen2
# fileopenTab2
# TwoDarrayTest2
# fopen_hash_fas3
# fopen_hash_multi2
# ConvertTwoDarray


# Send e-mail
sub qsub{
	my($hoge) = @_;
	system("echo $hoge | mail -s tlp01 xxx\@xxx.jp") ;
}


# Modify hash to cound
sub add_hash_count{
	use strict;
	my($hash, $key) = @_;
	if(exists$$hash{$$key}){
	    $$hash{$$key}++;
	}
	else{
		$$hash{$$key} = 1;
	}
}


# Change line-end
sub linend{ 
    use strict;
	my($directory) = @_; # Directory and file type (e.g. /path//*pl) 
	my@file_name = glob($directory);  # get file names
	foreach my$dir(@file_name){
		open (OUTA, "$dir")||die"$!";
		my@contents = <OUTA>;
		close (OUTA);
		open (OUTB, ">$dir")||die "$!";
		foreach my$tmp(@contents){
			$tmp =~ s/\r//g;
			$tmp =~ s/\n//g;
			print OUTB"$tmp\n";
		}
		close (OUTB);
	}
}


# Open file
sub fileopen{
	use strict;
	my($hoge) = @_;
	open(HOGA, "$hoge")||die"$!";
	my @foo = <HOGA>;
	close(HOGA);
	chomp@foo;
	return @foo;
}


# Open file ver. 2
# Use pass-by-reference 
sub fileopen2{
	use strict;
	my($hoge, $ad) = @_ ;
	open(HOGC, "$$ad")||die"$!";
	@$hoge = <HOGC>;
	close (HOGC);
	chomp@$hoge;
}


# Open file as 2D array, ver. 2
# Use pass-by-reference
sub fileopenTab2{
	use strict;
	my($hoge, $ad) = @_;
	my$of;
	my$d = 0;
    open(HOGA, "$$ad")||die"$!";
	while($of = <HOGA>){
	    chomp$of;
		@{$$hoge[$d]} = split (/\s+/, $of) ;
		$d++ ;
	} # while
	close (HOGA) ;
}


# Test column numbers in 2D array ver.2
sub TwoDarrayTest2{
	use strict;
	my($hoge) = @_;
	my$n_hoge = scalar@$hoge;
	my$c_hoge = scalar@{$$hoge[0]};
	for(my$i=1; $i < $n_hoge; $i++){
		my$ccc = scalar@{$$hoge[$i]};
		if($c_hoge != $ccc){
			print "Error in TwoDarrayTest2 in format\n" ;
			print "Error in $i th row.\n";
			exit;
		}
	}
}


# Open fasta file as hash
sub fopen_hash_fas3{
	use strict;
	my($hoge, $ad) = @_;
	open(DATH, "$$ad")|| die"$!";
	my$of;
	my$id;
	while($of = <DATH>){
	    chomp $of ;
		if($of =~ /^>(.+)$/){
		    $id = $1;
			$id =~ s/\s+//g ;
			if(exists$$hoge{$id}){
                print"Error in fopen_hash_fas3\n";
                print"Multiple $id sequences are found\n";
                exit;
            }
			$$hoge{$id} = "";
		}
		else{
			$id||die"Error in fopen_hash_fas3.\nInput is not fasta\n";
			$of =~ s/\s+//g;
			$$hoge{$id} .= $of;
		}
	}
	close(DATH);
}


# Open tab-column file
sub fopen_hash_multi2{
	use strict;
	my($hoge, $ad) = @_;
	open (DATH, "$$ad")||die"$!";
	my$of;
	while($of = <DATH>){
		chomp $of;
		my@tmp = split(/\s+/, $of);
		my$n_tmp = scalar@tmp;
		my$key = shift@tmp;
		my$val = join("\t", @tmp);
		$$hoge{$key} = $val;
	}
	close (DATH);
}


# Convert array into 2D array
sub ConvertTwoDarray{
    use strict;
	my($moto, $kekka) = @_;
	my$N_moto = scalar@$moto;
	for(my$i=0; $i < $N_moto; $i++){
		@{$$kekka[$i]} = split(/\s+/, $$moto[$i]);
	}
}


1;
