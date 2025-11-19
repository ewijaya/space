#!/usr/bin/perl -w
#===============================================================================
#
#
#===============================================================================
use strict;
use Data::Dumper;
use Carp;
use List::Util qw(max min sum);

my $filename = $ARGV[0]; 
my $in_lim   = 0.6;
my $top      = 10;
my $sbl_llim = 2;
my $sbl_ulim = 5000;

my %bhash = parse_to_aoa($filename);
my %hash = conv_to_ranked_hash(\%bhash);

# Only take Top 10 for finding best motif
while ( my ($k,$v) = each %hash ) {
    my ($rank) = (split (" - ", $k))[0]; 
    (my $rankno) = $rank =~ /(\d+) Rank/;
 
    if ( $rankno > $top || $v->[0] <= 0  ) {
        delete $hash{$k};
    }
    
}               

# Sort by Score
# Tie break use Rank in Param
my @motif_set = do {no warnings q/numeric/; 
                     sort { 
                     $hash{$b}->[0] <=> $hash{$a}->[0]
                                ||
                             $a <=> $b 
            }  (keys %hash)};


# Find siblings - Must be at least 1
# When both strings are of same or different length
# they must have at least 60% matching position

my %sibling_set;

foreach my $ms (@motif_set) {
    my ( $rank, $param, $motif ) = split( " - ", $ms );
    #print "$rank $param $motif\n";

    my @sbls = find_sibling( $motif, \@motif_set );
    push @{ $sibling_set{$ms} }, @sbls;
}






# Printing them
my $cnt;
foreach my $mset (@motif_set) {
    my $no_of_sbl = scalar( @{ $sibling_set{$mset} } );

    my ($mot_e_len) = ( split( " - ", $mset ) )[-1];
    my ($param)     = ( split( " - ", $mset ) )[1];
    my ($mot)       = ( split( " ",   $mot_e_len ) )[0];

    if ( $no_of_sbl >= $sbl_llim && $no_of_sbl < $sbl_ulim ) {
        $cnt++;

        if ( $cnt == 1 ) {
            print  "==== BEST Motif Seems To Be: === \n";

            my ( $mot, $e, $len ) = split( " ", $mot_e_len );


            print "MOTIF: ", change_motif_to_iupac($mot), "\n\n";

            my @red_mot;
            foreach my $sb ( @{ $sibling_set{$mset} } ) {

                my $motif = ( split( " - ", $sb ) )[1];
                my ( $m, $bm, $l ) = split( " ", $motif );
                push @red_mot, $m;
            }

            my @uredmot = uniq (@red_mot);
            my @iured_mot = map { change_motif_to_iupac($_) } @uredmot;

            my $wscr = @{ $hash{$mset} }[0];
            print "Significance Score : $wscr \n";

            my @inset = @{ $hash{$mset} }[ 1 .. $#{ $hash{$mset} } ];

            print join( "\n", @inset ), "\n";
            print "\n\n";

        }
        else {
            print "==== Other Potentially Good Motif: === \n";
            my ( $mot, $e, $len ) = split( " ", $mot_e_len );

            print  "MOTIF: ", change_motif_to_iupac($mot), "\n\n";


            my @red_mot;
            foreach my $sb ( @{ $sibling_set{$mset} } ) {

                my $motif = ( split( " - ", $sb ) )[1];
                my ( $m, $bm, $l ) = split( " ", $motif );
                push @red_mot, $m;
            }

            my @uredmot = uniq (@red_mot);
            my @iured_mot = map { change_motif_to_iupac($_) } @uredmot;

            print "Significance Score : $hash{$mset}->[0] \n";
            print join( "\n", @{ $hash{$mset} }[ 1 .. $#{ $hash{$mset} } ] ),"\n";
            print "\n";
        }
    }

}



#--------------------------------------------------
# Subroutines 
#-------------------------------------------------- 
sub srt_sel_uniq {

	my	@arr	= @_;
    
    my @res = map {[split(" -- ",$_)]} @arr;

    my @store;

    INS_LIST: 
    foreach my $res (
        do {
            no warnings qw /numeric/;
            sort { $a->[0] <=> $b->[0] || $b->[-1] <=> $a->[-1] }
                @res;
        }
        )
    {
        my ( $sqid, $pos, $str ) = split( ",", $res->[0] );
        print "$res->[0] -- $res->[-1]\n";
        push @store, $sqid;

    }
	
	return ;
}

sub check_appear_in_array {

    # Check if an elem appear in
    # an array
    my ( $elem, $ar ) = @_;
    my %elements;

    # Populate hash
    foreach ( @{$ar} ) {
        $elements{$_} = 1;
    }

    # Check them
    exists $elements{$elem} ? return 1: return 0;

}

sub find_sibling {

	my	($mt_we_len,$list)	= @_;
    my  ($mt, $we, $len) = split (" ",$mt_we_len);

    my @siblings;
    
    LS:
    foreach my $ls ( @{$list} ) {

        my ( $ls_rnk, $ls_param, $ls_mtf ) = split( " - ", $ls );
        my ( $lmt,    $lwe,      $llen )   = split( " ",   $ls_mtf );

        next LS if ( $ls_mtf eq $mt_we_len );
        if ( is_intersection( $mt, $lmt, $in_lim ) > 0 ) {
            push @siblings, $ls_param." - ".$ls_mtf;
        }
    }

	return @siblings;
}

sub find_sibling_with_rank {

	my	($mt_we_len,$list)	= @_;
    my  ($mt, $we, $len) = split (" ",$mt_we_len);

    my @siblings;
    
    LS:
    foreach my $ls ( @{$list} ) {

        my ( $ls_rnk, $ls_param, $ls_mtf ) = split( " - ", $ls );
        my ( $lmt,    $lwe,      $llen )   = split( " ",   $ls_mtf );

        next LS if ( $ls_mtf eq $mt_we_len );
        #if ( is_intersection( $mt, $lmt, $in_lim ) > 0 ) {
        if ( motif_similarity( $mt, $lmt) > $in_lim ) {
            push @siblings, $ls;
        }
    }

	return @siblings;
}

sub is_intersection {

    my ( $s1, $s2 , $int_lim) = @_;
    my $source;
    my $dest;

    if ( length($s1) <= length($s2) ) {
        $source = $s1;
        $dest   = $s2;
    }
    else {
        $source = $s2;
        $dest   = $s1;
    }

    my $len_source = length($source);
    my @dest_list  = substrings( $dest, $len_source );

    my $lim = ( $int_lim * length($source) ) + 1;

    my $exist_count = 0;
    foreach my $dest_list (@dest_list) {

        my $count = get_match_count( $source, $dest_list );
      
        if ( $count >= $lim ) {
            $exist_count++;
        }
    }

    $exist_count > 0 ? return 1: return 0;
}


sub get_match_count_no_N {

    my ( $s1, $s2 ) = @_;
    $s1 =~ tr/NS/ns/;
    my $merge     = $s1 ^ $s2;
    my $num_match = $merge =~ tr/\00//;

    return $num_match;
}

sub get_match_count {

    # bases compares with N is considered a match

    my ( $source_, $target_ ) = @_;
    my @sparts = ( $source_ =~ /(\[.*?\]|.)/g );

    my $source = $source_;
    $source =~ s/N/\[ATCG\]/g;

    my $target = $target_;
    $target =~ s/N/\[ATCG\]/g;

    my $sc       = bitwise($source);
    my $tg       = bitwise($target);
    my $mcount = matches( $sc, $tg );

    return $mcount;
}



sub matches {
    my ( $source, $target ) = @_;
    my $mmatch = ( $source & $target ) =~ tr/\0//;
    return (length($source) - $mmatch); 
}

sub bitwise {
    my $string = shift;
    my %bits = ( 'A' => 1, 'C' => 2, 'G' => 4, 'T' => 8 );
    join '', map {
        my $char = 0;
        $char |= $bits{$_} for /[ACGT]/g;
        chr($char)
    } $string =~ /(\[.*?\]|.)/g;
}

sub is_substring {

    my ( $source, $target ) = @_;

    # Return -1 if one is not substring
    # of the other, both of them must not be of the 
    # same length

    my $index = 0;
    if ( length $source < length $target ) {
        $index = index $target, $source;
    }
    elsif (length $source > length $target) {
        $index = index $source, $target;
    }

    return $index;
}

sub hd {
    my ( $str1, $str2 ) = @_;
    return ( $str1 ^ $str2 ) =~ tr/\001-\255//;
}

# End Motif intersection measure
sub motif_similarity {

    # dif len

	my	($s1,$s2)= @_;
    my $source;
    my $dest;
    if ( length($s1) <= length($s2) ) {
        $source = $s1;
        $dest   = $s2;
    }
    else {
        $source = $s2;
        $dest   = $s1;
    }

    my $best = -1;
    my $diff = 0;
    my $diff2 = 0;
    my $maxlen = 0;   

   foreach my $bp ( 0 .. ( length($dest) - 1 ) ) {
        my $fp   = $bp;

        IN:
        foreach my $bp2 ( 0 .. ( length($source) - 1 ) ) {
            my $fp2   = $bp2;
            my $len=length($source) - $bp2;
            if ( ( length($dest) - $bp  ) < ( length($source) - $bp2 ) )
            {
                $len = ( length($dest) - $bp  );
            }

            my $sstr = substr($dest,$fp,$len);
            my $sstr2 = substr($source,$fp2,$len);
            my $mcount = get_match_count_iupac( $sstr, $sstr2 );


            if ( $mcount > $best ) {
                $diff  = length($dest) - $len;
                $diff2 = length($source) - $len;
                $maxlen = $diff + $len + $diff2;
                #print "DEST: $bp $sstr  --- SOURCE: $bp2 $sstr2 -- LEN $len MAXLEN = $maxlen MATCH: $mcount\n";
                $best   = $mcount;
            }


        }

    }

    
    my $sim = $best/$maxlen;
    my $simn = sprintf("%.3f",$sim);
	return $simn;
}


sub get_match_count_iupac {

    # bases compares with N is considered a match

    my ( $source_, $target_ ) = @_;
    my @sparts = ( $source_ =~ /(\[.*?\]|.)/g );

    my $source = reverse_iupac($source_);
    my $target = reverse_iupac($target_);
    #print "$source $target\n";

    my $sc       = bitwise($source);
    my $tg       = bitwise($target);
    my $mcount = matches( $sc, $tg );
    #print "$source - $target - $mcount\n";

    return $mcount;
}

sub reverse_iupac {

	my	$amb_str	= shift;
    my @sparts = split(//,$amb_str);
    my %iupac = (
         'M'=> '[AC]',
         'R'=> '[AG]',
         'W'=> '[AT]',
         'S'=> '[CG]',
         'Y'=> '[CT]',
         'K'=> '[GT]',
         'V'=> '[ACG]',
         'H'=> '[ACT]',
         'D'=> '[AGT]',
         'B'=> '[CGT]',
         'N'=> '[ACGT]', 

    );

    my $final_str;
    foreach my $sparts (@sparts) {

        my $nsstr;

        if ( $iupac{$sparts} ) {
            $sparts = $iupac{$sparts};
        }

        $final_str .= $sparts;

    }

	
	return $final_str;
}



sub get_best_inst {

    my $ins_set   = shift;
    my @base_only = map { ( split( ",", $_ ) )[-1] } @{$ins_set};
    my %the_pwm   = compute_pwm(@base_only);


    my @new_set;
    foreach my $in ( @{$ins_set} ) {
        my $ionly = (split (",",$in))[-1];
        my $iscore = get_ins_score( $ionly, \%the_pwm );
            push @new_set, $in . " -- " . $iscore;
    }
    
    return @new_set;
}
             

sub get_ins_score {

    my ( $in, $pwm ) = @_;
    #print Dumper $pwm;

    my @bs = split( //, $in );

    # Find S value from the instance 
   
    my @insp;
    foreach my $p ( 0 .. $#bs ) {
        push @insp, @{ $pwm->{ $bs[$p] } }[$p];
    }
    my $sum_insp = sum(@insp);

    # Find Max Min from each column

    my %col_val;
    while ( my ($k,$v) = each %{$pwm} ) {
        foreach my $i ( 0 .. $#{$v} ) {
            push @{  $col_val{$i} }, $v->[$i];
        }

    }               

         my $max = 0;
         my $min = 0;
         while ( my ($kc,$vc) = each %col_val ) {
           $max += max (@{$vc});
           $min += min (@{$vc});
         }              

    my $ret_val = 0;
    $ret_val = ( ( log($sum_insp) - log( $min + 0.01 ) )
        / ( log($max) - log( $min + 0.01 ) ) ) * 100;
    return $ret_val;
}

sub compute_pwm {
    my @mi = @_;
    my $motif_count;
    my $L;

    #-------Beginning of PWM processing ----------------
    foreach my $mi (@mi) {
        chomp($mi);
        $mi =~ s/\s//g;
        $L = $mi;

        #for motif instances count
        my @words = split( /\W+/, $mi );
        $motif_count += @words;
    }
    $motif_count = 0;

    my $w = length($L);
    my @A = ();
    my @T = ();
    my @C = ();
    my @G = ();

    for ( my $j = 0; $j < $w; $j++ ) {

        # Initialize the base counts.
        my $an = 0;
        my $c = 0;
        my $g = 0;
        my $t = 0;

        foreach my $mi (@mi) {
            chomp($mi);
            my $L = $mi;
            my $sb = substr( $L, $j, 1 );
            while ( $sb =~ /a/ig ) { $an++ }
            while ( $sb =~ /t/ig ) { $t++ }
            while ( $sb =~ /c/ig ) { $c++ }
            while ( $sb =~ /g/ig ) { $g++ }

        }
        push( @A, $an );
        push( @T, $t );
        push( @C, $c );
        push( @G, $g );
    }

    my $sumA = sum(@A);
    my $sumT = sum(@T);
    my $sumC = sum(@C);
    my $sumG = sum(@G);

    my @m    = ();
    my @Anm1 = ();
    my @Tnm1 = ();
    my @Cnm1 = ();
    my @Gnm1 = ();
    my @sPos = ();

    #summing up A,T,C,G for all position
    for ( my $bn = 0; $bn < $w; $bn++ ) {
        my $sumPos = $A[$bn] + $T[$bn] + $C[$bn] + $G[$bn];
        push( @sPos, $sumPos );
        
        # Original version - Indp probabilities
        my $nm1A = $A[$bn]/$sumPos;
        my $nm1T = $T[$bn]/$sumPos;
        my $nm1C = $C[$bn]/$sumPos;
        my $nm1G = $G[$bn]/$sumPos;

        push( @Anm1, $nm1A );
        push( @Tnm1, $nm1T );
        push( @Cnm1, $nm1C );
        push( @Gnm1, $nm1G );
    }


    my %all_hash;
    $all_hash{'A'} = [@Anm1];
    $all_hash{'T'} = [@Tnm1];
    $all_hash{'C'} = [@Cnm1];
    $all_hash{'G'} = [@Gnm1];
    
    return %all_hash;
}

sub conv_to_ranked_hash {

	my	$bh	= shift;
    my %nh;
    my %done;

    OUT:
    foreach my $pg ( sort keys %{$bh} ) {

        my $ar = $bh->{$pg};

        if ( @{$ar} ) {
            #print "$pg\n";
            
            my $count = 0;

            IN:
            foreach my $elem ( @{$ar}  ) {
                $count++;
                
                #print "$count Rank\n";
                #print Dumper $elem ;
                
                my $motif = @{$elem}[0];
                my @rest = @{$elem}[1..$#{$elem}];

                my $dstr = $motif . "-" .$rest[0];

                
                #print "$dstr\n";
                #print Dumper \@rest ;
                 
                next IN if ($done{$dstr});
                push @{$nh{$count." "."Rank"." - ". $pg ." - ". $motif}}, @rest;
                $done{$dstr} = 1;
            }               

        }

    }

	
	return %nh;
}

sub parse_to_aoa {

	my	$fname	= shift;
    my $AoA;
    my $param;
    my $weescore;
    my %bighash;
    my @temp;

    my	$INFILE_file_name = $fname;		# input file name

    open ( INFILE, '<', $INFILE_file_name )
        or croak "$0 : failed to open input file $INFILE_file_name : $!\n";

    while ( <INFILE> ) {
        chomp;

        if ( /^(TAG:) \s (ParamGroup\d+) /xms ) {
            #print "$1 - $2\n";
            $bighash{$2} = $AoA = [];
        }
        elsif (/^MOTIF \s? : \s ([\[\]ATCG]+) \s (.*)/xms ) {
            #print "$1 - $2\n";
            my $motif = $1;
            #my $iumotif = change_motif_to_iupac($motif);
            push @temp, ($motif, $2);
        }
        elsif (/^\d+,-\d+,[ATCGN]+$/) {
            # print "$_\n";
           push @temp, $_;
        }
        elsif (/^=/) {
            push @{$AoA}, [@temp];
            @temp=();
        }

    }               

    close ( INFILE );			# close input file


	return %bighash;
}


sub change_motif_to_iupac {

	my	$amb_str	= shift;
    my @sparts = ( $amb_str =~ /(\[.*?\]|.)/g );


    my %iupac = (

        '[A]'  => 'A',
        '[T]'  => 'T',
        '[C]'  => 'C',
        '[G]'  => 'G',
        '[AC]'  => 'M',
        '[AG]'  => 'R',
        '[AT]'  => 'W',
        '[CG]'  => 'S',
        '[CT]'  => 'Y',
        '[GT]'  => 'K',
        '[ACG]' => 'V',
        '[ACT]' => 'H',
        '[AGT]' => 'D',
        '[CGT]' => 'B',
        '[ACGT]' => 'N',

    );

    my $final_str;
    foreach my $sparts (@sparts) {

        my $nsstr;

        if ( $sparts =~ /\[(.*?)\]/g ) {
            my $sstr = join( "", sort split( //, $1 ) );
            $nsstr = '[' . $sstr . ']';
            $sparts = $iupac{$nsstr};
        }

        $final_str .= $sparts;

    }
	
	return $final_str;

}


# max, min, sum now imported from List::Util (more efficient)

sub uniq {

	my	@array = @_;
    my %count = ();
    my @uniques = grep {! $count{$_}++} @array;
	
	return @uniques;
}

# From String::Substrings;

sub substrings($;$) {
    my ( $string, $length ) = @_;
    return undef unless defined $string;

    # paramter validation
    if ( my $r = ref($string) ) {
        croak "Please call me  `substrings STRING [,LENGTH]`
               but not with substrings $r, LENGTH";
    }
    if (defined($length)
        and ( my $r = ref($length)
            or $length !~ /^\d+/ )
        )
    {
        croak "Please call me substrings STRING [,LENGTH]`
        but not with  substrings '$string', $r";
    }

    my $strlength = length($string);
    my @s         = ();
    if ( defined $length ) {
        return @s
            if $length == 0;
        push @s,
            map { substr $string, $_, $length } ( 0 .. $strlength - $length );
    }
    else {
        foreach my $length ( 1 .. $strlength ) {
            push @s,
                map { substr $string, $_, $length }
                ( 0 .. $strlength - $length );
        }
    }
    return @s;
}
