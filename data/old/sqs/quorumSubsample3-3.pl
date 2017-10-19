# quorumSubsample 3.3 by John Alroy

# first draft 19-20.4.09
# version 1.0 release date 20.6.09
# version 2.0 release date 13.9.09
# new features: Pielou's J correction replaced with exclude-most-common-taxon
#  rule; throwback algorithm added; Gradstein stages scale can be used; TRIM
#  option added; outputs more stats, including occurrences drawn and counts of
#  taxa with one, two, or three occurrences; output file column names revised;
#  prints taxa-by-bins output file 
# version 3.0 release date 12.3.10
# new features: coverage now estimated with a version of Good's u that uses
#  counts of taxa found in a single reference instead of taxa with only one
#  occurrence; abundance data and Gradstein subepochs scale can be used;
#  outputs collection and occurrence totals; prints additional totals to screen
# version 3.1 release date 9.7.10
# new features: bug fixes in computation of quorumSubsampleCurves3.txt file
#  and output of two timer counts
# version 3.2 release date 26.7.10 
# new feature: INTERPOLATE option
# version 3.3 release date 11.12.10
# new feature: DISPERSE option

# location of input and output data files
$PATH = "./";
# options CR, UW, or O2W (default SQ)
$METHOD = "SQ";
# subsample abundance data (only works for CR and SQ; default no)
$ABUND = "";
# if multiple files are listed, subsampled diversity for each one will
#  be reported both jointly and separately
@FILES = ("Ammonoidea-occs.csv","Nautiloidea-occs.csv");
$TRIALS = 100;
# taxonomic rank (option species, default genus)
$RANK = "";
# time scale (options stages, subepochs, Peters, or a number; default
#  is 11 Myr bins)
$SCALE = "";
# a number in Ma is required if SCALE is numerical
$MAXMA = "";
# assign collections spanning multiple bins to the ones including more than
#  half of their age estimate limits (default no; options yes or a fraction)
$DEORPHAN = "";
if ( $DEORPHAN > 0 && $DEORPHAN < 0.5 || $DEORPHAN > 1 )	{
	print "\nExiting because DEORPHAN must be 'yes' or a number between 0.5 and 1.0\n\n";
	exit;
}
if ( $DEORPHAN =~ /y/i )	{
	$DEORPHAN = 0.5;
}
# use interpolated boundary estimates in PaleoDB download file (default no,
#  but recommended)
$INTERPOLATE = "yes";
if ( $SCALE > 0 && $MAXMA == 0 )	{
	print "\nExiting because MAXMA has not been set\n\n";
	exit;
}
# print data for bins sometimes under quota (default no)
$USEFAILED = "y";
# print data only for bins in the range of those that meet the subsampling
#  target (default no)
$TRIM = "";

if ( $METHOD eq "SQ" )	{
	$QUOTA = 0.5;
} elsif ( $METHOD eq "CR" )	{
	$QUOTA = 100;
} elsif ( $METHOD eq "UW" )	{
	$QUOTA = 1000;
} elsif ( $METHOD eq "O2W" )	{
	$QUOTA = 1000;
}

# reference quota (default none, and not recommended)
$REFQUOTA = "";
# disperse sampling among references with a throwback algorithm (optional but recommended)
$DISPERSE = "y";
# number of collections a bin must include to be analyzed
$MINCOLLS = 1;
# intervals must include as many collections as the maximum drawn
#  (default no, and not recommended)
$MATCHMAX = "";


my ($nbin,%lookup,@binlist,@newbinlist);
if ( $SCALE !~ /stages|epochs|peters/i && $SCALE == 0 )	{
	%BINS = ("Cambrian" => 4, "Ordovician" => 5, "Silurian" => 2, "Devonian" => 5, "Carboniferous" => 5, "Permian" => 4, "Triassic" => 4, "Jurassic" => 6, "Cretaceous" => 8, "Cenozoic" => 6);

	for my $p ( "Cambrian", "Ordovician", "Silurian", "Devonian", "Carboniferous", "Permian", "Triassic", "Jurassic", "Cretaceous", "Cenozoic" )	{
		for my $b ( 1..$BINS{$p} )	{
			$nbin++;
			$lookup{$p." ".$b} = $nbin;
			my $i = $p;
			if ( $i =~ /^C/ )	{
				$i =~ s/Cambrian/Cm/;
				$i =~ s/Carboniferous/C/;
				$i =~ s/Cretaceous/K/;
				$i =~ s/Cenozoic/Cz/;
			} else	{
				$i =~ s/(.).*/$1/;
			}
			$binlist[$nbin] = $i.$b;
		}
	}
} elsif ( $SCALE > 0 )	{
	my $ma = 0;
	while ( $ma <= $MAXMA )	{
		$nbin++;
		$ma += $SCALE;
		$binlist[$nbin] = $ma;
	}
	@binlist = sort { $b <=> $a } @binlist;
} elsif ( $SCALE =~ /stages|epochs|peters/ )	{
	if ( $SCALE =~ /stages/ )	{
		open IN,"<./Gradstein_stages";
	} elsif ( $SCALE =~ /epochs/ )	{
		open IN,"<./Gradstein_subepochs";
	} elsif ( $SCALE =~ /peters/ )	{
		open IN,"<./Peters_stages";
	}
	while (<IN>)	{
		s/\n//;
		$nbin++;
		$lookup{$_} = $nbin;
		$binlist[$nbin] = $_;
	}
	close IN;
	# intervals are in reverse order
	for my $i ( reverse 1..$nbin )	{
		$newbinlist[$nbin - $i + 1] = $binlist[$i];
	}
	@binlist = @newbinlist;
	for my $b ( keys %lookup )	{
		$lookup{$b} = $nbin - $lookup{$b} + 1;
	}
}

my (@shortnames,$fref,%incoll,%allcolls,$totaloccs,%bincolls,%binoccs,%binabund,%bybin,%bybinref,%id,%idname,%seen,%abund,%raw);
for my $file ( @FILES )	{
	my $shortname = $file;
	$shortname =~ s/\.(txt|csv|tab)//;
	push @shortnames , $shortname;
	if ( ! open IN,"<$PATH/$file" )	{
		print "\nExiting because $file couldn't be found\n\n";
		exit;
	}
	$_ = <IN>;
	s/\n//;
	my @f = split /,/,$_;
	# coll no is column 0
	my ($fgen,$fsp,$fabund,$fbin,$fbin2,$fmax,$fmin,$fextant);
	for my $i ( 0..$#f )	{
		if ( $f[$i] eq "occurrences.genus_name" )	{
			$fgen = $i;
		} elsif ( $f[$i] eq "occurrences.species_name" )	{
			$fsp = $i;
		} elsif ( $f[$i] eq "occurrences.abund_value" )	{
			$fabund = $i;
		} elsif ( $f[$i] eq "collections.stage" && $SCALE =~ /stages/i )	{
			$fbin = $i;
		} elsif ( $f[$i] eq "collections.subepoch" && $SCALE =~ /epochs/i )	{
			$fbin = $i;
		} elsif ( $f[$i] eq "collections.epoch" && $SCALE =~ /epochs/i )	{
			$fbin2 = $i;
		} elsif ( $f[$i] eq "Peters.interval" && $SCALE =~ /peters/i )	{
			$fbin = $i;
		} elsif ( $f[$i] eq "collections.10mybin" && $SCALE !~ /stages|epochs|peters/i )	{
			$fbin = $i;
		} elsif ( $f[$i] eq "collections.interpolated_base" && $INTERPOLATE =~ /y/i )	{
			$fmax = $i;
		} elsif ( $f[$i] eq "collections.interpolated_top" && $INTERPOLATE =~ /y/i )	{
			$fmin = $i;
		} elsif ( $f[$i] eq "collections.ma_max" && $INTERPOLATE !~ /y/i  )	{
			$fmax = $i;
		} elsif ( $f[$i] eq "collections.ma_min" && $INTERPOLATE !~ /y/i )	{
			$fmin = $i;
		} elsif ( $f[$i] eq "collections.reference_no" )	{
			$fref = $i;
		} elsif ( $f[$i] eq "extant" )	{
			$fextant = $i;
		}
	}

	if ( $SCALE > 0 && ! $fmax )	{
		print "\nExiting because $file doesn't have a maximum estimate field\n\n";
		exit;
	} elsif ( $SCALE > 0 && ! $fmin )	{
		print "\nExiting because $file doesn't have a minimum age estimate field\n\n";
		exit;
	} elsif ( $SCALE == 0 && ! $fbin )	{
		print "\nExiting because $file doesn't have the right time scale field\n\n";
		exit;
	}

	while (<IN>)	{
		s/\n//;
		my @f = split /,/,$_;
		$f[$fbin] =~ s/"//g;
		my $within = $lookup{$f[$fbin]};
		if ( ! $within && $lookup{$f[$fbin2]} )	{
			$within = $lookup{$f[$fbin2]};
		}
		if ( $SCALE > 0 )	{
			my $range = $f[$fmax] - $f[$fmin];
			for my $b ( 1..$nbin )	{
				my ($upper,$lower) = ($f[$fmin],$f[$fmax]);
				if ( $binlist[$b] > $f[$fmin] )	{
					$upper = $binlist[$b];
				}
				if ( $binlist[$b-1] < $f[$fmax] )	{
					$lower = $binlist[$b-1];
				}
				if ( $binlist[$b] <= $f[$fmin] && $binlist[$b-1] >= $f[$fmax] )	{
					$within = $b;
					last;
				} elsif ( $DEORPHAN >= 0.5 && $lower - $upper > $range * $DEORPHAN )	{
					$within = $b;
					last;
				}
			}
		}
		$allcolls{$f[0]}++;
		$totaloccs++;
		if ( $within > 0 && ( $ABUND !~ /y/i || $f[$abund] > 0 ) )	{
			$incoll{$f[0]}++;
			if ( $incoll{$f[0]} == 1 )	{
				$bincolls{$within}++;
				push @{$bybin{$within}} , $f[0];
				push @{$bybinref{$within}{$f[$fref]}} , $f[0];
			}
			$bin{$f[0]} = $within;
			my $name = $f[$fgen];
			if ( $RANK =~ /^s/i )	{
				# knocks out sp. spp. indet. most informals
				if ( $f[$fsp] =~ /[^a-z]/ )	{
					next;
				}
				$name .= " ".$f[$fsp];
			}
			if ( ! $id{$name} )	{
				$ngen++;
				$id{$name} = $ngen;
				$idname{$ngen} = $name;
				$infile{$ngen} = $shortname;
			}
			push @{$list{$f[0]}} , $id{$name};
			$binoccs{$within}++;
			if ( $ABUND =~ /y/i )	{
				$abund{$f[0]}{$id{$name}} = $f[$fabund];
				$binabund{$within} += $f[$fabund];
			}
			if ( ! $seen{$within}{$name} )	{
				$raw{$within}++;
				$seen{$within}{$name}++;
				if ( ! $first{$name} || $bin{$f[0]} < $first{$name} )	{
					$first{$name} = $bin{$f[0]};
				}
				if ( ! $last{$name} || $bin{$f[0]} > $last{$name} )	{
					$last{$name} = $bin{$f[0]};
				}
				$extant{$name} = $f[$fextant];
			}
		}
	}
	close IN;
}

print "\n";
if ( ! $fref )	{
	print "Warning: because the file doesn't include reference numbers, Good's u will be\n based on taxa found in one collection instead of one reference\n\n";
}

my @x = keys %allcolls;
printf "full data set: %d collections and $totaloccs occurrences\n",$#x + 1;
my @temp = keys %incoll;
my $temp = 0;
$temp += $binoccs{$_} foreach keys %binoccs;
printf "binned data set: %d collections and %d occurrences\n",$#temp+1,$temp;

my (%refs,%refseen);
for my $bin ( reverse 1..$nbin )	{
	my @byrefs = keys %{$bybinref{$bin}};
	$refs{$bin} = $#byrefs + 1;
	$refseen{$_}++ foreach keys %{$bybinref{$bin}};
}

if ( $fref > 0 )	{
	my @temp = keys %refseen;
	printf "binned data derive from %d references\n",$#temp+1;
}

my %rt;
my $totalextant;
for my $g ( keys %first )	{
	if ( $extant{$g} =~ /y/i )	{
		$last{$g} = $lookup{'Cenozoic 6'};
		$totalextant++;
	}
	for my $bin ( $first{$g}..$last{$g} )	{
		$rt{$bin}++;
	}
}
if ( $totalextant > 0 )	{
	print "$totalextant of the $ngen binned genera are extant\n\n";
} else	{
	print "none of the $ngen binned genera are extant\n\n";
}

my (%everseen,%failed,%lastfailed,@trialcolls,%collsused,%collsdrawn,%occsdrawn,%freq1,%freq2,%freq3,%ccurve,%atq,%bpdbybin,%qdrawn);
$|=1;
print "trials: ";
for my $t ( 1..$TRIALS )	{
	print "\rtrials: $t ";
	my %seen = ();
	my %lastused = ();
	my %usedintrial = ();
	for my $bin ( reverse 1..$nbin )	{
		if ( $bincolls{$bin} >= $MINCOLLS && ( ! $failed{$bin} || $USEFAILED =~ /y/i ) )	{
			my @colls = ();
			my %inref;
			if ( $REFQUOTA > 0 && $refs{$bin} > $REFQUOTA )	{
				my @byrefs = keys %{$bybinref{$bin}};
				my @subrefs = ();
				while ( $#subrefs + 1 < $REFQUOTA && $#byrefs > -1 )	{
					my $i = int( rand( $#byrefs + 1 ) );
					push @subrefs , $byrefs[$i];
					splice @byrefs , $i , 1;
				}
				for my $s ( @subrefs )	{
					push @colls , @{$bybinref{$bin}{$s}};
					for my $c ( @{$bybinref{$bin}{$s}} )	{
						$inref{$_}{$s}++ foreach @{$list{$c}};
					}
				}
			} else	{
				@colls = @{$bybin{$bin}};
				for my $r ( keys %{$bybinref{$bin}} )	{
					for my $c ( @{$bybinref{$bin}{$r}} )	{
						$inref{$_}{$r}++ foreach @{$list{$c}};
					}
				}
			}
			@trialcolls = @colls;
			my $nocc = 0;
			my $ncoll = $#colls + 1;
			my %freq = ();
			my ($f1,$f2,$f3) = (0,0,0);
			for my $c ( @colls )	{
				if ( $ABUND =~ /y/i )	{
					$nocc += $abund{$c}{$_} foreach @{$list{$c}};
					$freq{$_} += $abund{$c}{$_} foreach @{$list{$c}};
				} else	{
					$nocc += $#{$list{$c}} + 1;
					$freq{$_}++ foreach @{$list{$c}};
				}
			}
			for my $g ( keys %freq )	{
				if ( $freq{$g} == 1 )	{
					$freq1{$bin}++;
					$f1++;
				} elsif ( $freq{$g} == 2 )	{
					$freq2{$bin}++;
					$f2++;
				} elsif ( $freq{$g} == 3 )	{
					$freq3{$bin}++;
					$f3++;
				}
			}
			my @genera = ();
			my ($cdrawn,$odrawn) = 0;
			# classical rarefaction of occurrences
			if ( $METHOD eq "CR" )	{
				my @occs = ();
				my %drawn = ();
				for my $c ( @colls )	{
					if ( $ABUND =~ /y/i )	{
						for my $g ( @{$list{$c}} )	{
							for my $i ( 1..int($abund{$c}{$g}) )	{
								push @occs , $g;
							}
						}
					} else	{
						push @occs , $_ foreach @{$list{$c}};
					}
				}
				if ( $#occs + 1 >= $QUOTA )	{
					for $q ( 1..$QUOTA )	{
						my $i = int( rand( $#occs + 1 ) );
						my $g = $occs[$i];
						splice @occs , $i , 1;
						$drawn{$g}++;
						$odrawn++;
						my @temp = keys %drawn;
						$ccurve{$bin}{$q} += $#temp + 1;
						$atq{$bin}{$q}++;
					}
				} else	{
					$failed{$bin}++;
					$lastfailed{$bin} = $t;
					%drawn = ();
				}
				push @genera , $_ foreach keys %drawn;
			# unweighted list subsampling
			} elsif ( $METHOD eq "UW" )	{
				my %drawn = ();
				for $q ( 1..$QUOTA )	{
					if ( $#colls == -1 )	{
						$failed{$bin}++;
						$lastfailed{$bin} = $t;
						%drawn = ();
						last;
					}
					my $c = int( rand( $#colls + 1 ) );
					$drawn{$_}++ foreach @{$list{$colls[$c]}};
					$cdrawn++;
					$odrawn += $#{$list{$colls[$c]}};
					splice @colls , $c , 1;
					my @temp = keys %drawn;
					$ccurve{$bin}{$q} += $#temp + 1;
					$atq{$bin}{$q}++;
				}
				push @genera , $_ foreach keys %drawn;
			# occurrences squared-weighted list subsampling
			} elsif ( $METHOD eq "O2W" )	{
				my %drawn = ();
				my $sumo2 = 0;
				while ( $sumo2 < $QUOTA )	{
					if ( $#colls == -1 )	{
						$failed{$bin}++;
						$lastfailed{$bin} = $t;
						%drawn = ();
						last;
					}
					my $c = int( rand( $#colls + 1 ) );
					my $o2 = ( $#{$list{$colls[$c]}} + 1 )**2;
					if ( $o2 + $sumo2 > $QUOTA && $o2 + $sumo2 - $QUOTA > $QUOTA - $sumo2 )	{
						last;
					}
					$drawn{$_}++ foreach @{$list{$colls[$c]}};
					$cdrawn++;
					$odrawn += $#{$list{$colls[$c]}};
					$sumo2 += $o2;
					splice @colls , $c , 1;
				}
				push @genera , $_ foreach keys %drawn;
			# default is SQ
			} else	{
				my %drawn = ();
				my ($maxdrawn,$bpd) = (0,0);
				my ($maxocc,$maxc,$maxg) = (0,0,0,0);
				for my $c ( @colls )	{
					if ( $#{$list{$c}} + 1 > $maxocc )	{
						$maxocc = $#{$list{$c}} + 1;
						$maxc = $c;
					}
				}
				my %inmax = ();
				for my $g ( @{$list{$maxc}} )	{
					$inmax{$g}++;
				}
				delete $freq{''};
				for my $g ( keys %freq )	{
					my @refs = keys %{$inref{$g}};
					if ( $#refs > 0 || $inmax{$g} > 0 )	{
						$maxdrawn += $freq{$g};
						if ( $freq{$g} / $nocc > $bpd )	{
							$bpd = $freq{$g} / $nocc;
							$maxg = $g;
						}
					}
				}
				$maxdrawn -= $freq{$maxg};
				if ( $nocc > $freq{$maxg} )	{
					$maxdrawn /= $nocc - $freq{$maxg} ;
				}
				$bpdbybin{$bin} += $bpd;
				$maxbybin{$bin} += $maxdrawn;

				my $sumfreq = 0;
				my %collsinref = ();
				for my $r ( keys %{$bybinref{$bin}} )	{
					my $n = $#{$bybinref{$bin}{$r}} + 1;
					$collsinref{$_} = $n foreach @{$bybinref{$bin}{$r}};
				}
				while ( $sumfreq < $QUOTA )	{
					if ( $#colls == -1 || $nocc - $freq{$maxg} < 1 )	{
						$failed{$bin}++;
						$lastfailed{$bin} = $t;
						%drawn = ();
						$maxbybin{$bin} -= $maxdrawn;
						$bpdbybin{$bin} -= $bpd;
						last;
					}
					my $lastsumfreq = $sumfreq;
					my %lastdrawn = %drawn;
					my $c = -1;
					if ( $DISPERSE =~ /y/i )	{
						while ( $c == -1 )	{
							$c = int( rand( $#colls + 1 ) );
							if ( rand() > 1 / $collsinref{$colls[$c]} )	{						$c = -1;
							}
						}
					} else	{
						$c = int( rand( $#colls + 1 ) );
					}
					for my $g ( @{$list{$colls[$c]}} )	{
						if ( ! $drawn{$g} && $g ne $maxg )	{
							$sumfreq += $maxdrawn * $freq{$g} / ( $nocc - $freq{$maxg} );
						}
						if ( ! $drawn{$g} )	{
							$drawn{$g}++;
						}
					}
					# throwback algorithm
					if ( $sumfreq - $QUOTA > $QUOTA - $lastsumfreq )	{
						$sumfreq = $lastsumfreq;
						%drawn = %lastdrawn;
						last;
					}
					$cdrawn++;
					if ( $ABUND =~ /y/i )	{
						$odrawn += $abund{$colls[$c]}{$_} foreach @{$list{$colls[$c]}};
					} else	{
						$odrawn += $#{$list{$colls[$c]}} + 1;
					}
					splice @colls , $c , 1;
				}
				if ( $lastfailed{$bin} != $t )	{
					$qdrawn{$bin} += $sumfreq;
				}
				push @genera , $_ foreach keys %drawn;
			}
			if ( $lastfailed{$bin} != $t )	{
				for my $c ( @trialcolls )	{
					$collsused{$bin}++;
					for my $g ( @{$list{$c}} )	{
						if ( $lastused{$g} != $bin )	{
							$genused{$bin}++;
							$usedintrial{$bin}++;
							$lastused{$g} = $bin;
						}
					}
				}
				if ( $METHOD !~ /CR/ )	{
					$collsdrawn{$bin} += $cdrawn;
				}
				$occsdrawn{$bin} += $odrawn;
				for my $g ( @genera )	{
					if ( $g > 0 && ! $seen{$bin}{$g} )	{
						$sib{$bin}++;
						$byfile{$bin}{$infile{$g}}++;
					}
					if ( $g > 0 && ! $seen{$bin}{$g} && $lastfailed{$bin+1} != $t && $lastfailed{$bin+2} != $t )	{
						
						if ( $seen{$bin+1}{$g} )	{
							$two{$bin+1}++;
						}
						if ( $seen{$bin+1}{$g} && $seen{$bin+2}{$g} )	{
							$three{$bin+1}++;
							$sumthree++;
							$binsumthree{$bin+1}++;
						}
						if ( ! $seen{$bin+1}{$g} && $seen{$bin+2}{$g} )	{
							$part{$bin+1}++;
							$sumpart++;
							$binsumpart{$bin+1}++;
						}
					}
					$seen{$bin}{$g}++;
					$everseen{$g}{$bin}++;
				}
			} elsif ( $lastfailed{$bin} == $t && $USEFAILED !~ /y/i )	{
				$sib{$bin} = "";
				$two{$bin+1} = "";
				$three{$bin+1} = "";
				$sumthree -= $binsumthree{$bin+1};
				$part{$bin+1} = "";
				$sumpart -= $binsumpart{$bin+1};
			}
			if ( $lastfailed{$bin} == $t )	{
				my %freq;
				for my $c ( @{$bybin{$bin}} )	{
					$freq{$_}++ foreach @{$list{$c}};
				}
				$everseen{$_}{$bin}++ foreach keys %freq;
			}
		} else	{
			my %freq;
			for my $c ( @{$bybin{$bin}} )	{
				$freq{$_}++ foreach @{$list{$c}};
			}
			$everseen{$_}{$bin}++ foreach keys %freq;
		}
	}
}
print "\n";

my ($firstbin,$lastbin) = ($nbin,1);
if ( $TRIM =~ /y/i )	{
	($firstbin,$lastbin) = ("","");
	for $b ( reverse 1..$nbin )	{
		if ( $bincolls{$b} >= $MINCOLLS && ( $failed{$b} < $TRIALS || $USEFAILED =~ /y/i ) )	{
			if ( ! $firstbin )	{
				$firstbin = $lastbin;
			}
			$lastbin = $b;
		}
	}
	if ( ! $firstbin )	{
		$firstbin = $lastbin;
	}
}

if ( $METHOD =~ /CR|UW/ )	{
	open OUT,">./quorumSubsampleCurves3.txt";
	print OUT "items";
	for $b ( reverse $lastbin..$firstbin )	{
		print OUT "\t$binlist[$b]";
	}
	print OUT "\n";
	for my $q ( 1..$QUOTA )	{
		my ($step,$step2) = (0,1);
		if ( $q > 1 )	{
			$step = int( log( $q - 1 ) / log( 10 ) * 10 );
			$step2 = int( log( $q ) / log( 10 ) * 10 );
		}
		if ( $step < $step2 || $q == $QUOTA )	{
			print OUT "$q";
			for $b ( reverse $lastbin..$firstbin )	{
				if ( $atq{$b}{$q} > 0 )	{
					printf OUT "\t%.1f",$ccurve{$b}{$q}/$atq{$b}{$q};
				} else	{
					print OUT "\tNA";
				}
			}
			print OUT "\n";
		}
	}
	close OUT;
}

my %trials;
for my $bin ( reverse 1..$nbin )	{
	$trials{$bin} = $TRIALS;
	if ( $USEFAILED =~ /y/i )	{
		$trials{$bin} = $TRIALS - $failed{$bin};
	}
}

if ( $MATCHMAX =~ /y/i )	{
	my $maxdrawn;
	for my $bin ( reverse 1..$nbin )	{
		if ( $collsdrawn{$bin} / $trials{$bin} > $maxdrawn )	{
			$maxdrawn = $collsdrawn{$bin} / $trials{$bin};
		}
	}
	for my $bin ( reverse 1..$nbin )	{
		if ( $bincolls{$bin} < $maxdrawn || ( $collsused{$bin} > 0 && $collsused{$bin} / $trials{$bin} < $maxdrawn ) )	{
			$maxbybin{$bin} = "";
			$bpdbybin{$bin} = "";
			$qdrawn{$bin} = "";
			$sib{$bin} = "";
			$two{$bin+1} = "";
			$three{$bin+1} = "";
			$sumthree -= $binsumthree{$bin+1};
			$part{$bin+1} = "";
			$sumpart -= $binsumpart{$bin+1};
			$failed{$bin} = $TRIALS;
		}
	}
}

open OUT,">./quorumSubsample3.txt";
print OUT "Bin name\t";
if ( $SCALE > 0 )	{
	print OUT "Midpoint\t";
}
if ( $fref > 0 )	{
	print OUT "References\t";
}
print OUT "Collections\t";
if ( $REFQUOTA > 0 )	{
	print OUT "Collections used\t";
}
if ( $METHOD !~ /CR/ )	{
	print OUT "Collections drawn\t";
}
print OUT "Occurrences\t";
if ( $ABUND =~ /y/i )	{
	print OUT "Specimens\t";
	print OUT "Specimens drawn\t";
} else	{
	print OUT "Occurrences drawn\t";
}
print OUT "Range through\t";
if ( $RANK =~ /^s/i )	{
	print OUT "Raw species\t";
} else	{
	print OUT "Raw genera\t";
}
if ( $REFQUOTA > 0 && $RANK =~ /^s/i )	{
	print OUT "Species used\t";
} elsif ( $REFQUOTA > 0 )	{
	print OUT "Genera used\t";
}
print OUT "One occurrence\tTwo occurrences\tThree occurrences\tChao-2\t";
print OUT "Subsampled diversity\tThree timer diversity\t";
if ( $#FILES > 0 )	{
	for my $file ( @shortnames )	{
		print OUT "$file\t";
	}
}
print OUT "Two timers\tThree timers\tPart timers\tThree timer sampling";
print OUT "\tCoverage";
print OUT "\tDominance";
print OUT "\tQuorum drawn";
if ( $USEFAILED =~ /y/i )	{
	print OUT "\tFailed trials";
}
print OUT "\n";
my $grandu = 1;
if ( $sumthree > 0 )	{
	$grandu = $sumthree / ( $sumthree + $sumpart );
	printf "overall average three timer sampling stat: %.3f\n",$grandu;
}
my @under;
for $b ( reverse $lastbin..$firstbin )	{
	if ( $SCALE == 0 )	{
		print OUT "\"$binlist[$b]\"\t";
	} else	{
		printf OUT "%d\t",$binlist[$b-1];
		printf OUT "%.2f\t",( $binlist[$b] + $binlist[$b-1] ) / 2;
	}
	if ( ! $fref )	{
		printf OUT "$refs{$b}\t";
	} elsif ( $refs{$b} > 0 )	{
		printf OUT "$refs{$b}\t";
		print OUT "$bincolls{$b}\t";
	} else	{
		print OUT "NA\tNA\t";
	}
	if ( $bincolls{$b} >= $MINCOLLS && $failed{$b} < $TRIALS )	{
		if ( $REFQUOTA > 0 )	{
			$collsused{$b} /= $trials{$b};
			printf OUT "%.2f\t",$collsused{$b};
		}
		if ( $METHOD !~ /CR/ )	{
			$collsdrawn{$b} /= $trials{$b};
			printf OUT "%.2f\t",$collsdrawn{$b};
		}
		print OUT "$binoccs{$b}\t";
		if ( $ABUND =~ /y/i )	{
			print OUT "$binabund{$b}\t";
		}
		$occsdrawn{$b} /= $trials{$b};
		printf OUT "%.2f\t",$occsdrawn{$b};
		$genused{$b} /= $trials{$b};
		$freq1{$b} /= $trials{$b};
		$freq2{$b} /= $trials{$b};
		$freq3{$b} /= $trials{$b};
		$sib{$b} /= $trials{$b};
		$two{$b} /= $trials{$b};
		$three{$b} /= $trials{$b};
		$part{$b} /= $trials{$b};
		print OUT "$rt{$b}\t";
		print OUT "$raw{$b}\t";
		if ( $REFQUOTA > 0 )	{
			printf OUT "%.2f\t",$genused{$b};
			printf OUT "%.2f\t%.2f\t%.2f\t",$freq1{$b},$freq2{$b},$freq3{$b};
			if ( $freq2{b} > 0 )	{
				printf OUT "%.2f\t",$genused{$b} + ( $freq1{$b}**2 / ( 2 * $freq2{$b} ) );
			} else	{
				print OUT "NA\t";
			}
		} else	{
			printf OUT "%.2f\t%.2f\t%.2f\t",$freq1{$b},$freq2{$b},$freq3{$b};
			if ( $freq2{$b} > 0 )	{
				printf OUT "%.2f\t",$raw{$b} + ( $freq1{$b}**2 / ( 2 * $freq2{$b} ) );
			} else	{
				print OUT "NA\t";
			}
		}
		if ( $sib{$b} > 0 )	{
			printf OUT "%.2f\t",$sib{$b};
		} else	{
			print OUT "NA\t";
		}
		if ( $three{$b} > 0 && $sib{$b} > 0 )	{
			my $u = 1;
			if ( $three{$b} > 0 )	{
				$u = $three{$b} / ( $three{$b} + $part{$b} );
			}
		# SIB corrected with sampling stat
			printf OUT "%.2f\t",$sib{$b} * $grandu / $u;
			if ( $#FILES > 0 )	{
				for my $file ( @shortnames )	{
					$byfile{$b}{$file} /= $trials{$b};
					printf OUT "%.2f\t",$byfile{$b}{$file} * $grandu / $u;
				}
			}
			printf OUT "%.2f\t",$two{$b};
			printf OUT "%.2f\t",$three{$b};
			printf OUT "%.2f\t",$part{$b};
			printf OUT "%.3f",$u;
		} else	{
			print OUT "NA\t";
			if ( $#FILES > 0 )	{
				for my $file ( @shortnames )	{
					if ( $sib{$b} > 0 )	{
						$byfile{$b}{$file} /= $trials{$b};
						printf OUT "%.2f\t",$byfile{$b}{$file};
					} else	{
						print OUT "NA\t";
					}
				}
			}
			print OUT "NA\tNA\tNA\tNA";
		}
		$maxbybin{$b} /= $trials{$b};
		printf OUT "\t%.3f",$maxbybin{$b};
		$bpdbybin{$b} /= $trials{$b};
		$qdrawn{$b} /= $trials{$b};
		printf OUT "\t%.4f",$bpdbybin{$b};
		printf OUT "\t%.3f",$qdrawn{$b};
		if ( $USEFAILED =~ /y/i )	{
			printf OUT "\t%d",$failed{$b};
		}
		print OUT "\n";
	} else	{
		if ( $SCALE == 0 )	{
			push @under , $binlist[$b];
		} else	{
			push @under , $binlist[$b-1];
		}
		if ( $REFQUOTA > 0 )	{
			print OUT "NA\t";
		}
		if ( $METHOD !~ /CR/ )	{
			print OUT "NA\t";
		}
		if ( $binoccs{$b} > 0 )	{
			print OUT "$binoccs{$b}\t";
		} else	{
			print OUT "NA\t";
		}
		if ( $ABUND =~ /y/i && $binabund{$b} > 0 )	{
			print OUT "$binabund{$b}\t";
		} elsif ( $ABUND =~ /y/i )	{
			print OUT "NA\t";
		}
		print OUT "NA\t";
		if ( $raw{$b} > 0 )	{
			print OUT "$rt{$b}\t$raw{$b}\t";
		} else	{
			print OUT "NA\tNA\t";
		}
		# genera used
		if ( $REFQUOTA > 0 )	{
			print OUT "NA\t";
		}
		# f1, f2, f3, ICE
		print OUT "NA\tNA\tNA\tNA";
		# SIB and u
		print OUT "\tNA\tNA";
		if ( $#FILES > 0 )	{
			for my $file ( @shortnames )	{
				print OUT "\tNA";
			}
		}
		# 2T, 3T, PT, sampling
		print OUT "\tNA\tNA\tNA\tNA";
		# coverage, dominance, quorum drawn
		print OUT "\tNA\tNA\tNA";
		if ( $USEFAILED =~ /y/i )	{
			print OUT "\tNA";
		}
		print OUT "\n";
	}
}
close OUT;

open OUT,">./quorumSubsampleRanges3.txt";
print OUT "genus";
for $b ( $lastbin..$firstbin )	{
	print OUT "\t$binlist[$b]";
}
print OUT "\n";
my @genera = keys %everseen;
@genera = sort { $idname{$a} cmp $idname{$b} } @genera;
for my $g ( @genera )	{
	print OUT $idname{$g};
	for $b ( $lastbin..$firstbin )	{
		if ( $everseen{$g}{$b} == 0 )	{
			$n = "0";
		} elsif ( $TRIALS <= 10 )	{
			$n = sprintf "%.1f",$everseen{$g}{$b} / $TRIALS;
		} else	{
			$n = sprintf "%.2f",$everseen{$g}{$b} / $TRIALS;
		}
		print OUT "\t$n";
	}
	print OUT "\n";
}
close OUT;

if ( $#under == 0 )	{
	print "$under[0] is under quota\n";
} elsif ( $#under == 1 )	{
	print "$under[0] and $under[1] are under quota\n";
} elsif ( $#under > 0 )	{
	$under[$#under] = "and ".$under[$#under];
	printf "%d intervals are under quota (".join(', ',@under).")\n",$#under + 1;
}

print "\n";

