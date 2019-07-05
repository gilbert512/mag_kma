#   Copyright (c) 2007 Boulder Real Time Technologies, Inc.           
#                                                                     
#   This software module is wholly owned by Boulder Real Time         
#   Technologies, Inc. This software may be used freely in any 
#   way as long as the copyright statement above is not removed.
#
#
#   Modified by Sheen (2018.10.10)
#
#   Compute local magnitude from vertical component based on Sheen et al. (2018) 
#   Apply Station corrections, giveny by Sheen et al. (2018)

package Mlkma ;

use lib "$ENV{ANTELOPE}/data/evproc" ;

our @ISA= ( "Magnitude" ) ;

use evproc ;

use strict ;
use warnings ;

use lib "$ENV{ANTELOPE}/data/perl" ;

use Datascope ; 


# SHEEN station correction
# a station with large amplification should be assigned with a negative value
 our %stncrr = (
   'AMD',   -0.0353,
   'BAR',   -0.0595,
   'BGD',   -0.1176,
   'BOSB',   0.0310,
   'BRD',   -0.1225,
   'BRD',   -0.1700,
   'BURB',  -0.0005,
   'BUS',    0.0075,
   'BUS2',   0.0239,
   'BUYB',   0.0613,
   'CHC',   -0.0460,
   'CHC2',   0.0276,
   'CHJ',    0.0551,
   'CHJ2',   0.0549,
   'CHNB',   0.0595,
   'CHS',   -0.0351,
   'CHYB',   0.1038,
   'CIGB',   0.1491,
   'DACB',   0.1146,
   'DAG',   -0.1285,
   'DAG2',  -0.1161,
   'DGY',    0.0474,
   'DGY2',   0.0676,
   'DUSB',  -0.0053,
   'ECDB',   0.0963,
   'EMSB',  -0.0473,
   'EURB',   0.1097,
   'EUSB',   0.0578,
   'GAHB',  -0.0787,
   'GAPB',   0.0567,
   'GKP1',  -0.1264,
   'GOCB',   0.0416,
   'GSU',    0.0012,
   'GWYB',   0.0421,
   'HALB',  -0.1267,
   'HAMB',   0.0545,
   'HANB',  -0.0640,
   'HAWB',   0.1190,
   'HKU',   -0.0891,
   'HSB',   -0.0556,
   'HWCB',   0.0250,
   'HWSB',  -0.1447,
   'IMWB',   0.1314,
   'JECB',   0.0917,
   'JEO',   -0.0196,
   'JEO2',  -0.0233,
   'JJB',   -0.1511,
   'JJU',   -0.2899,
   'JRB',   -0.0336,
   'JSB',   -0.0270,
   'KNUD',   0.0296,
   'KOHB',  -0.1311,
   'KRN',    0.0136,
   'KSA',    0.0680,
   'KWJ',   -0.0160,
   'KWJ2',   0.0337,
   'MGB',    0.0087,
   'MKL',    0.0356,
   'MND',   -0.0189,
   'NAWB',   0.0765,
   'NPR',    0.0392,
   'OKCB',  -0.0043,
   'OKEB',  -0.0145,
   'OYDB',   0.0936,
   'SEHB',  -0.0036,
   'SEO',    0.0394,
   'SEO2',  -0.0309,
   'SES',    0.0102,
   'SES2',   0.0035,
   'SHHB',   0.2612,
   'SMKB',   0.0062,
   'SND',   -0.0512,
   'SNU',   -0.0388,
   'TJN',   -0.0169,
   'UCN',   -0.0313,
   'ULJ',   -0.0595,
   'ULJ2',  -0.0269,
   'WSN',   -0.1500,
   'YGN',   -0.0544,
   'YKB',    0.0733,
   'YKDB',   0.1133,
   'YNCB',   0.0629,
   'YNDB',   0.0260,
   'YOCB',  -0.0502,
   'YODB',   0.1344,
   'YOJB',   0.1166,
   'YPDB',   0.0458,
   'YSAB',   0.0559,
   'YSB',   -0.1655,
);



sub compml {
	my $self = shift ;
	my $sta = shift ;
	my $millimeters = shift ;

	my $distance = $self->{stations}{$sta}{delta}*111.11 ;

	if ($distance < 0.0 || $distance > 600.0) {return;}
	if ( $millimeters <= 0.0 ) {return;}

#   SHEEN Vertical component (Sheen et al., 2018)
	my $ml = log($millimeters)/log(10) + 0.5107*log($distance/100.)/log(10) + 0.001699*($distance-100.)+3.;

#   SHEEN Station correction;
	if ( defined $stncrr{$sta} ) {
		$ml += $stncrr{$sta};
        }


	return $ml ;
}

sub new {
	return Magnitude::new @_ ;
}

sub getwftimes {
	my $self = shift ;

	my $ret = setup_processes $self ;

	if ($ret ne "ok" ) { return makereturn ( $self, $ret ) ; }

	$self->{stations} = {} ;

	my ($otime,$odepth,$oauth) = dbgetv ( @{$self->{dbo}}, "time", "depth", "auth" ) ;
	my $date = yearday ( $otime ) ;

	if ( defined $self->{params}{auth_accept} ) {
		my $ok = dbex_eval ( @{$self->{dbo}}, "auth =~ /$self->{params}{auth_accept}/" ) ;
		if ( ! $ok ) {
			addlog ( $self, 1, "wrong origin auth " . $oauth ) ;
			return makereturn ( $self, "skip" ) ; 
		}
	}

	my $event_tend = -1.e20 ;
	for ($self->{dba}[3] = 0; $self->{dba}[3] < $self->{nassoc}; $self->{dba}[3]++) {
		my ($sta, $delta) = dbgetv ( @{$self->{dba}} , "sta", "delta" ) ;
		if ( defined $self->{stations}{$sta} ) { next ; }

		my $process ;
		my $channels = {};
		my $ndbv ;

		($ret, $process, $channels, $ndbv) = match_sta ($self, $sta, $otime) ;
		if ( $ret ne "ok" ) { next; }

		if ($delta*111.1 > 600.0) {
			addlog ( $self, 1, $sta . ": station too far away" ) ;
			next ;
		}

		my $pt = dbex_eval ( @{$self->{dbo}}, "ptime(" . $delta . "," . $odepth . ")" ) ;
		my $st = dbex_eval ( @{$self->{dbo}}, "stime(" . $delta . "," . $odepth . ")" ) ;

		my $twin = $process->{signal_twin} ;
		if ( substr($process->{signal_twin}, 0, 1) eq "f") {
			my $fac = substr($process->{signal_twin}, 1) ;
# SHEEN signal window for S or Lg wave
#			$twin = 1.1 * $fac * ($st - $pt) ;
			$twin = $fac * ($st - $pt) ;
		}

		my $noise_twin = $process->{noise_twin};
		if ($process->{noise_twin} eq "tproc") {
			$noise_twin = $twin ;
			if ($noise_twin > 60.0) {$noise_twin = 60.0 ;}
		}

		my $noise_tstart = $otime + $pt - $noise_twin - $process->{noise_toffset} ;
		my $noise_tend = $noise_tstart + $noise_twin ;
# SHEEN signal window for S or Lg wave
		my $signal_tstart = $otime + $st  - $process->{signal_toffset} ;
		my $signal_tend = $signal_tstart + $twin ;

		my $tstart = $noise_tstart - 100.0 ;
		my $tend = $signal_tend ;

		my $hash = {
			"chan_expr" => $process->{chan_expr},
			"delta" => $delta,
			"tstart" => $tstart,
			"tend"	=> $tend,
			"noise_tstart" => $noise_tstart,
			"noise_tend"	=> $noise_tend,
			"signal_tstart" => $signal_tstart,
			"signal_tend"	=> $signal_tend,
			"noise_twin" => $noise_twin,
			"snr_thresh" => $process->{snr_thresh},
			"tupdate" => $self->{params}{update_time},
			"nchans" => $ndbv,
			"channels" => $channels,
			"disposition" => "DataNotReady",
		} ;
		if ( defined $process->{clip_upper} && defined $process->{clip_lower} ) {
			$hash->{clip_upper} = $process->{clip_upper} ;
			$hash->{clip_lower} = $process->{clip_lower} ;
		}
		if ( defined $process->{filter} ) {
			$hash->{filter} = $process->{filter} ;
			if ($hash->{filter} eq "auto") {
				my $expr = sprintf 
					'sta == "%s" && chan =~ /%s/ && %.3f >= time && ( %.3f <= endtime || endtime == null("endtime") )',
						$sta, $process->{chan_expr}, $otime, $otime ;
				my @dbv = dbsubset ( @{$self->{dbc}}, $expr ) ;
				my $ndbv = dbquery ( @dbv, "dbRECORD_COUNT" ) ;
		
				if ($ndbv < 1) {
					addlog ( $self, 0, "station ". $sta . ": no channel matches to "
								. $process->{chan_expr} . " in calibration table" ) ;
					undef $hash->{filter} ;
				} else {
					$dbv[3] = 0;
					my $segtype = dbgetv (@dbv, "segtype");
					if ($segtype eq "V") {
						$hash->{filter} = "WAV" ;
					} elsif ($segtype eq "A") {
						$hash->{filter} = "WAA" ;
					} else {
						addlog ( $self, 0, "station ". $sta . 
							" Cannot determine auto filter for segtype " . $segtype ) ;
						undef $hash->{filter} ;
					}
				}
				dbfree @dbv ;
			} elsif ($hash->{filter} eq "autosp") {
				my $expr = sprintf 
					'sta == "%s" && chan =~ /%s/ && %.3f >= time && ( %.3f <= endtime || endtime == null("endtime") )',
						$sta, $process->{chan_expr}, $otime, $otime ;
				my @dbv = dbsubset ( @{$self->{dbc}}, $expr ) ;
				my $ndbv = dbquery ( @dbv, "dbRECORD_COUNT" ) ;
		
				if ($ndbv < 1) {
					addlog ( $self, 0, "station ". $sta . ": no channel matches to "
								. $process->{chan_expr} . " in calibration table" ) ;
					undef $hash->{filter} ;
				} else {
					$dbv[3] = 0;
					my $segtype = dbgetv (@dbv, "segtype");
					if ($segtype eq "V") {
						$hash->{filter} = 'INT s0.2;G 2080.0 1.e-6' ;
					} elsif ($segtype eq "A") {
						$hash->{filter} = 'INT2 s0.2;G 2080.0 1.e-6' ;
					} else {
						addlog ( $self, 0, "station ". $sta . 
							" Cannot determine auto filter for segtype " . $segtype ) ;
						undef $hash->{filter} ;
					}
				}
				dbfree @dbv ;
			}
		}
		$self->{stations}{$sta} = $hash ;
		if ( $signal_tend > $event_tend ) { $event_tend = $signal_tend; }
	}

#	display $self ;

	if ( scalar ( keys ( %{$self->{stations}} ) ) < 1 ) {
		addlog ( $self, 0, "No channels to process" ) ;
		return makereturn ( $self, "ok" ) ; 
	}

	if ( defined $self->{params}{maximum_bad_fraction} ) {
		$self->{maximum_bad_fraction} = $self->{params}{maximum_bad_fraction} ;
	} else {
		$self->{maximum_bad_fraction} = 0.0;
	}

	if ( defined $self->{params}{maximum_wait_time} ) {
		$self->{expire_time} = $event_tend + $self->{params}{maximum_wait_time} ;
		my $now_time = now() + $self->{params}{maximum_wait_time} ;
		if ( $now_time > $self->{expire_time} ) {
			$self->{expire_time} = $now_time ;
		}
	}

	if ( defined $self->{expire_time} ) {
		return makereturn ( $self, "ok", "stations" => $self->{stations},
				"expire_time" => $self->{expire_time} ) ;
	} else {
		return makereturn ( $self, "ok", "stations" => $self->{stations} ) ;
	}
}

sub process_channel {
	my $self = shift ;
	my $ret = $self->SUPER::process_channel(@_) ;

	if ( $ret->{disposition} ne "channeldone" 
		&& $ret->{disposition} ne "stationdone"
		&& $ret->{disposition} ne "processdone" ) {return $ret;}

	my $sta = $ret->{sta} ;
	my $chan = $ret->{chan} ;

	if ( defined $self->{stations}{$sta}{channels}{$chan}{is_nullcalib} 
			&& $self->{stations}{$sta}{channels}{$chan}{is_nullcalib} ) {
		addlog ( $self, 1, "%s: %s: Channel mag not computed because of null calib",
 						$sta, $chan )  ;
		$self->{stations}{$sta}{disposition} = "NullCalib" ;
		return $ret ;
	}
	if ( defined $self->{stations}{$sta}{channels}{$chan}{is_clipped} 
			&& $self->{stations}{$sta}{channels}{$chan}{is_clipped} ) {
		addlog ( $self, 1, "%s: %s: Channel mag not computed because of clipped data",
 						$sta, $chan )  ;
		$self->{stations}{$sta}{disposition} = "DataClipped" ;
		return $ret ;
	}
	if ( ! defined $self->{stations}{$sta}{channels}{$chan}{signal_amax} ) {
		addlog ( $self, 1, "%s: %s: Channel mag not computed because of no data",
 						$sta, $chan )  ;
		return $ret ;
	}
 	if ( defined $self->{stations}{$sta}{channels}{$chan}{snr} ) {
		if ( $self->{stations}{$sta}{snr_thresh} < 1.0
				|| $self->{stations}{$sta}{channels}{$chan}{snr}
					> $self->{stations}{$sta}{snr_thresh} ) {
			my $millimeters =
 				$self->{stations}{$sta}{channels}{$chan}{signal_amax} ;
			if ( $self->{stations}{$sta}{snr_thresh} >= 1.0 ) {
				$millimeters -= 
 					$self->{stations}{$sta}{channels}{$chan}{noise_std} ;
			} 
 			$self->{stations}{$sta}{channels}{$chan}{m} = compml ( 
				$self, $sta, $millimeters ) ;
 			$self->{stations}{$sta}{channels}{$chan}{m_time} = 
				$self->{stations}{$sta}{channels}{$chan}{signal_tmax} ;
 			$self->{stations}{$sta}{channels}{$chan}{m_snr} = 
				$self->{stations}{$sta}{channels}{$chan}{snr} ;
 			$self->{stations}{$sta}{channels}{$chan}{m_val1} = $millimeters ;
 			$self->{stations}{$sta}{channels}{$chan}{m_units1} = "mmwa" ;
			addlog ( $self, 1, "%s: %s: Channel mag = %.3f",
 					$sta, $chan,
 					$self->{stations}{$sta}{channels}{$chan}{m} ) ;
		} else {
# SHEEN print snr
			addlog ( $self, 1, "%s: %s: Channel mag not computed because of low snr, %6.4f",
 						$sta, $chan, $self->{stations}{$sta}{channels}{$chan}{snr}  )  ;
			$self->{stations}{$sta}{disposition} = "LowSnr" ;
				
		}
	} else {
 		$self->{stations}{$sta}{channels}{$chan}{m} = compml ( $self, $sta, 
 				$self->{stations}{$sta}{channels}{$chan}{signal_amax} ) ;
 		$self->{stations}{$sta}{channels}{$chan}{m_time} = 
				$self->{stations}{$sta}{channels}{$chan}{signal_tmax} ;
 		$self->{stations}{$sta}{channels}{$chan}{m_snr} = 
				$self->{stations}{$sta}{channels}{$chan}{snr} ;
 		$self->{stations}{$sta}{channels}{$chan}{m_val1} = 
 				$self->{stations}{$sta}{channels}{$chan}{signal_amax} ;
 		$self->{stations}{$sta}{channels}{$chan}{m_units1} = "mmwa" ;
		addlog ( $self, 2, "%s: %s: Channel mag = %.3f",
 					$sta, $chan,
 					$self->{stations}{$sta}{channels}{$chan}{m} ) ;
	}

	return $ret ;
}

sub process_network {
	my $self = shift ;
	my $ret = $self->SUPER::process_network(@_) ;

	my @dborigin = @{$self->{dbo}} ;
	$dborigin[3] = 0 ;
	my $auth = dbgetv ( @dborigin, "auth" ) ;
	dbputv ( @dborigin, "auth", $auth . "Ml" ) ;

##########
# SHEEN magnitude from a trimmed mean
#
	if (defined $self->{t_mean} ) {
		dbputv ( @dborigin, "ml", $self->{t_mean}, "mlid", $self->{magid} ) ;
	}

# Original magnitude from a median
#	if (defined $self->{m_median} ) {
#		dbputv ( @dborigin, "ml", $self->{m_median}, "mlid", $self->{magid} ) ;
#	}

	return $ret ;
}

