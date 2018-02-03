#!/usr/bin/perl -w
use strict;
use Getopt::Long;

# input
my $system    = $^O;
my $query_sequence			= "";			# the query protein sequence
my $ssite_output_path 	= ""; 		# which is the output file of ssite 
my $jsd_score_path 			= "";			# which is the output file of JSDcalculator
my $bsfpt 							= "";			# which is the LIGAND figure profile library
my $libbsr 							= "";			# which is the BISTE Information library
my $clustering_exe_path	=	"";			# the path of the clustering exe, which you can find in COACH program at ZhangLab

## output
my $finger_dis_dat = "";
my $saveBsites_dat = "";
my $saveBsites_prob_dat = "";
my $saveBsites_clr = "";

##############################################################
### input parameters
unless (lc($system) eq "linux")
{
    printf("Your operating system $system is unsupported at this time\n");
    printf("Currently only Linux is supported\n");
    exit();
}
GetOptions('query_sequence:s' => \$query_sequence, 'ssite_output_path:s' => \$ssite_output_path, 'jsd_score_path:s' => \$jsd_score_path, 'bsfpt:s' => \$bsfpt, 'libbsr:s' => \$libbsr, 'clustering_exe_path:s' => \$clustering_exe_path, 'finger_dis_dat:s' => \$finger_dis_dat, 'saveBsites_dat:s' => \$saveBsites_dat, 'saveBsites_prob_dat:s' => \$saveBsites_prob_dat, 'saveBsites_clr:s' => \$saveBsites_clr);

if(!$query_sequence || !$ssite_output_path || !$jsd_score_path || !$bsfpt || !$libbsr || !$clustering_exe_path || !$finger_dis_dat || !$saveBsites_dat || !$saveBsites_prob_dat || !$saveBsites_clr)
{
    print "\nVotedSelectLBS usage:\nVotedSelectLBS.pl arguments\n";
    print "====================\n";
    print "Mandatory arguments:\n";
    print "====================\n";
    print "-query_sequence query_sequence\n";
    print "-ssite_output_path ssite_output_path\n";
    print "-jsd_score_path jsd_score_path \n";
    print "-bsfpt bsfpt\n";
    print "-libbsr libbsr\n";
		print "-clustering_exe_path clustering_exe_path\n";
		print "-finger_dis_dat save_finger_dis_dat path\n";
		print "-saveBsites_dat saveBsites_dat path\n";
		print "-saveBsites_prob_dat saveBsites_prob_dat path\n";
		print "-saveBsites_clr saveBsites_clr path\n";

    exit();
}

###############################################################
## Initial Global Parameters
my $Lch=length $query_sequence;
my %Qtemplate=();
my $tan_cut        = 1;   #(1-TANI fingerprint) clustering cutoff
my $tan_cut1       = 0.5; #(1-TANI fingerprint) clustering cutoff
my $nTemplates     = 1;           #how many top templates to select
my $nTemplates1    = 10;          #if not enough templates satisfying cutoff, how many top templates to use
my $q_cut          = 0.5;         #overal quality cutoff for selecting templates
my $zcut           = 0.5;         #z-score cut off to be a binding site residue



my

 @AA=qw(A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V );
my @AA1=qw(A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V B Z -);

my %AA2index=
    (
     'A'=>'1', 'R'=>'2', 'N'=>'3', 'D'=>'4', 'C'=>'5', 'Q'=>'6', 'E'=>'7', 'G'=>'8', 'H'=>'9', 'I'=>'10', 'L'=>'11', 'K'=>'12',
     'M'=>'13', 'F'=>'14', 'P'=>'15', 'S'=>'16', 'T'=>'17', 'W'=>'18', 'Y'=>'19', 'V'=>'20', 'B'=>'21', 'Z'=>'22', '-'=>'23');

my %AA2=
    (
     'A'=>'1', 'R'=>'2', 'N'=>'3', 'D'=>'4', 'C'=>'5', 'Q'=>'6', 'E'=>'7', 'G'=>'8', 'H'=>'9', 'I'=>'10', 'L'=>'11', 'K'=>'12',
     'M'=>'13', 'F'=>'14', 'P'=>'15', 'S'=>'16', 'T'=>'17', 'W'=>'18', 'Y'=>'19', 'V'=>'20');
my %SS2num=('C'=>1, 'H'=>2, 'E'=>4);

##################################################################
## load BSITE.bsr file from $libbsr , which is the Binding sites file
my %multi=();
my %bslib=();
my %bsligname=();
open(LIB, "$libbsr");
while(my $line = <LIB>)
{
    chomp($line);
    if(substr($line, 11, 1) eq 'P')
    {
        last;
    }
    if($line=~/(\S+)\s+\:(.+)\:(\S+)\:(\d+)\:(\S+)/)
    {
        my $rec =$1;
				my $ligname=$2;
        my $bsno=$3;
        my $bres=$5;
        my $site=$rec . "_" . $bsno;
        $bslib{$site}=$bres;
				push(@{$multi{$rec}}, $site);
				$ligname =~ s/\s+//g;
				$bsligname{$site}=$ligname;
    }
}
close(LIB);

##################################################################
## load BSITE.ftp file from $bsfpt
my %libfpt=();
open(FPT, "$bsfpt");
while(my $line = <FPT>)
{
    chomp($line);
    if($line=~/(\S+)\s+(\S+)/)
    {
        my $site      = $1;
        my $fpt       = $2;
        $libfpt{$site} = $fpt;
    }
}
close(FPT);


###############################################################
#step 4: select templates based on searching result
###############################################################

my %quality=();
my %pbsr=();
my %feature=();

#####################################################################
### Load the $jsd_score_path file content to %JSD
open(FIN, "$jsd_score_path") or die "No $jsd_score_path\n";
my %JSD=();
while(my $line = <FIN>)
{
    if($line =~ /(\d+)\s+\S+\s+(\S+)/)
    {
			my $res=$1+1;
			my $score=$2;
			$JSD{$res}=$score;
    }
}
close(FIN);

######################################################################
### Select templates from the output file ($ssite_output_path) of ssite cpp program
print "Select templates from the output file ($ssite_output_path) of ssite cpp program\n";

&select_templates("$ssite_output_path");
my @keys = keys %quality;
foreach my$k(@keys)
{
    print "$k $feature{$k} $pbsr{$k}\n", ;
}
my $n_selected=@keys;
my $relieve_flag=0;
if(@keys<$nTemplates)
{
    $relieve_flag=1;
    %quality=();
    print "selected number of templates are too small ($n_selected), use top $nTemplates1\n";  
    &select_templates1("$ssite_output_path");
}
@keys = reverse sort{$quality{$a} <=> $quality{$b}} keys %quality;
$n_selected=@keys;
if($n_selected<1)
{
    die "No templates selected!\n";
}
my @selected=@keys;

my @rst=();
### clustering based on fingerprint similarity
if($relieve_flag==1)
{
    ###compute the pairwise similarity between the fingerprint of the selected template ligand
    my @dis=();
#    open(SIM, ">$$finger_dis_dat");
    open(SIM, ">$finger_dis_dat");
    for(my $i=0; $i<@keys; $i++)
    {            	    
			my $k=$keys[$i];
			#printf SIM "$k\t";
			my @residues=split(/,/, $pbsr{$k});    
			for(my $j=0; $j<@keys; $j++)
			{	
	    	if($i==$j)
	    	{
					$dis[$i][$j]=0;
	    	}
	    	elsif($j<$i)
	    	{
					$dis[$i][$j]=$dis[$j][$i];
	    	}
	    	else
	    	{
					my $k1   = $keys[$j];
					my $fpt1 = $libfpt{$k};
					my $fpt2 = $libfpt{$k1};
		
					$dis[$i][$j] = 1.0 - &tannimoto($fpt1, $fpt2);
	    	}
	    	printf SIM "%.3f ", $dis[$i][$j];
			}
			printf SIM "\n";	
    }
    close(SIM);

    $tan_cut=$tan_cut1;
    @rst=`$clustering_exe_path $finger_dis_dat $n_selected $tan_cut`;    
}
else #no clustering
{
    my $info="cluster 0:";
    for(my $i=0; $i<$n_selected; $i++)
    {
			$info .= "$i ";
    }    
    push(@rst, $info);
}

my %cscore=(); ##total c-score of each cluster
my %clusters=();
my $n_cluster=-1;
my %prob=();
my $bsrncut=0.03*$Lch;
if($bsrncut>10)
{
    $bsrncut=10;
}

&vote(\@rst, "$saveBsites_dat", "$saveBsites_prob_dat", "$saveBsites_clr");
&reRank("$saveBsites_dat", "$saveBsites_prob_dat", "$saveBsites_clr");

print "Done!\n";


##############################################################################
##### Sub function procedure

sub reRank
{
		my ($saveBsites_dat, $saveBsites_prob_dat, $saveBsites_clr)=@_;
	
    my %cscore=();
    my %pred=();
    my $n=0;
    open(IN, "$saveBsites_dat");
    while(my $line=<IN>)
    {
			chomp($line);
			if($line =~ /(\S+)\t/)
			{
	    	$cscore{$n}=$1;
	    	$pred{$n}=$line;
	    	$n++;
			}
    }
    close(IN);

    my @keys = sort{$cscore{$b} <=> $cscore{$a}} keys %cscore;
    my %liginf=();
    my %ligcount=();
    my %clr=();
    open(IN, "$saveBsites_clr") or die "clr is missed\n";
    $n=-1;  
    while(my $line=<IN>)
    {
			#print "$line";
			chomp($line);
	
			if($line =~ /cluster/)
			{
	    	$n++;
				if($n>0)
	    	{
					my @keys = sort {$ligcount{$b}<=>$ligcount{$a}} keys %ligcount;
					my $inf="";		
					foreach my $k(@keys)
					{
		    		$inf .= "$k $ligcount{$k} ";
					}
					my $m=$n-1;
					$inf =~ s/\s+$//;
					$liginf{$m}=$inf;		
	    	}
				%ligcount=();
	    	#print "$n\n";
			}
			else
			{
	    	push(@{$clr{$n}}, $line);
	    	if($line =~ /(\S+)/)
	    	{
					my $site = $1;
					my $lig  = $bsligname{$site};

					if(exists $ligcount{$lig})
					{
		    		$ligcount{$lig}++;		    
					}
					else
					{
		    		$ligcount{$lig}=1;
					}
	    	}
			}
    }
    close(IN);

    if($n>=0)
    {
			my @keys = sort {$ligcount{$b}<=>$ligcount{$a}} keys %ligcount;
			my $inf="";		
			foreach my $k(@keys)
			{
	    	$inf .= "$k $ligcount{$k} ";
			}
			my $m=$n;
			$inf =~ s/\s+$//;
			$liginf{$m}=$inf;		
    }

    open(OUT, ">$saveBsites_dat");
    foreach my $k(@keys)
    {
			print OUT $pred{$k} . "\t" . $liginf{$k} . "\n";
    }
    close(OUT);   

    open(OUT, ">$saveBsites_clr");
    foreach my $k(@keys)
    {
			my @a=@{$clr{$k}};
			print OUT "cluster $k\n";
			foreach my $j(@a)
			{
	    	if($j =~ /(\S+)\t(.+)/)
	    	{
					my $site=$1;
					my $inf=$2;
					my $lig=$bsligname{$site};
					print OUT "$site\_$lig\t$inf\n";
	    	}
			}
    }
    close(OUT);

    my @prob=();
    open(IN, "$saveBsites_prob_dat");
    while(my $line=<IN>)
    {
			chomp($line);
			if($line =~ /\d+/)
			{
	    	push(@prob, $line);
			}
    }
    close(IN);

    open(OUT, ">$saveBsites_prob_dat");
    foreach my $k(@keys)
    {
			print OUT "$prob[$k]\n";
    }
    close(OUT); 
}

### Ligand Chemical similarity
sub tannimoto
{
    my($f1, $f2)=@_; 
    my @f1bin = split(//, $f1);
    my @f2bin = split(//, $f2);
    my $bins  = $#f1bin;
    $bins     = $#f2bin if($#f2bin < $bins);
    my $NA=0; 
    my $NB=0; 
    my $NAB=0;  
    my $Tan=0;
    for(my $i=0; $i<=$bins; $i++)
    {
			$NA++  if($f1bin[$i] ne '0');
			$NB++  if($f2bin[$i] ne '0');
			$NAB++ if(($f1bin[$i] eq $f2bin[$i]) && ($f1bin[$i] ne '0'))
    }
    my $d = ($NA+$NB-$NAB);
    $Tan = sprintf("%4.3f",($NAB/$d)) if($d > 0) ;
    
    return $Tan;
}

sub select_templates
{
    my ($filename)=@_;   

    open(FIN, $filename) or die "$filename is missed!\n";    
    while(my $line1 = <FIN>)
    {
			chomp($line1);
		
			if($line1 =~ /(\S+)/)
			{	   	    
		    my @inf=split(/:/, $line1, -1);
		    my $site=$inf[0];
				my $pdb=substr($site, 0, 5);
	    	my $bres=$bslib{$site};    	   
				my $raw     = $inf[1];
		    my $sia     = $inf[2];
		    my $cov     = $inf[3];
		    my $sib     = $inf[4];
		    my $lcov    = $inf[5];
		    my $predBSR = $inf[6];

		    if($predBSR !~ /\d+/) 
	  	  {
					#print "ignore $site $predBSR\n";
					next;
	    	}

	    	$predBSR =~ s/,$//;
	    	my @pred=split(/,/, $predBSR);
	    	my $cons=0;
	    	foreach my $p(@pred)
	    	{
					if(!exists $JSD{$p})
					{
		    		print "wrong res $p\n";
		    		next;
					}
					$cons += $JSD{$p};
	    	}
	    if(@pred>0)
	    {
				$cons=$cons/@pred;
	    }

	    my $q = 2.0/(1.0 + exp( -0.5*$raw*$cov - 0.5*$sib*$lcov - 0.2*$cons - 0.00*$sia) )-1.0;
	    $Qtemplate{$site}=$q;
      if($q>$q_cut)
	    {   
				my @pBSR        = split(/,/, $predBSR);
				$quality{$site} = $q;
				$pbsr{$site}    = $predBSR;
				$feature{$site} = sprintf "%.2f %.2f %.2f %.2f %.2f %.2f %.2f", $q, $raw, $sia, $cov, $sib, $lcov, $cons;
	    }
		}
	}
  close(FIN);
}

sub select_templates1
{
    my ($filename)=@_;   

    my @allt = sort {$Qtemplate{$b} <=> $Qtemplate{$a}} keys %Qtemplate;
    my %topT=();
    for(my $i=0; $i<$nTemplates1; $i++)
    {
			$topT{$allt[$i]}=1;
    }

    open(FIN, $filename) or die "$filename is missed!\n";    
    while(my $line1 = <FIN>)
    {
			chomp($line1);

			if($line1 =~ /(\S+)/)
			{	   	    
	    	my @inf=split(/:/, $line1, -1);

	    	my $site=$inf[0];
	    	
	    	if(!exists $topT{$site})
	    	{
					next;
	    	}
	   
		    my $pdb=substr($site, 0, 5);
		    my $bres=$bslib{$site};    	   
				my $raw     = $inf[1];
		    my $sia     = $inf[2];
		    my $cov     = $inf[3];
		    my $sib     = $inf[4];
		    my $lcov    = $inf[5];
		    my $predBSR = $inf[6];

		    if($predBSR !~ /\d+/) 
		    {
					#print "ignore $site $predBSR\n";
					next;
		    }

		    $predBSR =~ s/,$//;
		    my @pred=split(/,/, $predBSR);
		    my $cons=0;
		    foreach my $p(@pred)
		    {
					if(!exists $JSD{$p})
					{
		    		print "wrong res $p\n";
		    		next;
					}
					$cons += $JSD{$p};
	    	}
	    	if(@pred>0)
	    	{
					$cons=$cons/@pred;
	    	}

	    	my $q=$Qtemplate{$site};         
	    	my @pBSR        = split(/,/, $predBSR);
	    	$quality{$site} = $q;
	    	$pbsr{$site}    = $predBSR;
	    	$feature{$site} = sprintf "%.2f %.2f %.2f %.2f %.2f %.2f %.2f", $q, $raw, $sia, $cov, $sib, $lcov, $cons;	
			}
    }
    close(FIN);
}

sub frquency
{
    my ($ALN_ref, $Nseq_ref, $AA_ref, $SW_ref, $pseudocount)=@_;
    my %align   = %$ALN_ref; 
    my $Nseq    = $$Nseq_ref;
    my %weights = %$SW_ref;  
    my %AA2in   = %$AA_ref;
    
    my @Qres   = split(//, $align{1, 2});
    my $Ncol   = $#Qres;
    my %res_count=();

   
    my $Qresno=0;
    my %Qmapping=();
    for(my $j=0; $j<=$#Qres; $j++)
    {
        $res_count{$j}=0;
        if($Qres[$j] ne '-')
        {
            $Qresno++;
            $Qmapping{$Qresno}=$j;
        }
    }

    
    my @ARR=();
    for(my $i=1; $i<=$Nseq; $i++)
    {
        my @res=split(//, $align{$i, 2});
        for(my $j=0; $j<=$#res; $j++)
        {
            $ARR[$i][$j]=$res[$j];
        }
    }
    my $AAcount = keys %AA2;
    my %AA_freq=();
    my %sum_seq_weights=();
    my $k=0;

    for(my $j=0; $j<=$Ncol; $j++)
    {
				if($Qres[$j] eq '-')
				{
	    			next;
				}
				$k++;	
        foreach my $key (@AA)
        {
            $AA_freq{$k, $key}=0;
        }
				my $w=0;
        for(my $i=1; $i<=$Nseq; $i++)
        {
	    			if(!exists $AA2{$ARR[$i][$j]})
	    			{
							next;
	    			}
	    			$w += $weights{$i};	    
            $AA_freq{$k, $ARR[$i][$j]} += $weights{$i}; ##weighted frequency in clolumn $j
        }
        foreach my $key (@AA)
        {
            $AA_freq{$k, $key} = ($AA_freq{$k, $key}+ $pseudocount)/($w + $AAcount * $pseudocount); ##I change it here	    
        }
    }

    return %AA_freq;
}

# Get Seqeuence Weights based on Henikoff-Henikoff schema
sub load_sequence_weights
{
    my ($ALN_ref, $Nseq_ref, $AA_ref, $seq_weights_ref)=@_;

    my %align  = %$ALN_ref;
    my $Nseq   = $$Nseq_ref;
    my %AA2in  = %$AA_ref;
    
    my $weights = 0;

    my %NRC=();my %RC=();my %seen=();
    for(my $i=1; $i<=$Nseq; $i++)
    {
        my @res=split(//, $align{$i, 2});

        for(my $j=0; $j<=$#res; $j++)
        {
            my $AAN=$AA2in{$res[$j]};
            $RC{$j, $AAN}++;      #number of times $AAN appears at the jth column of the alignment
            if(exists $seen{$j, $AAN})
            {
                next;
            }
            else
            {
                $seen{$j, $AAN}=1;
                $NRC{$j}++;  #total number of AA type in the jth column
            }
        }
    }

    for(my $i=1; $i<=$Nseq; $i++)
    {
        my @res=split(//,$align{$i,2});

        my $ali_L=@res;
        for(my $j=0; $j<=$#res; $j++)
        {
            my $AAN=$AA2in{$res[$j]};
            $$seq_weights_ref{$i} += 1/($NRC{$j} * $RC{$j, $AAN});
        }
        $$seq_weights_ref{$i} = sprintf("%6.4f", $$seq_weights_ref{$i}/$ali_L); ##normalization
        $weights += $$seq_weights_ref{$i};
    }
    return $weights;
}


sub smooth_ss
{
    my ($ssref, $len)=@_;

    my %sec=%$ssref;
    
    #smooth single  --x-- => -----
    for(my $i=3; $i<=$len-2; $i++)
    {
				if($sec{$i}==2 || $sec{$i}==4)
				{
	    		my $j=$sec{$i};
	    		if($sec{$i-2} != $j)
	    		{
						if($sec{$i-1} != $j)
						{
		    			if($sec{$i+1} != $j)
		    			{
								if($sec{$i+2} != $j)
								{
			    					$sec{$i}=1;
								}
		    			}
						}
	    		}
				}
    }

    #   smooth double 
    #   --xx-- => ------
    for(my $i=1; $i<=$len-5; $i++)
    {
				#helix
				if($sec{$i} != 2)
				{
	    		if($sec{$i+1} != 2)
	    		{
						if($sec{$i+2} == 2)
						{
		    			if($sec{$i+3} == 2)
		    			{
								if($sec{$i+4} != 2)
								{
								    if($sec{$i+5} != 2)
								    {
												$sec{$i+2}=1;
												$sec{$i+3}=1;
								    }
								}
		    			}
						}
	    		}
				}
	
				#beta
				if($sec{$i} != 4)
				{
				    if($sec{$i+1} != 4)
				    {
								if($sec{$i+2} ==4)
								{
								    if($sec{$i+3} == 4)
								    {
											if($sec{$i+4} != 4)
											{
											    if($sec{$i+5} != 4)
											    {
															$sec{$i+2}=1;
															$sec{$i+3}=1;
											    }
											}
					    			}
									}
	    				}
					}
    }
    

    #smooth connect
    for(my $i=1; $i<=$len-2; $i++)
    {
			if($sec{$i} == 2)
			{
			    if($sec{$i+1} != 2)
			    {
						if($sec{$i+2} == 2)
						{
						    $sec{$i+1}=2;
						}
			    }
			}
			elsif($sec{$i} == 4)
			{
			    if($sec{$i+1} != 4)
			    {
						if($sec{$i+2} == 4)
						{
						    $sec{$i+1}=4;
						}
			    }
			}
    }  

    
    return %sec;
}

sub z_scores
{
     my($score_ref, $z_score_ref)=@_;
     my %r_score  =%$score_ref;

     my @keys = keys %r_score;
     my $Lch = @keys;
     my $z1_a  = 0; my $z1_a2=0; my %zscore=();
     foreach my $i(@keys)
     {
         $z1_a  += $r_score{$i};
         $z1_a2 += $r_score{$i}**2;
     }
     $z1_a/=$Lch;
     $z1_a2/=$Lch;
     my $z1_sd = sqrt($z1_a2 - $z1_a**2);

     foreach my $i(@keys)
     {
         $$z_score_ref{$i}=($r_score{$i}-$z1_a)/$z1_sd;
     }
}

sub vote
{
    my ($rst_ref, $saveBsites_dat, $saveBsites_prob_dat, $saveBsites_clr)=@_;

    my @rst=@$rst_ref;



    open(BOUT1, ">$saveBsites_dat");
    open(BOUT2, ">$saveBsites_prob_dat");
    open(BOUT3, ">$saveBsites_clr");
    foreach my $r (@rst)
    {
			chomp($r);
			if($r =~ /^cluster\s+(\d+):(.+)/)
			{
			    $n_cluster++;
			    my $c=$1;
			    my $line=$2;
			    $clusters{$n_cluster}=$line;
			    my @a=split(/\s+/, $line);
			    
			    #make consensus prediction in the cluster @a
			    for(my $k=1; $k<=$Lch; $k++)
			    {
						$prob{$k}=0.0;
			    }
			    
			    my $score=0;
			    my $q0=-1;
			    my $pcount=0;
			    print BOUT3 "cluster $n_cluster\n";
	    		foreach my $i(@a)
	    		{
							$pcount++;
							
							my $site = $keys[$i];
							my @pred = split(/,/, $pbsr{$site});
							print BOUT3 "$site\t$feature{$site} @pred\n";
							
							if($quality{$site}>$q0)
							{
							    $q0=$quality{$site};
							}
							$score+=$quality{$site};
		
							for(my $k=0; $k<@pred; $k++)
							{
							    $prob{$pred[$k]} += 1;
							}
	    		}
			    $score = $score/@a;
			    for(my $k=1; $k<=$Lch; $k++)
			    {
						$prob{$k} /= $pcount if $pcount>0;
						printf BOUT2 "%.2f ", $prob{$k};
			    }
	   			 print BOUT2 "\n";
	    
	    
			    my @pres=();
			    my %zscore=();
			    &z_scores(\%prob, \%zscore);
	    
	    		my @keyr = reverse sort{$zscore{$a} <=> $zscore{$b}} keys %zscore;
	    		my $count=0;
	    		my $jsd=0;
	    
	    		foreach my $k(@keyr)
	    		{   	
						if($zscore{$k}>$zcut)
						{		
		    			push(@pres, $k);
		    			$jsd += $JSD{$k};
		    			$count++;
		    			if($count>=$bsrncut ) #for hard target, make at most min(10, 0.03*$Lch) predictions
		    			{
								#last;
		    			}
						} 
	    		}
			    if($count==0)
			    {
						$zcut=0;
						foreach my $k(@keyr)
						{   	
						    if($zscore{$k}>$zcut)
						    {		
									push(@pres, $k);
									$jsd += $JSD{$k};
									$count++;
									if($count>=3) #for hard target, make at most min(10, 0.03*$Lch) predictions
									{
									    last;
									}	 
						    }    
						}
			    }
	    
	   			$jsd = $jsd/$count if $count>0;
	    
	    		my $size  = @a/@keys;
	    		my $size1 = log(1+@a ** (1.0/3));
	    		$cscore{$n_cluster} = 2/(1+exp( -$size*$q0 - 0.1*$size1 -0.2*$jsd)  ) -1;
	    
	    		
	    		@keyr= sort {$a <=> $b} @pres;
	    		my $ns=@a;
	    		printf BOUT1 "%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t", $cscore{$n_cluster}, $q0, $size, $size1, $jsd;
	    		print BOUT1 join(",", @keyr) . "\n";
			}
    }
    close(BOUT1);
    close(BOUT2);
    close(BOUT3);    
}
