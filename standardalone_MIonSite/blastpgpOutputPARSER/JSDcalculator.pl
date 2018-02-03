#!/usr/bin/perl -w
use strict;
use Getopt::Long;

# input
my $system    = $^O;
my $query_sequence			= "";			# the query protein sequence
my $blast_fasta_path 		= ""; 		# which is the output file of parse_blast.pl 
my $save_jsd_score_path = "";			# where do you want to save the JSD_score

##############################################################
### input parameters
unless (lc($system) eq "linux")
{
    printf("Your operating system $system is unsupported at this time\n");
    printf("Currently only Linux is supported\n");
    exit();
}
GetOptions('query_sequence:s' => \$query_sequence, 'blast_fasta_path:s' => \$blast_fasta_path, 'save_jsd_score_path:s' => \$save_jsd_score_path);

if(!$query_sequence || !$blast_fasta_path || !$save_jsd_score_path)
{
    print "\nJSDcaculator usage:\nJSDcaculator.pl arguments\n";
    print "====================\n";
    print "Mandatory arguments:\n";
    print "====================\n";
    print "-query_sequence query_sequence\n";
    print "-blast_fasta_path blast_fasta_path (where you want to save the blast.fasta file)\n";
    print "-save_jsd_score_path saveQsave_jsd_score_pathueyPSFM (where you want to save the save_jsd_score_path file)\n";

    exit();
}

#####################################################
# Initial SETTINGS (Do not change below this line)
my %seqQ=();
my $Lch=length $query_sequence;  
for(my $i=1;$i<=$Lch;$i++)
{
    my $seq1=substr($query_sequence, $i-1, 1);
    $seqQ{$i}=$seq1;
}


my @AA=
    (
     "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET","PHE",
     "PRO", "SER", "THR", "TRP", "TYR", "VAL", "ASX", "GLX", "UNK");

my %AA2index=
    (
     'A'=>'1', 'R'=>'2', 'N'=>'3', 'D'=>'4', 'C'=>'5', 'Q'=>'6', 'E'=>'7', 'G'=>'8', 'H'=>'9', 'I'=>'10', 'L'=>'11', 'K'=>'12', 
     'M'=>'13', 'F'=>'14', 'P'=>'15', 'S'=>'16', 'T'=>'17', 'W'=>'18', 'Y'=>'19', 'V'=>'20', 'B'=>'21', 'Z'=>'22', '-'=>'23');	
my %index2AA=
    (
     '1'=>'A','2'=>'R','3'=>'N','4'=>'D','5'=>'C','6'=>'Q','7'=>'E','8'=>'G','9'=>'H','10'=>'I','11'=>'L','12'=>'K', 
     '13'=>'M','14'=>'F','15'=>'P','16' =>'S','17'=>'T','18'=>'W','19'=>'Y','20'=>'V','21'=>'B','22'=>'Z','23'=>'-');	
	
my %ts=
    (
     'ALA'=>'A', 'ARG'=>'R', 'ASN'=>'N', 'ASP'=>'D', 'CYS'=>'C', 'GLN'=>'Q', 'GLU'=>'E', 'GLY'=>'G',
     'HIS'=>'H', 'ILE'=>'I', 'LEU'=>'L', 'LYS'=>'K', 'MET'=>'M', 'PHE'=>'F', 'PRO'=>'P', 'SER'=>'S', 
     'THR'=>'T', 'TRP'=>'W', 'TYR'=>'Y', 'VAL'=>'V', 'ASX'=>'N', 'GLX'=>'Q', 'UNK'=>'G',
     'A'=>'ALA', 'R'=>'ARG', 'N'=>'ASN', 'D'=>'ASP', 'C'=>'CYS', 'Q'=>'GLN', 'E'=>'GLU', 'G'=>'GLY',
     'H'=>'HIS', 'I'=>'ILE', 'L'=>'LEU', 'K'=>'LYS', 'M'=>'MET', 'F'=>'PHE', 'P'=>'PRO', 'S'=>'SER', 
     'T'=>'THR', 'W'=>'TRP', 'Y'=>'TYR', 'V'=>'VAL', 'B'=>'ASX', 'Z'=>'GLX', 'X'=>'UNK');
	
     # A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X
my @BLOSUM62 = 
    (
     [ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,-2,-1, 0],
     [-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3, 0, 1,-1],
     [-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3, 6, 0,-1],
     [-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3, 1, 0,-1],
     [ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2],
     [-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2, 0, 5,-1],
     [-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2, 0, 2,-1],
     [ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3, 0,-2,-1],
     [-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3, 1, 0,-1],
     [-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-3,-3,-1],
     [-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-3,-2,-1],
     [-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2, 0, 1,-1],
     [-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-2, 0,-1],
     [-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-3,-3,-1],
     [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-2,-1,-2],
     [ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2, 1, 0, 0],
     [ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0, 0,-1, 0],
     [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,-4,-2,-2],
     [-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-2,-1,-1],
     [ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-3,-2,-1],
     [-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3, 6, 0,-1],
     [-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2, 0, 5,-1],
     [ 0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-1,-1,-1]);

# This is the BLOSUM62 background distribution.
my %blsm_bck=
    (
     'A'=>0.078, 'R'=>0.051, 'N'=>0.041, 'D'=>0.052, 'C'=>0.024, 'Q'=>0.034, 'E'=>0.059, 'G'=>0.083, 
     'H'=>0.025, 'I'=>0.062, 'L'=>0.092, 'K'=>0.056, 'M'=>0.024, 'F'=>0.044, 'P'=>0.043, 'S'=>0.059, 
     'T'=>0.055, 'W'=>0.014, 'Y'=>0.034, 'V'=>0.072, 'B'=>0.041, 'Z'=>0.034, 'X'=>0.083,);
	 
my @BLOSM_FREQ = 
    ( [0.0215,0.0023,0.0019,0.0022,0.0016,0.0019,0.0030,0.0058,0.0011,0.0032,0.0044,0.0033,0.0013,0.0016,0.0022,0.0063,0.0037,0.0004,0.0013,0.0051,0.0019,0.0019,0.0058],
      [0.0023,0.0178,0.0020,0.0016,0.0004,0.0025,0.0027,0.0017,0.0012,0.0012,0.0024,0.0062,0.0008,0.0009,0.0010,0.0023,0.0018,0.0003,0.0009,0.0016,0.0020,0.0025,0.0017],
      [0.0019,0.0020,0.0141,0.0037,0.0004,0.0015,0.0022,0.0029,0.0014,0.0010,0.0014,0.0024,0.0005,0.0008,0.0009,0.0031,0.0022,0.0002,0.0007,0.0012,0.0141,0.0015,0.0029],
      [0.0022,0.0016,0.0037,0.0213,0.0004,0.0016,0.0049,0.0025,0.0010,0.0012,0.0015,0.0024,0.0005,0.0008,0.0012,0.0028,0.0019,0.0002,0.0006,0.0013,0.0037,0.0016,0.0025],
      [0.0016,0.0004,0.0004,0.0004,0.0119,0.0003,0.0004,0.0008,0.0002,0.0011,0.0016,0.0005,0.0004,0.0005,0.0004,0.0010,0.0009,0.0001,0.0003,0.0014,0.0004,0.0003,0.0008],
      [0.0019,0.0025,0.0015,0.0016,0.0003,0.0073,0.0035,0.0014,0.0010,0.0009,0.0016,0.0031,0.0007,0.0005,0.0008,0.0019,0.0014,0.0002,0.0007,0.0012,0.0015,0.0073,0.0014],
      [0.0030,0.0027,0.0022,0.0049,0.0004,0.0035,0.0161,0.0019,0.0014,0.0012,0.0020,0.0041,0.0007,0.0009,0.0014,0.0030,0.0020,0.0003,0.0009,0.0017,0.0022,0.0035,0.0019],
      [0.0058,0.0017,0.0029,0.0025,0.0008,0.0014,0.0019,0.0378,0.0010,0.0014,0.0021,0.0025,0.0007,0.0012,0.0014,0.0038,0.0022,0.0004,0.0008,0.0018,0.0029,0.0014,0.0378],
      [0.0011,0.0012,0.0014,0.0010,0.0002,0.0010,0.0014,0.0010,0.0093,0.0006,0.0010,0.0012,0.0004,0.0008,0.0005,0.0011,0.0007,0.0002,0.0015,0.0006,0.0014,0.0010,0.0010],
      [0.0032,0.0012,0.0010,0.0012,0.0011,0.0009,0.0012,0.0014,0.0006,0.0184,0.0114,0.0016,0.0025,0.0030,0.0010,0.0017,0.0027,0.0004,0.0014,0.0120,0.0010,0.0009,0.0014],
      [0.0044,0.0024,0.0014,0.0015,0.0016,0.0016,0.0020,0.0021,0.0010,0.0114,0.0371,0.0025,0.0049,0.0054,0.0014,0.0024,0.0033,0.0007,0.0022,0.0095,0.0014,0.0016,0.0021],
      [0.0033,0.0062,0.0024,0.0024,0.0005,0.0031,0.0041,0.0025,0.0012,0.0016,0.0025,0.0161,0.0009,0.0009,0.0016,0.0031,0.0023,0.0003,0.0010,0.0019,0.0024,0.0031,0.0025],
      [0.0013,0.0008,0.0005,0.0005,0.0004,0.0007,0.0007,0.0007,0.0004,0.0025,0.0049,0.0009,0.0040,0.0012,0.0004,0.0009,0.0010,0.0002,0.0006,0.0023,0.0005,0.0007,0.0007],
      [0.0016,0.0009,0.0008,0.0008,0.0005,0.0005,0.0009,0.0012,0.0008,0.0030,0.0054,0.0009,0.0012,0.0183,0.0005,0.0012,0.0012,0.0008,0.0042,0.0026,0.0008,0.0005,0.0012],
      [0.0022,0.0010,0.0009,0.0012,0.0004,0.0008,0.0014,0.0014,0.0005,0.0010,0.0014,0.0016,0.0004,0.0005,0.0191,0.0017,0.0014,0.0001,0.0005,0.0012,0.0009,0.0008,0.0014],
      [0.0063,0.0023,0.0031,0.0028,0.0010,0.0019,0.0030,0.0038,0.0011,0.0017,0.0024,0.0031,0.0009,0.0012,0.0017,0.0126,0.0047,0.0003,0.0010,0.0024,0.0031,0.0019,0.0038],
      [0.0037,0.0018,0.0022,0.0019,0.0009,0.0014,0.0020,0.0022,0.0007,0.0027,0.0033,0.0023,0.0010,0.0012,0.0014,0.0047,0.0125,0.0003,0.0009,0.0036,0.0022,0.0014,0.0022],
      [0.0004,0.0003,0.0002,0.0002,0.0001,0.0002,0.0003,0.0004,0.0002,0.0004,0.0007,0.0003,0.0002,0.0008,0.0001,0.0003,0.0003,0.0065,0.0009,0.0004,0.0002,0.0002,0.0004],
      [0.0013,0.0009,0.0007,0.0006,0.0003,0.0007,0.0009,0.0008,0.0015,0.0014,0.0022,0.0010,0.0006,0.0042,0.0005,0.0010,0.0009,0.0009,0.0102,0.0015,0.0007,0.0007,0.0008],
      [0.0051,0.0016,0.0012,0.0013,0.0014,0.0012,0.0017,0.0018,0.0006,0.0120,0.0095,0.0019,0.0023,0.0026,0.0012,0.0024,0.0036,0.0004,0.0015,0.0196,0.0012,0.0012,0.0018],
      [0.0019,0.0020,0.0141,0.0037,0.0004,0.0015,0.0022,0.0029,0.0014,0.0010,0.0014,0.0024,0.0005,0.0008,0.0009,0.0031,0.0022,0.0002,0.0007,0.0012,0.0141,0.0015,0.0029],
      [0.0019,0.0025,0.0015,0.0016,0.0003,0.0073,0.0035,0.0014,0.0010,0.0009,0.0016,0.0031,0.0007,0.0005,0.0008,0.0019,0.0014,0.0002,0.0007,0.0012,0.0015,0.0073,0.0014],
      [0.0058,0.0017,0.0029,0.0025,0.0008,0.0014,0.0019,0.0378,0.0010,0.0014,0.0021,0.0025,0.0007,0.0012,0.0014,0.0038,0.0022,0.0004,0.0008,0.0018,0.0029,0.0014,0.0378]);


my %AA2sim=();
for(my $i=0;$i<23;$i++)
{
    for(my $j=0;$j<23;$j++)
    {
			my $x= $i+1;my $y= $j+1;
			$AA2sim{$index2AA{$x},$index2AA{$y}}= $BLOSUM62[$i][$j];
    }
}
my $PSEUDOCOUNT= 0.0000001;
my $lambda     = 0.5;
my $window     = 3;

#############################################################################
## Load the blast_out alignment profile to ALN
my %ALN=();
my $Pcount=0;
open(ALN,"<$blast_fasta_path") || die "Cant open $blast_fasta_path";
while(my $line=<ALN>)
{
    chomp($line);
    if($line =~ /^>(\S+)/)
    {
			my $Pname=$1;
			$Pcount++;
			my $Evalue= $1 if($line =~ /E=(\S+)/);	
			$ALN{$Pcount, 0}=$Pname;
			$ALN{$Pcount, 1}=$Evalue;
    }
    else
    {
			$line =~ s/X/-/g;  ###replace X by gap
			$line=~s/J/-/g;
			$line=~s/O/-/g;
			$line=~s/U/-/g;

			$ALN{$Pcount, 2}=$line;	
    }    
}
close(ALN);

#########################################################################
## Compute JSD score, then saving $save_jsd_score_path
printf "computing JSD socre...\n";
open(OUT1,">$save_jsd_score_path");
if($Pcount > 1)
{
    $Pcount=1000 if ($Pcount > 1000);
    my %seq_weights = ();
    my $tot_weights =   &load_sequence_weights(\%ALN, \$Pcount,\%AA2index,\%seq_weights); 
    
    my %jsd_scores = ();my %zscore1=(); my %wscore1=();      
    my %res_count = &js_divergence_score(\%ALN, \$Pcount, \%blsm_bck, \%AA2index, \%seq_weights, $PSEUDOCOUNT, $lambda, \%jsd_scores, $tot_weights);
   
    my @Qres   =split(//, $ALN{1, 2});
    my $Qresno=0;
    for(my $j=0; $j<=$#Qres; $j++)
    {
			if($Qres[$j] ne '-')
			{
	    	$Qresno++;
	    	printf OUT1 "%-4d\t%s\t%.6f\t%d/%d\n",$Qresno-1, $seqQ{$Qresno}, $jsd_scores{$Qresno}, $res_count{$j}, $Pcount;
			}
    }

    if($Qresno!=$Lch)
    {
			print "warning: length is not equal ($Qresno!=$Lch)";
    }
}else{
   	for(my $j=1; $j<=$Lch; $j++)
   	{
    	printf OUT1 "%-4d\t%s\t%.6f\t%4d/%4d\n",$j, $seqQ{$j}, 1.0, 1, 1;
   	}
}
close(OUT1);


################################################################################################
## sub function procedure

sub js_divergence_score
{
    my ($ALN_ref, $Nseq_ref, $dist_ref, $AA_ref, $SW_ref, $pseudocount, $lambda, $JSD_ref, $tot_weights)=@_;
    my %align   = %$ALN_ref;
    my $Nseq    = $$Nseq_ref;
    my %distr   = %$dist_ref;
    my %weights = %$SW_ref;  
    my %AA2in   = %$AA_ref;    

    my @Qres   =split(//,$align{1, 2});
    my $Ncol   =$#Qres;
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
	    	if($res[$j] eq $Qres[$j])
	    	{
					$res_count{$j}++;
	    	}
			}
    }
    my $AAcount = keys %AA2in;
    my %AA_freq=();my %sum_seq_weights=();my %R=();my %JSD_scores=();
    for(my $j=0; $j<=$Ncol; $j++)
    {
			foreach my $key (sort {$AA2in{$a} <=> $AA2in{$b}} keys %AA2in)
			{
	    	$AA_freq{$j, $key}=0;
	    	$R{$j, $key}      =0;
			}
			for(my $i=1; $i<=$Nseq; $i++)
			{
	    	$AA_freq{$j, $ARR[$i][$j]} += $weights{$i}; ##weighted frequency in clolumn $j
			}
			foreach my $key (sort {$AA2in{$a} <=> $AA2in{$b}} keys %AA2in)
			{	   
	    	$AA_freq{$j, $key} = ($AA_freq{$j, $key}+$pseudocount)/($tot_weights + $AAcount * $pseudocount); ##I change it here
	    	if($key ne '-')
	    	{
					$R{$j, $key}=($lambda*$AA_freq{$j,$key}) + ((1-$lambda)*$distr{$key});	   
	    	}
			}
    }
    for(my $j=0;$j<=$Ncol;$j++)
    {
			my $JSD_score=0;my $weighted_gap_penalty=0;my $gap_sum=0;
			foreach my $key (sort {$AA2in{$a} <=> $AA2in{$b}} keys %AA2in)
			{
	    	if($key eq '-'){next;} ##I add $lambda below
	    	if($R{$j, $key} != 0.0)
	    	{
					if($AA_freq{$j, $key} == 0.0)
					{
		    		print "warning AA frequency = 0!\n";
		    		$JSD_score += (1-$lambda)*($distr{$key} * log2($distr{$key}/$R{$j,$key}));
					}
					elsif($distr{$key} == 0.0)
					{
		    		print "warning backround frequency = 0!\n";
		    		$JSD_score += $lambda * ($AA_freq{$key} * log2($AA_freq{$key}/$R{$j,$key}));
					}
					else
					{
		    		$JSD_score += $lambda * ($AA_freq{$j,$key} * log2($AA_freq{$j,$key}/$R{$j,$key})) + (1-$lambda) * ($distr{$key} * log2($distr{$key}/$R{$j,$key}));
					}
	    	}
			}	
			##	$JSD_score = $JSD_score/2;
			for(my $i=1; $i<=$Nseq; $i++)
			{
	    	if($ARR[$i][$j] eq '-')
	    	{
					$gap_sum +=   $weights{$i}; 
	    	}
			}
			$weighted_gap_penalty = (1-  ($gap_sum/$tot_weights));
			$JSD_scores{$j}=($JSD_score*$weighted_gap_penalty);	 ##a simple process of gap penalty
    }
    for(my $j=1;$j<=$Qresno;$j++)
    {
			$$JSD_ref{$j}= sprintf("%5.4f",$JSD_scores{$Qmapping{$j}});
    }

    return %res_count;
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
	    	$RC{$j, $AAN}++;	  #number of times $AAN appears at the jth column of the alignment
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


sub log2 
{
    my $n = shift;
    return log($n)/log(2);
}

