#!/usr/bin/perl -w
use strict;
use Getopt::Long;

# input
my $system    = $^O;
my $protfasta_path = "";
my $sequence = "";
my $blast_out_path = "";
my $psipredSS_path = "";
my $parse_blast_pl_path = "";
# output
my $blast_fasta = "";   
my $saveQueyPSFM = "";

##############################################################
### input parameters
unless (lc($system) eq "linux")
{
    printf("Your operating system $system is unsupported at this time\n");
    printf("Currently only Linux is supported\n");
    exit();
}
GetOptions('protfasta_path:s' => \$protfasta_path, 'sequence:s' => \$sequence, 'blast_out_path:s' => \$blast_out_path, 'psipredSS_path:s' => \$psipredSS_path, 'parse_blast_pl_path:s' => \$parse_blast_pl_path, 'blast_fasta:s' => \$blast_fasta, 'saveQueyPSFM:s' => \$saveQueyPSFM);

if(!$protfasta_path || !$sequence || !$blast_out_path || !$psipredSS_path || !$parse_blast_pl_path || !$blast_fasta || !$saveQueyPSFM)
{
    print "\nprepareQueryPSFM usage:\nprepareQueryPSFM.pl arguments\n";
    print "====================\n";
    print "Mandatory arguments:\n";
    print "====================\n";
    print "-protfasta_path protfasta_path\n";
    print "-sequence sequence\n";
    print "-blast_out_path blast_out_path (blastpgp> output)\n";
    print "-psipredSS_path psipredSS_path (psipred output file)\n";
    print "-parse_blast_pl_path parse_blast_pl_path (where the parse_blast.pl file is located)\n";
		print "-blast_fasta blast_fasta (where you want to save the blast.fasta file)\n";
		print "-saveQueyPSFM saveQueyPSFM (where you want to save the saveQueyPSFM file)\n";

    exit();
}

###############################################################################
## Global parameters
my @AA=qw(A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V );
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

my $Lch=length $sequence;

###############################################################################
## parse the psipred predicted secondary structure information

open(SS, "$psipredSS_path");
my $j=0;
my %sec=();
#<SS>; 
<SS>;
while(my $line=<SS>)
{
   if($line=~/\s*(\d+)\s+\S+\s+(\S+)/)
   {
			my $j=$1;
			my $ps=$2;
			$sec{$j}=$SS2num{$ps};        
   }
}
close(SS);
%sec = &smooth_ss(\%sec, $Lch);


###############################################################################
## compute query protein frequency (PSFM)
printf "compute frequency .....\n";
my %close_freq = &wfreq("$blast_out_path");
open(PRO, ">$saveQueyPSFM");
for(my $i=1; $i<=$Lch; $i++)
{
    my $res=substr($sequence, $i-1, 1);
    printf PRO "$res $sec{$i} ";

    foreach my $A(@AA)
    {
			printf PRO "%.4f ", $close_freq{$i, $A};
    }

    printf PRO "\n";
}
close(PRO);


###############################################################################
## 


####################################################################################
### sub function
sub wfreq
{	
    my ($blast)=@_;

    `$parse_blast_pl_path -fas $blast_out_path -q $protfasta_path $blast_fasta`;
    my %ALN=();
    my $Pcount=0;
    open(ALN,"$blast_fasta") || die "Cant open $blast_fasta";
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

    my %freq=();
    if($Pcount > 1)
    {
        $Pcount=1000 if ($Pcount > 1000);
        my %seq_weights = ();
				&load_sequence_weights(\%ALN, \$Pcount, \%AA2index, \%seq_weights);
    

				##compute weighted frequency now
				my $PSEUDOCOUNT= 0.0000001;
				%freq = &frquency(\%ALN, \$Pcount, \%AA2index, \%seq_weights, $PSEUDOCOUNT);
    }
    else
    {
	 		my @Qres   = split(//, $ALN{1, 2});
	 		for(my $j=0; $j<@Qres; $j++)
	 		{
	    	foreach my $key (@AA)
	     	{
		 			$freq{$j+1, $key}=0;
	     	}
	 		}
    }

    return %freq;
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
