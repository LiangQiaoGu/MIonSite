# MIonSite
MIonSite: an accurate metal ion binding site predictor

For each type of metal ions, the training dataset and independent validation dataset are organized as follows:
	
Taking Zn2+ as an example, its training dataset consists of two files:

        Training\ZN_seq.fasta 
	Training\ZN_label.fasta
and its corresponding independent validation dataset also consists of two files:

 	Validation\ZN_seq.fasta 
	Validation\ZN_label.fasta	

In all the label files, '0' represents a non-binding site, while '1' represents a metal ion binding site.
