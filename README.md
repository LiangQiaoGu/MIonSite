# MIonSite
MIonSite: an accurate metal ion binding site predictor

For each type of metal ions, the training dataset and independent validation dataset are organized as follows:
	
Taking Zn2+ as an example, its training dataset consists of two files:

	Training_dataset\ZN_seq.fasta 
	Training_dataset\ZN_label.fasta
and its corresponding independent validation dataset also consists of two files:

 	Indepedent_validation_dataset\ZN_seq.fasta 
	Indepedent_validation_dataset\ZN_label.fasta	

In all the label files, '0' represents a non-binding site, while '1' represents a metal ion binding site.
