# MIonSite
MIonSite: an accurate metal ion binding site predictor

There are three folders, i.e., "standardalone_MIonSite", "Training_dataset", and "Independent_validation_dataset", which is the standardalone program of MIonSite, the training protein dataset, and the independent testing protein dataset, respectively.

The steps of installing the MIonSite program are clearly described in the README.md file of the "standardalone_MIonSite" folder.


For each type of metal ions, the training dataset and independent validation dataset are organized as follows:
	
Taking Zn2+ as an example, its training dataset consists of two files:

	Training_dataset\ZN_seq.fasta 
	Training_dataset\ZN_label.fasta
and its corresponding independent validation dataset also consists of two files:

 	Indepedent_validation_dataset\ZN_seq.fasta 
	Indepedent_validation_dataset\ZN_label.fasta	

In all the label files, '0' represents a non-binding site, while '1' represents a metal ion binding site.
