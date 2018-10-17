# MIonSite
MIonSite: an accurate predictor for targeting the binding sites of Zn2+, Ca2+, Mg2+, Mn2+, Fe3+, Cu2+, Fe2+, Co2+, Na+, K+, Cd2+, and Ni2+ from protein sequence.

There are four folders, i.e., "standardalone_MIonSite", "Training_dataset", "Independent_validation_dataset", and "BTD", which are the standardalone program of MIonSite, the training protein dataset, the independent testing protein dataset, and the blind test dataset, respectively.

The steps of installing the MIonSite program are clearly described in the README.md file of the "standardalone_MIonSite" folder.


For each type of metal ions, the training dataset and independent validation dataset are organized as follows:
	
Taking Zn2+ as an example, its training dataset consists of one file:

	Training_dataset\ZN_seq.fasta 
	
and its corresponding independent validation dataset also consists of two files:

 	Indepedent_validation_dataset\ZN_seq.fasta 



For each protein in BTD, we only know its metal ion-binding site, but we do not know the type of metal ion.
