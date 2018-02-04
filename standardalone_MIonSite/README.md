MIonSite: Identifying Metal Ion-Binding Sites via Enhanced AdaBoost Algorithm with Protein Sequence Information

Note: The MIonSite can be only used to predict the Zn, Ca, Mg, Mn, and Fe binding sites from protein sequence.

Installation:
	
	1. Download the MIonSite package and decompression 

				$ tar xvzf MIonSite.tar.gz
			
	   Let $HOME_FOLDER to the absolute path of the MIonSite folder 
	
	2. Install the SANN software, which is available on $HOME_FOLDER or http://lee.kias.re.kr/~protein/wiki/doku.php?id=sann:download:sann
			
	3. Configure the $HOME_FOLDER/java/Config.properties file, as follows:
	
				BLAST_BIN_DIR=$HOME_FOLDER/blast-2.2.26
				BLASTPGP_EXE_PATH=$HOME_FOLDER/blast-2.2.26/blastpgp
				BLASTPGP_DB_PATH=$HOME_FOLDER/nr/nr
				PSIPRED321_FOLDER_DIR=$HOME_FOLDER/psipred321
				BLASTPGP_OUTPUT_PARSER_DIR=$HOME_FOLDER/blastpgpOutputPARSER
				SANN_RUNNER_PATH=$HOME_FOLDER/SANN/sann/bin/sann.sh
				MODEL_FOLDER=$HOME_FOLDER/models


Example:
	
   			$ cd $HOME_FOLDER/java 
   			$ java MIonSite $HOME_FOLDER/example/example.fasta $HOME_FOLDER/example/output ZN

References
[1] xxx
[2] xxx
[4] Keehyoung Joo, Sung Jong Lee and Jooyoung Lee, Sann: Solvent accessibility prediction of proteins by nearest neighbor method, PROTEINS-STRUCTURE FUNCTION AND BIOINFORMATICS, 80, 1791-1797 (2012)
