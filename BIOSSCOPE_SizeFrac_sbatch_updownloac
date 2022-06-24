# These commands were typed directly into the terminal command line

sbatch \
	--job-name=BIOS_dada2 \
	--nodes=1 \
	--tasks-per-node=32 \
	--cpus-per-task=1 \
	--mem=120G \
	--time=24:00:00 \
	--output=dada2_out \
	--error=dada2_err \
	--wrap="Rscript BIOS_dada2.R"

scp -r fastq-20220601T002923Z-001.zip jcomstock@pod.cnsi.ucsb.edu:/home/jcomstock/BIOS_Frac/
scp -r fastq-20220601T002923Z-002.zip jcomstock@pod.cnsi.ucsb.edu:/home/jcomstock/BIOS_Frac/
scp -r fastq-20220601T003324Z-001.zip jcomstock@pod.cnsi.ucsb.edu:/home/jcomstock/BIOS_Frac/
scp -r fastq-20220601T003435Z-001.zip jcomstock@pod.cnsi.ucsb.edu:/home/jcomstock/BIOS_Frac/
scp -r silva_nr99_v138.1_wSpecies_train_set.fa.gz jcomstock@pod.cnsi.ucsb.edu:/home/jcomstock/BIOS_Frac/

#to unzip all .zip files
unzip '*.zip'


#to import to carlsonlab cluster
scp -r fastq-20220601T002923Z-001.zip carlsonlab@pod.cnsi.ucsb.edu:/home/carlsonlab/BIOS_Frac/
scp -r fastq-20220601T002923Z-002.zip carlsonlab@pod.cnsi.ucsb.edu:/home/carlsonlab/BIOS_Frac/
scp -r fastq-20220601T003324Z-001.zip carlsonlab@pod.cnsi.ucsb.edu:/home/carlsonlab/BIOS_Frac/
scp -r fastq-20220601T003435Z-001.zip carlsonlab@pod.cnsi.ucsb.edu:/home/carlsonlab/BIOS_Frac/
scp -r silva_nr99_v138.1_wSpecies_train_set.fa.gz carlsonlab@pod.cnsi.ucsb.edu:/home/carlsonlab/BIOS_Frac/

# download files to the local computer
scp -r jcomstock@pod.cnsi.ucsb.edu:/home/jcomstock/BIOS_Frac/dada2_output/ /home/mobaxterm/Desktop/


#import specialized fasta files
scp -r BIOSFrac_SAR11.fasta jcomstock@pod.cnsi.ucsb.edu:/home/jcomstock/Phyloassigner/PhyloAssigner_python_UCSB-main/
scp -r BIOSFrac_SAR202.fasta jcomstock@pod.cnsi.ucsb.edu:/home/jcomstock/Phyloassigner/PhyloAssigner_python_UCSB-main/
scp -r BIOSFrac_Cyano.fasta jcomstock@pod.cnsi.ucsb.edu:/home/jcomstock/Phyloassigner/PhyloAssigner_python_UCSB-main/

# download phyloassigner otput files to the local computer
scp -r jcomstock@pod.cnsi.ucsb.edu:/home/jcomstock/Phyloassigner/PhyloAssigner_python_UCSB-main/BIOSSCOPE_Frac/ /home/mobaxterm/Desktop/dada2_output
