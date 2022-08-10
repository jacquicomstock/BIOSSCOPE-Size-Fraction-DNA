#create shell script to run phyloassigner
vi phylo_plastid.sh

#paste the following into the created shell script
python pythonassigner_v0.9.py \
    --out_dir plastid_output/ \
    --ref_align databases/Cyanos_and_Plastid_refDB/ref.aln \
    --ref_tree databases/Cyanos_and_Plastid_refDB/ref_tree.txt \
    --query_seqs plastid.fasta \
    --mapping databases/Cyanos_and_Plastid_refDB/edge.mapping \
    --threads 32 \
    --placer pplacer

#move appropriate files into the folder
mv databases/ BIOSSCOPE_Frac/
mv pythonassigner_v0.9.py BIOSSCOPE_Frac/

#for plastid data there are too many ASVs (1500) to run upfront on the cluster, need to submit job to slurm
#this is the code to submit a batch job
sbatch \
	--job-name=plastid_phyloassigner \
	--nodes=1 \
	--tasks-per-node=32 \
	--cpus-per-task=1 \
	--mem=60G \
	--time=00:30:00 \
	--output=plastid_out \
	--error=plastid_err \
	--wrap="bash BIOSFrac_plastid.sh"
  
#can do the above steps for all taxa types I want further classification with through Phyloassigner
  
#move files back
mv databases/ 
mv pythonassigner_v0.9.py 

#move phyloassigner export files to local community
scp -r carlsonlab@pod.cnsi.ucsb.edu:/home/carlsonlab/PhyloAssigner_python_UCSB-main/plastid_output/ /home/Mobaxterm/Desktop/
