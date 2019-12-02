

# READ MAPPING for the LyzMcre transgene (cre recombinase):

## OVERVIEW
To validate that the Cpt2KO genotype is indeed set-up correctly we can look at read mapping and expression of the LyzM gene (synonomous to Lys2). This is straight forward to achieve from the standard "ExonCoverage" analysis and can also be assessed from the featureCounts output, as the Lys2 gene is present in the mouse reference assembly. However, the transgene 'cre recombinase' must be checked separately, as it is a foreign gene from enterobacteria phage 1, that is inserted after LyzM to enable conditional gene knockouts in myeloid tissue (the mouse lineage with this transgene is referred to as LyzMcre, and it can be crossed with floxed mice lineages to achieve tissue specific gene deletion). To check for expression of cre recombinase, we need to create a new dB (reference) file to map to, or append the cre gene to an existing dB file. Rather than mapping the entire RNAseq dataset to the new/altered dB, we can extract the portion of reads that failed to map to the mouse assembly and attempt to remap them with the new dB. This will drastically reduce processing time, as only a fraction of the reads did not map.

## ANTICIPATED WORKFLOW:

Extract unmapped reads from .bam and/or .sam files (samtools view) >>
Convert .bam/.sam to fasta or fastq, or whatever file type that BLASTn accepts >>
Download and install BLAST+ and magicBLAST software >>
Download cre recombinase fasta file >>
Create BLAST dB from cre fasta file with NCBI tool makeblastDb >>
map reads to cre gene per sample with magicBLAST tool (which accepts fastq or fasta)

## RESOURCES:
* samtools

		See these pages for help configuring the command:
		http://www.htslib.org/doc/samtools.html # commands and parameters
		http://broadinstitute.github.io/picard/explain-flags.html # use FLAGs to subset mapped and unmapped data
		

* BLAST software - https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download

* Cre recombinase sequence - https://www.ncbi.nlm.nih.gov/nuccore/NC_005856.1?report=fasta&from=436&to=1467



### DETAILS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

```bashscript

# creating cre recombinase dB on local environment (since magicBLAST software is not on server) 
 
 cd set/working/directory/creRecombinase
 
 makeblastdb -in cre.fa -input_type fasta -parse_seqids -title creRecombinase -logfile cre_makeblsdB_log -dbtype nucl


# extracting unmapped RNAseq reads (on server)
  
find . -type f -name *.sorted.bam | xargs -I % bsub -q research-hpc -M 61000000 -R 'rusage[mem=63000]' -n 8 -a 'docker(EDIT)' -oo %.unmapped.bsub.log 'samtools view -b -f 4 -o %.unmapped.bam %'

bjobs | cut -b -7 | xargs echo -n # print current jobs


# converting .bam back to fastq format
  
find . -type f -name '*.unmapped.bam' | xargs -I % bsub -q research-hpc -M 61000000 -R 'rusage[mem=63000]' -n 8 -a 'docker(EDIT)' -o %.fastq 'samtools fastq %'

# NOTE: for some reason bsub does not handle this command like others. For example -oo option for a bsub.log output does not produce the stder file. If you make the option -o then it will produce the actual .fastq.  The command 'samtools fastq' itself has no standard output option such as -o, instead they say to us the '>' to write to a file, and there are lots of options for splitting up or merging .bam files. So the above command works fine, but no log files are produced. LATER learned, the stder information is included at the beginning, end, and a little bit within the body of the fastq file. This information needs to be removed for downstream processing of the fastq file.


# sending just the .fastq files to local environment
  
rsync -zt --progress remote/server/location/*.fastq' /local/location/fastq/
  

#  testing magicBLAST with a fastq
  
magicblast -query reads.fastq -db my_reference -infmt fastq -outfmt tabular -no_unaligned -out 
```

 ### NOTE: The resulting fastq is almost ready for magicblast but needs to be cleaned in three ways
* remove all header lines associated with the bsub run
* somewhere near the end of the file, in a different place each time, an insertion of [M::bam2fq_mainloop] processed XXX reads] occurs. This insertion breaks up the fastq formating and causes a parsing error. Remove just this insertion and the read will run perfectly fine
		
#### EXAMPLE OF INSERTION:
		
		DDDCDIIIIHHHIHIIIIHIDHGFHIIIIIIIIIIIHHIIIIIIHIHIIH
		@SN1063:704:H5FJHBCX2:1:2216:9907:22278
		TCCTAGGTTTTTTGTT[M::bam2fq_mainloop] processed 549952 reads
		ATTCCAGATGAATTTGCAAATTGCTCCTTCTAAT
		+
		
		EXAMPLE OF APPROPRIATE CORRECTION:
		
		DDDCDIIIIHHHIHIIIIHIDHGFHIIIIIIIIIIIHHIIIIIIHIHIIH
		@SN1063:704:H5FJHBCX2:1:2216:9907:22278
		TCCTAGGTTTTTTGTTATTCCAGATGAATTTGCAAATTGCTCCTTCTAAT
		+
		
* remove the last lines of the file pertaining to bsub run reporting
