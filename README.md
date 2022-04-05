# HUPAN: HUman Pan-genome ANalysis

---

 **1. Introduction**
 
The human reference genome is still incomplete, especially for those population-specific or individual-specific regions, which may have important functions. It encourages us to build the pan-genome of human population. Previously, our team developed a "map-to-pan" strategy, [EUPAN][1], specific for eukaryotic pan-genome analysis. However, due to the large genome size of individual human genome, [EUPAN][2] is not suit for pan-genome analysis involving in hundreds of individual genomes. Here, we present an improved tool, HUPAN (HUman Pan-genome ANalysis), for human pan-genome analysis.

The HUPAN homepage is http://cgm.sjtu.edu.cn/hupan/

The HUPAN paper is available at https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1751-y

**2. Installation**

**Requirements** 

 - R 3.1 or later (https://www.r-project.org/)
    
    R is utilized for visualization and statistical tests in HUPAN
    toolbox. Please install R first and make sure R and Rscript are
    under your PATH.

 - R packages 

    Several R packages are needed including ggplot2, reshape2
    and ape packages. Follow the Installation step, or you can install the packages by yourself.

**Installation procedures** 

 - Download the HUPAN toolbox from [github][3]:

    `git clone git@github.com:SJTU-CGM/HUPAN.git`
 
 - Alternatively, you also could obtain the toolbox in the [HUPAN][4]
   website and uncompress the HUPAN toolbox package:
 
    `tar zxvf HUPAN-v**.tar.gz`

 - Install necessary R packages:

    `cd HUPAN & Rscript installRPac`
 
 - Compile necessary tools:
 
    `make`

    You will find executable files: *ccov*, *bam2cov* and *hupan* et al. in bin/ directory.

 - Add bin/ to PATH and add lib/ to LD_LIBRARY_PATH. To do this, add the
   following text to ~/.bash_profile:
  
       export PATH=$PATH:/path/to/HUPAN/bin:
       export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/HUPAN/lib/:
       export PERL5LIB=$PERL5LIB:/path/to/HUPAN/lib/:

 - and run: 
   
    `source ~/.bash_profile`

 - Test if HUPAN toolbox is installed successfully: `hupan`. If you see
   the following content, congratulations! HUPAN toolbox is successfully
   installed. If not, see if all the requirements are satisfied or
   contact the authors for help.

   
     Usage: hupan <command> ...
        
        Available commands (in order of use):
         qualSta      	View the overall sequencing quality of a large number of fastq files
         trim         	Trim or filter low-quality fastq reads parallelly
         assemble     	Assemble reads parallelly
         alignContig   Align assembled reads to the reference genome using nucmer
         extractSeq   	Extract potentially unaligned contigs parallelly
         assemSta     	Checks statistics of potentially unaligned contigs
         getUnalnCtg  	Extract the unaligned contigs from nucmer alignment (processed by quast)
         mergeUnalnCtg Merges unaligned contigs from multiple samples
         rmRedundant  	Remove redundant contigs within a fasta file using CD-HIT
         blastAlign   	Align sequences to nucleotide database using BLAST
         getTaxClass  	Obtain the taxonomic classification of sequences
         rmCtm        	Detect and discard the potential contamination
         splitSeq     	Split sequence file into multiple smaller sized files
         genePre      	Ab initio gene prediction of the novel sequences using RNA and protein evidence
         mergeNovGene 	Merge MAKER results from multiple result files
         filterNovGene	Filter the novel precited genes
         pTpG         	Get the longest transcripts to represent genes
         alignRead    	Map reads to a reference parallelly
         alignContig  	Map assembled genomes to the referenece using BWA or Bowtie2
         sam2bam      	Convert alignments (.sam) to sorted and indexed .bam files
         geneCov      	Calculate gene body coverage and CDS coverage
         geneExist    	Determine gene presence-absence based on gene body coverage and CDS coverage
         sim          	Simulate and plot the pan-genome and the core genome
         
        Other available commands not utilised in the main pipeline (but that can still be used for extra functionality):
         fastaSta     	Simple script to calculate number of sequences and bases within a single fasta file
         simSeq       	Simulate and plot the total size of novel sequences
         bamSta       	Calculates the coverage of the genomes using .bam files and Qualimap
         bam2bed      	Calculate genome region presence-absence from .bam
         subSample    	Select subset of samples from gene PAV profile
         gFamExist    	Determine gene family presence-absence based on gene presence-absence
        	

**3.	Main analysis procedures**

Below, please find the procedure and exact HUPAN commands used for my (Jess's) research. The original instructions from the HUPAN authors can be found on the HUPAN GitHub page. 
For my research, we used a SLURM high-performance computing facility, which is one of the systems HUPAN is designed for.

In this research, multiple different ancestral populations, consisting of 25 - 35 samples each, were run separately through the pipeline to obtain population-specific non-reference sequences. These sequences were then merged in an additional redundancy-removing step to obtain the non-redundant set of non-reference sequences from all the populations present.

Each ancestral population was contained in its own directory and each command below was run in each directory respectively, until the merging step.

Notes:
* The trimming and quality control checking of the fastq files for each sample was performed separately to the pipeline.
* For most steps, the number of threads used can be specified with the `-t` flag. I usually utilised 8 or 16, depending on the state of the cluster.
* Before starting the pipeline, I created a `HUPANdatabases` directory that contains all the necessary databases and reference sequences and annotations needed for the pipeline. In this directory there is:
  * The human reference genome fasta file saved in `ref/ref.fa` (downloaded from https://gatk.broadinstitute.org/hc/en-us/articles/360035890951-Human-genome-reference-builds-GRCh38-or-hg38-b37-hg19)
  * The human reference genome transcript file saved in `ref/ref.transcripts.fa` (downloaded from https://www.gencodegenes.org/human/releases.html, release 35)
  * The human reference genome annotation GTF file saved in `ref/ref.genes.gtf` (downloaded from https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.chr_patch_hapl_scaff.annotation.gtf.gz, release 38)
  * NCBI's BLAST non-redundant nucleotide database saved as `/blastindex/data/blastdb.nal`, downloaded and set up as follows:
    ```
    mkdir ncbidownloads & cd ncbidownloads
    wget https://ftp.ncbi.nih.gov/blast/db/FASTA/nt.gz | gunzip & cd ..
    hupanSLURM blastAlign mkblastdb ncbidownloads blastindex /cbio/bin  # /cbio/bin is where the BLAST executable is
    ```
  * NCBI's taxonomy database saved in `taxonomyinfo/nucl_gb.accession2taxid`, downloaded and set up as follows:
    ```
    cd ncbidownloads
    wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid
    wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz & tar -zvxf new_taxdump.tar.gz
    mkdir taxonomyinfo && mv ncbidownloads/nucl_gb.accession2taxid taxonomyinfo && mv ncbidownloads/new_taxdump/rankedlineage.dmp taxonomyinfo
    ```
    
**(1) *De novo* assembly of individual genomes**

This step assembles the fastq files of each genome into a whole genome using SGA. The assembled sequences are saved in a directory called `01_assembled`.

This step has been changed from the original HUPAN process; a Nextflow pipeline is now used, developed by Gerrit Botha and myself. The Nextflow pipeline repository can be found at: https://github.com/h3abionet/dec-2020-hackathon-stream1/tree/main/assembly

**(2) Aligning samples to the reference genome**

* This step aligns the assembled genome to the reference genome using the nucmer tool in MUMmer. 
* The aligned sequences are saved in a folder called `02_aligned`.
* `/cbio/bin` is the path to the MUMmer, nucmer and show-cords executable files, which all run through a MUMmer Singularity container.
```
hupanSLURM alignContig -t 16 01_assembled 02_aligned /cbio/bin /cbio/projects/015/HUPANdatabases/ref/ref.fa
```
* This will result in `.delta` and `.coords` files for each sample.
* hupanSLURM alignContig code is found in `/cbio/projects/015/HUPANdatabases/HUPAN/lib/HUPANrmHighSLURM.pm`

**(3) Extracting contigs with a low similarity to the reference genome**

* This step identifies and discards all those contigs with > 95% identity and contig length alignment to the reference genome.
* The sequences with low similarity (and that are potentially non-reference sequences) are saved in a folder called `03_candidate`.
* The final argument (`02_aligned`) is the path to the MUMmer/nucmer results, created in the previous step.
* This step is run interactively and produces output in real time, so first run `screen` and then use an interactive node: `srun --nodes=1 --ntasks=1 --mem=50g -t 24:00:00 --pty bash`
```
hupanSLURM extractSeq 01_assembled 03_candidate 02_aligned
```
* This will result in a `.candidate.unaligned.contig` file for each sample.
* hupanSLURM extractSeq code is found in `/cbio/projects/015/HUPANdatabases/HUPAN/lib/HUPANrmHighSLURM.pm`

**(4) Assessing the candidate sequences using QUAST**

* This steps maps the assembled contigs in the candidate folder against the reference and checks the statistics of these contigs using QUAST.
* The QUAST results of the candidate contigs are saved in a folder called `04_quastresult`.
* `/cbio/projects/015/HUPANdatabases/images/` is where the QUAST exectuable is located, which runs a QUAST Singularity container.
```
hupanSLURM assemSta -t 16 03_candidate/data 04_quastresult /cbio/projects/015/HUPANdatabases/images/ /cbio/projects/015/HUPANdatabases/ref/ref.fa
```
hupanSLURM assemSta code is found in `/cbio/projects/015/HUPANdatabases/HUPAN/lib/HUPANassemStaSLURM.pm`

**(5) Collecting unaligned (both fully and partially) contigs**

* This steps collects all unaligned (both fully and partially) contigs from the 03_candidate folder using information from the `04_quastresult` folder.
* The final unaligned sequences will be saved in folder called `05_unaligned`.
* The HUPAN notes on GitHub say to use `-p` to specify .contig suffix, but the actual .pm file says to use `-s`.
```
hupanSLURM getUnalnCtg -s .contig 03_candidate/data 04_quastresult/data 05_unaligned
```
This will result in `.fully.contig`, `partially.contig` and `partially.coords` files for each sample.
hupanSLURM getUnalnCtg code is found in `/cbio/projects/015/HUPANdatabases/HUPAN/lib/HUPANunalnCtgsSLURM.pm`

**(6) Merging all the unaligned contigs**

* This step merges all the unaligned contigs in `05_unaligned`.
* The resulting files will be saved in the folder `06_mergeunaligned`.
* Partially and fully unaligned contigs are still kept separate, and will each result in one file each.
```
hupanSLURM mergeUnalnCtg 05_unaligned/data 06_mergeunaligned
```
This will result in total.fully.fa, total.partilly.coords and total.partilly.fa files with all the merged data.
hupanSLURM mergeUnalnCtg code is found in `/cbio/projects/015/HUPANdatabases/HUPAN/lib/HUPANunalnCtgsSLURM.pm`

**(7) Removing redundancy**

* This step removes redundancy (i.e. similar sequences) using CD-HIT. 
* It does this by forming 'clusters' of similar contigs at a defined identity threshold and choosing the longest contig in that cluster as the representative contig. 
* The command must be performed twice; once for fully unaligned sequences, and once for partially unaligned sequences. 
* They must be run in a separate directory, called `07_rmredundant`. The restults will then be stored in directories within this directory.
* The default % identity for clustering using CD-HIT in HUPAN is 90%, but you can change it by using the `-c` flag, e.g. `-c 0.85`
* `/cbio/bin` is the location of the cd-hit-est exectuable, which runs a Singularity container.
```
mkdir 07_rmredundant && cd 07_rmredundant
hupanSLURM rmRedundant cdhitCluster -t 8 ../06_mergeunaligned/total.fully.fa rmredundant.fully /cbio/bin
hupanSLURM rmRedundant cdhitCluster -t 8 ../06_mergeunaligned/total.partilly.fa rmredundant.partially /cbio/bin/
```
* This will result in `cluster_info.txt` and `non-redundant.fa` files each for both fully and partially unaligned sequences.
* hupanSLURM rmRedundant code is found in `/cbio/projects/015/HUPANdatabases/HUPAN/lib/HUPANrmRdtSLURM.pm`

**(8) Identifying contamination in the unaligned sequences**

* This step identifies contamination by comparing all unaligned sequences to the BLAST nucleotide database. 
* The database has been previously downloaded and set up according to HUPAN GitHub instructions (found above).
* The aligned results will be saved in a directory called `08_blastresult`. Within this directory, there will be a `rmredundant.fully` and a `rmredundant.partially` directory.
* `/cbio/bin` is where the BLAST executable file is found.
```
hupanSLURM blastAlign blast -t 16 07_rmredundant 08_blastresult /cbio/projects/015/HUPANdatabases/blastindex/data/blastdb /cbio/bin
```
* This will result in `non-redundan.blast` files.
* hupanSLURM blastAlign code is found in `/cbio/projects/015/HUPANdatabases/HUPAN/lib/HUPANblastAlignSLURM.pm`

**(9) Obtaining the taxonomic classifications of the sequences**

* This steps looks at all the BLAST hits and classifies them into their correct taxonomies (if there are any).
* It uses two databases downloaded from NCBI: a taxonomy database, and a database that converts accession numbers to taxonomic IDs. Both have already been downloaded and set up (instructuions found above).
* The command has to be run twice, once each for partially and fully unaligned sequences.
* Each command will create its own directory; `09.1_taxclassfully` for fully unaligned reads and `09.2_taxclasspartially` for partially unaligned reads.
```
hupanSLURM getTaxClass 08_blastresult/data/rmredundant.fully/non-redundan.blast /cbio/projects/015/HUPANdatabases/taxonomyinfo 09.1_taxclassfully
hupanSLURM getTaxClass 08_blastresult/data/rmredundant.partially/non-redundan.blast /cbio/projects/015/HUPANdatabases/taxonomyinfo 09.2_taxclasspartially
```
* hupanSLURM getTaxClass code is found in `/cbio/projects/015/HUPANdatabases/HUPAN/lib/HUPANgetTaxClassSLURM.pm`

**(10) Removing contaminating sequences**

* This step removes the contaminating sequences (microbiological and non-primate eukaryotes) from the non-reference sequences.
* It produces fasta files of the sequences removed, so we can get information on what contaminants were there.
* The command has to be run twice, once eact for partially and fully unaligned sequences.
* The `-i 60` flag is the percentage identity limit for the contamination. Default is 90% for HUPAN, so we must specify.
* The first file is the non-redundant sequence file, second is the BLAST alignment file, third is the accession file for the taxonomy, and the fourth is the output directory.
```
hupanSLURM rmCtm -i 60 07_rmredundant/rmredundant.fully/non-redundant.fa 08_blastresult/data/rmredundant.fully/non-redundan.blast 09.1_taxclassfully/accession.name 10.1_rmctmfully
hupanSLURM rmCtm -i 60 07_rmredundant/rmredundant.partially/non-redundant.fa 08_blastresult/data/rmredundant.partially/non-redundan.blast 09.2_taxclasspartially/accession.name 10.2_rmctmpartially
```
* hupanSLURM rmCtm code is found in `/cbio/projects/015/HUPANdatabases/HUPAN/lib/HUPANrmContaminateSLURM.pm`

**(11) Combining fully and partially unaligned reads**

* This steps removes redundancy again by merging the final fully and partially unaligned sequences, again using CD-HIT.
* A directory must first be made and then the fully and partially unaligned sequences combined and moved into this directory. This is done in the first command below.
* This combined file is then run through CD-HIT to remove redundancy (at 90% identity) and the final non-redundant non-reference contigs are saved in `12_finalpangenome`.
* `/cbio/bin` is the location of the CD-HIT executable.
```
mkdir 11_nonreference && cat 10.1_rmctmfully/novel_sequence.fa 10.2_rmctmpartially/novel_sequence.fa > 11_nonreference/nonrefernce.before.fa
hupanSLURM rmRedundant cdhitCluster -t 8 11_nonreference/nonrefernce.before.fa 12_finalpangenome /cbio/bin/
```
* This will result in a `tmp` directory and `cluster_info.txt` and `non-redundant.fa` files.
* hupanSLURM rmRedundant code is found in `/cbio/projects/015/HUPANdatabases/HUPAN/lib/HUPANrmRdtSLURM.pm`

**(12) Merging the different populations**

**At this point, each population has its own non-redundant non-reference sequences. Now we can combine them and remove redundancy.**

**By doing this, we are ending up with the same sequences we would have had if we had run all the populations together, but now we have the population-specific information separated and can easily perform downstream or upstream analysis on population differences.**
* To do this, use `cat` to combine all the sequences into one file, and then use CD-HIT to remove redundancy at 90%.
* The population directories should all be in the same location and can therefore be accessed using *.
* The final non-reference sequences are saved in `allpopulations/`
```
cat */12_finalpangenome/non-redundant.fa > ../allpopulations/11_nonreference/all.nonreference.before.fa
cd ../allpopulations
hupanSLURM rmRedundant cdhitCluster -t 8 11_nonreference/all.nonreference.before.fa 12_finalpangenome /cbio/bin/
```
* This will result in the final non-redundant non-reference sequence file in `allpopulations/12_finalpangenome/non-redundant.fa` for all populations combined.
* All steps from this point are performed in the `allpopulations/` directory.

**(13) Splitting the sequences into smaller files**

* This step splits the non-redundant non-reference sequences into smaller files so that MAKER can handle them easier and quicker.
* It takes the fasta file and saves it as smaller files in the second directory given.
* It splits the fasta file into files that contain roughly 2 Mbp, set using the `-m 2000000` flag
* The output directory is `13_genepredinput`.
```
hupanSLURM splitSeq -m 2000000 12_finalpangenome/non-redundant.fa 13_genepredinput
```
* This will result in a number of `part*.fa` files.
* The number of files will depend on how many sequences were in the `non-redundant.fa` file.
* hupanSLURM splitSeq code is found in `/cbio/projects/015/HUPANdatabases/HUPAN/lib/HUPANsplitSeqSLURM.pm`

**(14) Predicting genes from the non-reference sequences**

* This step takes the split up non-reference sequence files and inputs them into MAKER.
* The config is the MAKER configuration file where we have specified all of the variables and datasets to use. A copy of this config file can be found in this repository.
* `/cbio/bin` is the location of the MAKER executable, which runs a MAKER Singularity container.
```
hupanSLURM genePre -t 16 13_genepredinput 14_genepred /cbio/projects/008/jess/HUPANrun/AfricanPanGenome/config /cbio/bin
```
* This will result in various MAKER predictions for each `part*.fa` file.
* hupanSLURM genePre code is found in `/cbio/projects/015/HUPANdatabases/HUPAN/lib/HUPANgenePreSLURM.pm`

**(15) Merging the novel predictions**

* This step merges the novel predictions from the different parts.
* It takes in the directory where all the different parts' results are (`14_genepred`) and places them into single files in `15_genepredmerge`.
* `/cbio/bin` is where the MAKER executable is.
* This step is run interactively and produces output in real time (takes about 10 minutes), so first run `screen` and then use a minimal interactive node: `srun --pty bash`
```
hupanSLURM mergeNovGene 14_genepred/result/ 15_genepredmerge /cbio/bin
```
* This will result in `all.maker.map`, `combine.all.maker.gff`, `combine.all.maker.proteins.fasta` and `combine.all.maker.transcripts.fasta` files.
* hupanSLURM mergeNovGene code is found in `/cbio/projects/015/HUPANdatabases/HUPAN/lib/HUPANmergeNovGene.pm`

**(16) Filtering the gene predicitons**

* This step filters and refines the novel predicted genes based on HUPAN specifications. These specifications can be found in the HUPAN paper supplementary methods.
* The first command copies the `non-redundant.fa` sequence file from `12_finalpangenome` to the `15_genepredmerge` directory, because the script needs it to work. Call it `novel.fa`.
* The `ref/` directory contains exactly what the HUPAN script looks for with the exact names, i.e. `ref.fa` and `ref.transcripts.fa` (explained above).
* The last three inputs are (in order) the locations of the BLAST, CD-HIT and RepeatMasker executables respectively.
```
cp 12_finalpangenome/non-redundant.fa 15_genepredmerge/novel.fa
hupanSLURM filterNovGene -t 16 15_genepredmerge 16_genepredfilter /cbio/projects/015/HUPANdatabases/ref /cbio/bin/ /cbio/bin/ /cbio/bin/
```
* hupanSLURM filterNovGene code is found in `/cbio/projects/015/HUPANdatabases/HUPAN/lib/HUPANfilterNovGeneSLURM.pm`

**(17) Assembling the pan-genome**

* This step has does two things. First:
  * It takes in the human reference annotation file in `ref/ref.genes.fa` and retains only the longest transcript of each annotation while discarding others, like mRNA, lncRNA, etc.
  * The `-f` flag just tells the script that the file is in gtf format, not gff.
  * It saves the annotations in a file called `ref.genes-ptpg.gtf`, found in the `/cbio/projects/015/HUPANdatabases/ref` directory. 
  * The annotation file was downloaded from https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.chr_patch_hapl_scaff.annotation.gtf.gz (release 38).
    * The annotation file downloaded MUST be later than release 21, otherwise the HUPAN script won't work as it will be in a different format (format changed after release 21).
    * However, you also have to remove any annotations from patch/alt/decoy/unplaced/mitochondrial DNA contigs, as the HUPAN pipeline seems to be unable to recognise them. So remove from the annotation file all lines after "chrY" annotations. I have not shown this step below, as the exact line numbers may be different.
  * In summary: run the `pTpG` script on the annotation file, and then remove the extra annotations and save the final working file as `ref.genes-ptpg-primaryseqs.gtf`.
* Then secondly, we combine this annotation file and the novel sequence annotation (from step 16) into one annotation file
  * You must first convert the novel annotation file into GTF format too, or otherwise start with GFF format for the reference annotation.
* Both steps done in real time (interactively) so use `screen` and then start an interactive node with extra memory: `srun --mem=15g --pty bash`
```
hupanSLURM pTpG -f /cbio/projects/015/HUPANdatabases/ref/ref.genes.gtf /cbio/projects/015/HUPANdatabases/ref/ref.genes-ptpg.gtf
# insert step to remove any additional annotations and save the primary sequence annotations as ref.genes-ptpg-primaryseqs.gtf
# insert step to convert the Final.gff file to GTF format.
mkdir 17_pan && cat /cbio/projects/015/HUPANdatabases/ref/ref.genes-ptpg-primaryseqs.gtf 16_genepredfilter/data/Final/Final.gtf > 17_pan/pan.gtf
```
* hupanSLURM pTpG code is found in `/cbio/projects/015/HUPANdatabases/HUPAN/lib/HUPANpTpGSLURM.pm`

**(18) Aligning the samples to the pan-genome**

* This step does three things:
  * First, it combines the novel sequences with the human reference sequences into one file called `pan.fa`. This is the final pan-genome file.
  * It then uses the `bowtie2-build` function to build a Bowtie index from the pan-genome sequences. 
    * This is done interactively and requires a lot of memory, so a `screen` session and an interactive node are required for the second command: `srun --mem=50g --time=24:00:00 --pty bash`
* Finally, it then aligns all the fastq files (of each genome) to the pan-genome file using Bowtie 2. The fastq files of all samples are found in `00_fastq`.
* The `-f` flag specifies Bowtie2 (not BWA), and the `-s` and `-k` flags specify the naming convention of the fastq files given (e.g. `LZP0B_R1.fastq.gz`)
* `18_map2pan` is the output directory.
* `/cbio/projects/015/HUPANdatabases/images/` is where I have stored the Bowtie 2 Singularity container.
```
cat /cbio/projects/015/HUPANdatabases/ref/ref.fa 12_finalpangenome/non-redundant.fa > 17_pan/pan.fa
cd 17_pan && /cbio/projects/015/HUPANdatabases/images/bowtie2-build pan.fa pan && cd ..
hupanSLURM alignRead -f bowtie2 -t 8 -s .fastq.gz -k _R 00_fastq 18_map2pan /cbio/projects/015/HUPANdatabases/images/ 17_pan/pan
```
* The final command will result in a `.sam file` for each fastq file (genome/individual) given.
* hupanSLURM alignRead code is found in `/cbio/projects/015/HUPANdatabases/HUPAN/lib/HUPANmapSLURM.pm`

**(19) Converting to SAM format**

* This step uses Samtools to convert the `.sam` files to `.bam` format and then sort and index the `.bam` files.
* It takes in the `.sam` files in the `18_map2pan` directory and outputs the sorted and indexed `.bam` files to the `19_panBam` directory.
* `/cbio/bin` is where the Samtools executable is found.
```
hupanSLURM sam2bam -t 8 18_map2pan/data 19_panBam /cbio/bin
```
* This will result in a `.bam` and `.bai` file for each fastq file (genome/individual) given.
* hupanSLURM sam2bam code is found in `/cbio/projects/015/HUPANdatabases/HUPAN/lib/HUPANsamToBamSLURM.pm`

**(20) Calculating gene coverage**

* This step calculates the coverage of each gene or CDS in all of the genomes used.
* It does this by looking at the alignment information in the indexed and sorted `.bam` files from the previous step and identifying whether these alignments fall within the annotations in the `pan.gtf` file.
* The results will be output to the `20_geneCov` directory.
```
hupanSLURM geneCov -t 8 19_panBam/data 20_geneCov/ 17_pan/pan.gtf
```
* This will result in a `.sta` file for each genome.
* hupanSLURM geneCov code is found in `/cbio/projects/015/HUPANdatabases/HUPAN/lib/HUPANgeneCovSLURM.pm`

**(21) Determing gene presence-absence variance profile**

* This steps checks whether a gene is present or absent (based on a defined coverage threshold) for each sample.
* The gene coverage information from the `20_geneCov` directory is first combined into two files, one recording gene coverage and one recording CDS coverage.
  * The first input is the prefix you want to use for the merged file, and the second input is that data directory where the `.sta` files are.
* The final HUPAN command takes in the two files outputted from the previous command and creates the gene presence-absence variance profile.
* It does this by checking the gene or CDS coverage at a certain threshold. HUPAN sets theirs at 95% gene coverage.
* For this command, the two integer values are (in order) the minimum gene coverage and the minimum CDS coverage.
```
cd 20_geneCov && hupanSLURM mergeGeneCov summary_ data && cd ..
mkdir 21_geneExist && hupanSLURM geneExist 20_geneCov/summary_gene.cov 20_geneCov/summary_cds.cov 0 0.95 > 21_geneExist/cds95.gene.exist
```
* The first command will result in `geneCovmergedgene.cov` and `geneCovmergedcds.cov` files.
* hupanSLURM mergeGeneCov code is found in `/cbio/projects/015/HUPANdatabases/HUPAN/lib/HUPANgeneCovSLURM.pm`
* The second command will result in the `gene.exist file`.
* hupanSLURM geneExist code is found in `/cbio/projects/015/HUPANdatabases/HUPAN/lib/HUPANgeneExistSLURM.pm`

And now you're finally done!

**Performing simulations using the HUPAN suite**

* The HUPAN pipelines offers two simulation outputs on your data:
  * The `sim` function simulates the size of the core and pan-genomes from the gene presence-absence data generated in step 21.
  * The `simSeq` function simulates and plots the total amount of novel sequence as more individuals are added.

* For the `sim` function:
  * A `screen` session and a minimum interactive node are required: `srun --pty bash`
  * You must choose whichever `.gene.exist` file is most apprpropriate for analysis of the pan-genome. Here, gene coverage of 95% is chosen.
```
mkdir 22_simulations
hupanSLURM sim 21_geneExist/gene95.gene.exist 22_simulations/PAV
```
* This will result in a `sim_out.txt` file and a `sim_out.txt_PAV_plot.pdf` file with the simulated data and plots in them.

* For the `simSeq` function:
  * The first command performs the simulation on the sample-specific unaligned data in `05_unaligned`. 
  * `/cbio/bin` is the location of the CD-HIT executable which is used for the simulation.
  * This simulation will take a few days to run depending on the number of samples.
  * The second command just plots the result and is interactive, so a `screen` session and a minimum interactive node are required: `srun --pty bash`
```
hupanSLURM simSeq simNovelSeq -t 8 05_unaligned/data/ 22_simulations/novelSeqs /cbio/bin
hupanSLURM simSeq plotNovelSeq 22_simulations/novelSeqs/ outputPlot
```
* This will result in a `simNovelSeq.pdf` file with the simulated plot.
