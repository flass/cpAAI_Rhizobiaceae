# cpAAI_Rhizobiaceae
A pipeline and reference protein sequence data for generating core-proteome alignment and computing core-proteome amino-acid identity (cpAAI) values from bacterial genomes of the *Rhizobiaceae* family, to assess their membership to an existing or a new genus.

## Citation
These software/data are in support of this paper:

Nemanja Kuzmanović, Camilla Fagorzi, Alessio Mengoni, Florent Lassalle, George C diCenzo. 2021. **Taxonomy of Rhizobiaceae revisited: proposal of a new framework for genus delimitation.** *BioRxiv* 2021.08.02.454807.
[https://doi.org/10.1101/2021.08.02.454807]

## Dependencies
- Python 3 (recommended version >= 3.8.5), with package:
	- Bio (BioPython) (recommended version >= 1.78)
- NCBI Blast+ (recommended version >= 2.6.0+)
- Clustal Omega (recommended version >= 1.2.3)  
  or
- MAFFT (recommended version >= v7.487) (preferred)
- R (recommended version >= 4.0.2), with package:
	- phangorn

## Description
The `data` folder contains reference marker protein sequence sets and their corresponding alignments for 170 non-recombining core protein markers conserved among a set of 97 *Rhizobiaceae* genomes, including all available genomes of type strains of species described in the *Rhizobiaceae* family (as of April 2021).  
The script `genome2cpAAI.py` in the `pipeline` folder can be used to obtain homolous sequence to the reference marker set described above, or to any analogous dataset, for instance a core-proteome dataset defined for another prokaryotic taxon.

## Pipeline usage
The script `genome2cpAAI.py` takes as input a list of complete genomic nucleotide sequence (contig) files from the query organisms (argument `-q`), a list of files containing the reference marker protein sequence sets (argument `-p`) and a path for the folder where the output will be written (argument `-o`). Optionally, a list of files containing the alignments of the reference marker protein sequence sets can be provided (argument `-A`). All sequence and alignment files should be in Fasta format.  
This script will locate in the input genome the loci coding for the protein homologs of the reference marker proteins using `tblastn`. It will then extract these sequences, translate them, and align them with the reference sequence set using `mafft` (the aligner `clustalo` can be used alternatively if specified through the argument `--aligner`); if the reference alignments were provided, these are used as a alignmentguide using `mafft -add` algorithm, allowing faster processing.  
Finally, the resulting core protein alignments are concatenated into a single file.

Here is an example of how to use the script `genome2cpAAI.py`:
```sh
# generate list of reference marker protein sequence and alignment files
# lists of files for the 170 Rhizobiaceae core protein markers are already available in data/
ls data/protein_sequences/*.faa > protein_sequences_list
ls data/protein_alignments/*.aln > protein_alignments_list
mkdir run_genome2cpAAI/ tmp/
# run the pipeline
genome2cpAAI.py -q testgenomelist -p protein_sequences_list \
 -A data/protein_alignments_list -o run_genome2cpAAI \
 --threads 8 --tmp tmp
```

Then, the output concatenated alignment `concatenated_marker_proteins.aln` (located in the specified result folder, here `run_genome2cpAAI/`) can be used to compute a core-proteome tree using the phylogenetic program of your choice (task not included in this package) and to compute the **cpAAI values** between query and reference strains.  
This can be done with the following `R` script:
```R
library(phangorn)
nfconcatprotaln = "run_genome2cpAAI/concatenated_marker_proteins.aln"
concatprotaln = as.matrix(read.FASTA(nfconcatprotaln, type='AA'))
cpaai = 100 * (1 - as.matrix(dist.aa(concatprotaln, scaled=TRUE)))
write.table(cpaai, "cpAAI_matrix.txt", sep='\t', col.names=NA, row.names=T)
# the tree obtained based on the concatenate can be used to order the matrix
nfrootedtree = "yourNewicktreefile"
rootedconcattree = read.tree(nfrootedtree)
rcpaai = cpaai[rootedconcattree$tip.label,rootedconcattree$tip.label]
write.table(rcpaai, "orderred_cpAAI_matrix.txt", sep='\t', col.names=NA, row.names=T)
```
## Visualisation
...