# ShortBRED: Short, Better Representative Extract Dataset

ShortBRED is a pipeline to take a set of protein sequences, reduce them to a set of unique identifying strings ("markers"), and then search for these markers in metagenomic data and determine the presence and abundance of the protein families of interest.

For more information on the technical aspects to this program, or to cite ShortBRED, please use the following citation:

James Kaminski, Molly K. Gibson, Eric A. Franzosa, Nicola Segata, Gautam Dantas, and Curtis Huttenhower. High-Specificity Targeted Functional Profiling in Microbial Communities with ShortBRED." PLoS Computational Biology 11.12 (2015).

http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004557

## Installing ShortBRED

To install ShortBRED, simply download the file below, and extract its contents. 

[https://github.com/biobakery/shortbred/archive/0.9.4.zip](https://github.com/biobakery/shortbred/archive/0.9.4.zip)

You may also install it with Git by using the "Clone" command in the menu to your left.


To run ShortBRED, you will need the following dependencies:

* Python 2.7.9 
* Biopython v1.65
* ncbi-blast-2.2.28+
* usearch v6.0.307 (Please make sure this is up to date. Earlier versions of usearch use a different command for making database than what is expected by ShortBRED.)
* MUSCLE v3.8.31
* CD-HIT version 4.6


## Program Structure

ShortBRED consists of two main scripts:.

**ShortBRED-Identify** - This takes a FASTA file of amino acid sequences, searches for overlap among itself and against a separate reference file of amino acid sequences, and then produces a FASTA file of markers.

**ShortBRED-Quantify** - This takes the FASTA file of markers and quantifies their relative abundance in a FASTA file of nucleotide metagenomic reads.


## Input Files

**Proteins of Interest** - These are the sequences you are trying to find in the metagenomic data, saved as a FASTA file of amino acid sequences. ShortBRED will cluster them into families, and then create short markers for the families.

Note: Using a database of more than 2,500 AA sequences may require a large amount of time to run. Though the ARDB has ~8,000 sequences, it clusters down to ~900 sequences at 95% id and ShortBRED  typically needs to run overnight to search these against a large reference database (~500,000 sequences) using 10 threads on a typical cluster.

**Reference Set of Proteins** - ShortBRED compares your proteins of interest to this set of proteins, and eliminates short regions of overlap. ShortBRED removes these sections so that what remains from 

**Short Nucleotide Reads** - This is the metagenomic dataset you wish to analyze, typically a large FASTA file of containing millions of short nucleotide reads. 

---
## Running ShortBRED Identify to Create New Markers
## Method 1: Create new markers from a set of proteins of interest and a reference set of proteins.
 
To create markers for the sample data included with ShortBRED, set your current working directory to where you unpacked ShortBRED and type: 

```
$ ./shortbred_identify.py --goi  example/input_prots.faa --ref example/ref_prots.faa --markers mytestmarkers.faa --tmp example_identify
```

The sample data included with ShortBRED is quite small, so this command should run in less than a minute on a typical machine.  It creates a set of markers ("mytestmarkers.faa") that represents typical ShortBRED-Identify ouptut, but on a smaller scale.

## Specifying paths to dependencies

ShortBRED-Identify requires several dependencies to run. Since many users may have these programs saved under different names and in different paths, we include command line arguments to supply the path to these programs when they differ from the default.  

```
Programs:
  --usearch STRUSEARCH  Provide the path to usearch. Default call will be "usearch".
  --muscle STRMUSCLE    Provide the path to muscle. Default call will be "muscle".
  --cdhit STRCDHIT      Provide the path to usearch. Default call will be "cd-hit".
  --blastp STRBLASTP    Provide the path to blastp. Default call will be "blastp".
  --makeblastdb STRMAKEBLASTDB
                        Provide the path to  makeblastdb. Default call will be to "blastp".
```

## Other options for input to ShortBRED-Identify

ShortBRED-Identify can build markers from other inputs. The options additional options are described below, and are often used when a researcher is rerunning ShortBRED, and wishes to skip redundant parts of the pipeline. 

```
Input:
  --goi SGOIPROTS       Enter the path and name of the proteins of interest file.
  --ref SREFPROTS       Enter the path and name of the file containing reference protein sequences.
  --refdb DIRREFDB      Can be specified in place of reference proteins [--ref]. Enter the path and name for a blastdb of reference proteins.
  --goiblast SGOIBLAST  Used when modifying existing ShortBRED-Identify results. Enter the path and name of the blast results from the goi-to-goi search.
  --refblast SREFBLAST  Used when modifying existing ShortBRED-Identify results. Enter the path and name of the blast results from the goi-to-ref search.
  --goiclust SCLUST     Used when modifying existing ShortBRED-Identify results. Enter the path and name of the clustered genes of interest file.
  --map_in SMAPIN       Used when modifying existing ShortBRED-Identify results. Enter the path and name of the two column file connecting proteins to families.
```

## Method 2: Create ShortBRED markers using a set of proteins of interest, and a blast database of reference proteins.

In our original example, ShortBRED called blastp to build the reference database. If you have already created a reference database using ShortBRED or blastp, you can skip that step and work with it directly. For example, if you ran the first test command for ShortBRED-Identify, you can then enter the command below to use the reference database you built.

```
./shortbred_identify.py --goi  example/input_prots.faa --refdb example_identify/refdb/refdb --markers mytestmarkers_prebuiltdb.faa --tmp example_identify_prebuiltdb
```


## Method 3: Create a new set of ShortBRED markers using a previous set of blast search results.

It may be the case that you have run ShortBRED to create a set of markers, and wish to change a parameter such as marker length. If you are using a large reference database and a large set of protein sequences, it can be time-consuming to run the full search again.

You can rerun ShortBRED by pointing it to the 1) original faa file of clustered proteins of interest (--goiclust), 2) results from proteins of interest vs proteins of interest blast search (--goiblast), 3) results from proteins of interest vs reference proteins blast search (--refblast), and 4) the two column "map" file which maps each original protein of interest to a cluster of proteins.

```
./shortbred_identify.py --goiclust  example_identify/clust/clust.faa --goiblast example_identify/blastresults/selfblast.txt --refblast example_identify/blastresults/refblast.txt --map_in example_identify/clust/clust.map --markers mytestmarkers_previous_search_results.faa
```

Please note that you will not be able to change the cluster id and other parmeters upstream of the blast search using this method.

## Parameter options

Users may alter the parameters ShortBRED uses to build markers.

```
Parameters:
  --clustid DCLUSTID    Enter the identity cutoff for clustering the genes of interest. Examples: .90, .85, .10,...
  --qclustid DQCLUSTID  Enter the identity cutoff for clustering the quasi-markers. Examples: .90, .85, .10,...
  --consthresh DCONSTHRESH
                        Enter the consensus threshold for assigning AA's in the family alignments to the consensus sequences. The default is .70. Examples: .60, .70, .80,...
  --threads ITHREADS    Enter the number of threads to use.
  --id DID              Enter the identity minimum for a short, high-identity region. Examples: .90, .85, .10,...
  --len DL              Enter the length maximum for a short, high-identity region. l=(length hit region)/(length query gene) Examples: .30, .20, .10,... 
  --minAln ILENMIN      Enter the minimum for a short, high-identity region. Examples: 10, 20, 30,... 
  --markerlength IMLENGTH
                        Enter the minimum marker length.
  --totlength ITOTLENGTH
                        Enter the maximum length for the combined markers for a gene. Default is 200
  --qthresh ITHRESH     Enter a maximum quasi-score.
  --qmlength IQMLENGTH  Enter a minimum length for QM's (quasi-markers).
  --xlimit IXLIMIT      Enter the number of Xs to allow in JMs (the number of ambiguous amino acids in junction markers)
```

## Output options

```
Output:
  --markers SMARKERS    Enter name and path for the marker output file
  --map_out SMAP        Enter name and path for the output map file
  --tmpdir STMP         Set directory for temporary output files.
```
---
## Running ShortBRED-Quantify to Profile WGS data and isolate genomes.

## Method 1: Profile a WMS (whole metagenome sample) using ShortBRED-Quantify and a set of markers.

If you would like to test ShortBRED-Quantify using your new markers, enter the following command:

```
$ ./shortbred_quantify.py --markers mytestmarkers.faa --wgs example/wgs.fna  --results exampleresults.txt --tmp example_quantify
```

This command should also run quickly, as there are only 100 nucleotide reads in example/wgs.fna. You can then open up results.txt and see the ShortBRED counts for each protein family, which provides the relative abundance of the protein families in the wgs data. 

---

## Method 2: Profile a set of ORFs from isolate genome (an "annotated" genome) using ShortBRED-Quantify and a set of markers.

ShortBRED can search a file of amino acid ORF's for hits to amino acid markers. The command is similar to running ShortBRED on wgs files, but uses the --genome flag to specify a genome instead of a wgs file.

```
$ python shortbred_quantify.py --genome isolate_ORFs.faa --markers sample_markers/ShortBRED_ABR_101bp_markers.faa --tmp bacterialgenome --usearch /path/to/usearch --maxhits 0 --maxrejects 0
```

There are two additional flags, "--maxhits" and  "--maxrejects" which instruct usearch to allow multiple markers to map to each ORF, and try multiple markers. If you set these to 0, usearch will test the full database of ORFs against each of the markers.


As with ShortBRED-Identify, there are many settings available, which can be viewed using "-h".

---
## Understanding the Output from ShortBRED Quantify

By default, the main results from ShortBRED-Quantify are placed in a file called “results.tab” in the same directory where you ran ShortBRED-Quantify. (You can change this using the "--results" argument.) This lists each family, along with the normalized abundance for the family in RPKM (reads per kilobase of reference sequence per million sample reads).

The following files can be found in the tmp folder for a run of ShortBRED-Quantify:

**(markers).faa.log**

This is a short log file from the ShortBRED-Quantify run, and contains information on the parameters you used for ShortBRED-Quantify along with some basic information on the wms file used. 

**(markers or ORF file).faa.udb**

This is the database created by USEARCH. In the default mode (searching a metagenome file of nucleotide reads), ShortBRED calls USEARCH to build a database of the amino acid markers, and then searches each nucleotide read in the metagenome file against it.

If you are searching a file of ORFs from an isolate genome, ShortBRED will construct the USEARCH database from the amino acid ORFs and search the markers against it.

**(“wgs” or “fullgenome”)_results.tab**

This contains the raw output from the USEARCH search. ShortBRED-Quantify processes this file to identify valid, long length hits and generate the final results.

**SBhits.txt**
 
This file is a subset of *_results.tab described above. It contains only the hits that ShortBRED deemed to be valid given the parameters used in ShortBRED-Quantify. For example, there may be a high identity hit in *_results.tab that will not show up in this file because its length was too low.

**markers.tab**

This presents the results from ShortBRED-Quantify at the level of each individual marker, and contains the following fields:

* **Family** - The family/cluster to which the marker belongs.
* **Marker** - The fasta ID for the marker.
* **Normalized Count** - This is the normalized count for the *marker* described in the section “Profiling protein family metagenomic abundance with ShortBRED-Quantify” of the [ShortBRED paper](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004557). It is expressed in units of RPKMs (reads per kilobase of reference sequence per million sample reads).
* Please note that the “final” result ShortBRED calculates is at the *family* level. It takes the median of these marker normalized counts as the family level value.
* **Hits** - The number of valid ShortBRED hits to the marker.
* **MarkerLength** - The length of the marker in amino acids.
* **ReadLength** - The average read length in the WGS data.
* **HitSpace** - This is the adjusted marker length L’ listed in the [ShortBRED paper](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004557) in the section “Profiling protein family metagenomic abundance with ShortBRED-Quantify”. Please note that this value is in nucleotides, so we convert MarkerLength (in amino acids) to nucleotides by multiplying by 3. For example, using a marker with a MarkerLength of 31 (amino acids) and average WMS read length of 100 (bp), we get: 100-(31*3)-1= 6.

## Contributions ##
Thanks go to these wonderful people:

