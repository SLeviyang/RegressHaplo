# The RegressHaplo R Package

This package reconstructs haplotypes from a BAM file using a penalized regression approach.   The algorithm assumes that the BAM file reads are drawn from a low diversity population, roughly <2% diversity.  BAM files reflecting high diversity will lead to long run times and bad performance.

The penaliezd regression has the form

min_x ||y-Px||^2 + rho*J(x)

where J(x) is a quadratic penalty term.   A link with further details regarding the regression will be posted **here** shortly.  

####If you have any questions or encounter a bug, please open an issue through github or e-mail Sivan Leviyang sr286@georgetown.edu.  

## Installation

It's easiest to install RegressHaplo using the R package devtools.  After installing devtools, e.g. `install.packages("devtools")`, execute the following commands in R

```r
library(devtools)
install_github("SLeviyang/RegressHaplo")
```

RegressHaplo depends on the following R packages: **igraph, plyr, dplyr, ggplot2, reshape2, Rsamtools, rmutil, Biostrings, GenomicAlignments**.  These packages must be installed prior to using RegressHaplo.  All are CRAN packages, which can be installed using the `install.packages` command or through the Rstudio GUI, except for Biostrings and GenomicAlignments, which are Bioconductor packages and can be installed as follows:

```r
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
biocLite("GenomicAlignments")
```

## Running RegressHaplo

Once RegressHaplo is installed, you must load it into your R session:

```r
library("RegressHaplo")
```
The RegressHaplo analysis is split into components comprising the following pipeline:


**BAM file --> variant calls --> read\_table --> loci --> potential haplotypes --> parameters --> solutions --> haplotypes --> fasta file**

The intial input to the pipline is a BAM file.  **RegressHaplo assumes that the BAM file has been indexed and that the index file name is the same as bam\_file, but with a .bai appended**.  If your bam file does not have an index, use samtools with the following command

```
samtools index [bam_file]
```

The full pipeline can be run using the command (**the directory out\_dir must exist, RegressHaplo will not create it**):

```r
full_pipeline(bam_file, out_dir, start_pos=NULL, end_pos=NULL, num_trials=700)
```
where 

* `bam_file` - path to the BAM file/
* `out_dir` - directory in which output datafiles are written. 
* `start_pos`,`end_pos` - variation will only be considered within the region between start\_pos and end\_pos, where position is relative to the reference.  All haplotypes returned will have the same length as the reference, but positions outside the specified region will be given consensus values.
* `num_trials` - number of times the penalized regression will be solved using different starting points and different values for rho.  (Values for rho can be specified, see below.)  Raising num\_trials improves reconstruction, but increases run time.


`full_pipeline` has several optional parameters that allow for greater control of the pipeline, see R help through `help("full_pipeline")` for further details.  RegressHaplo will overwrite any files in the output directory during execution of the pipeline. 

**After the pipeline is run, the final haplotype reconstruction is written to the file `final_haplo.fasta` in the specified output directory**

## Parsing RegressHaplo Output

The RegresHaplo haplotype reconstruction can be accessed directly through the `final_haplo.fasta` file.  In the file, each haplotype is of the form **haplotype#\_freq**, e.g. haplotype1\_0.35, where freq gives the reconstructed frequency of the haplotype.

Alternatively, several functions provide an interface to pipeline results:

1. `get_variable_positions.pipeline(out_dir)` - returns the positions on the reference that were assumed variable during the reconstruction, all other positions take on consensus values.
2. `get_fasta.pipeline(out_dir)` - returns a list with two elements: haplotypes and freq.  The element haplotypes is a character.vector containing the reconstructed haplotypes, each equal in length to the reference.  The element freq is a numeric.vector containing the reconstructed frequencies.
3. `get_solutions_summary.pipeline(out_dir)` - returns a data.frame with columns rho (giving the rho value of the particular solution), K (the number of haplotypes reconstruted), and fit (the fit of haplotype reconstruction, lower being better).  

The example dataset, discussed directly below, shows how each of these functions is used.

## Example

The best way to understand RegressHaplo is to look through an example.

See the **example/** folder for an example dataset and output.   See **example/example\_readme.pdf** for the code used to generate the example output and further details.  


## Pipeline Components

Greater control of RegressHaplo can be acheived by running the pipeline components separately.  Each pipeline component takes the output of the previous pipeline components as its input, thereby allowing a user to replace any given pipeline component with their own analysis.  For example, variant calling can be done with another tool.  

### Pipeline Functions

1.  `bam_to_variant_calls.pipeline(bam_file, out_dir)` - outputs **variant\_calls.csv**.
2. `variant_calls_to_read_table.pipeline(bam_file, out_dir)` - outputs **read\_table.csv**.
3. `read_table_to_loci.pipeline(out_dir)` - outputs **loci.csv**.
4. `loci_to_haplotypes.pipeline(out_dir)` - outputs **h.csv**.
5. `haplotypes_to_parameters.pipeline(out_dir)` - outputs **y.csv, P.csv**
6. `parameters_to_solutions.pipeline(out_dir, num_trials=100)` - outputs **solutions.csv**
7. `solutions_to_haplotypes.pipeline(out_dir)` - outputs **final_haplo.csv**
8. `haplotypes_to_fasta.pipeline(bam_file, out_dir)` - outputs **final_haplo.fasta**

**The pipeline component functions have further parameters that allow for greater control of the output.   Use the R help command to see details for each function.**


### Pipeline Output Files

The following files are generated by the RegressHaplo pipeline.  The folder **example/output** contains an example of each of these files.


1. `variant_calls.csv`  A file containing a data.frame with columns: pos, A, C, G, T, -.   pos gives the position number relative to the reference.  The A,C,G,T,- columns have entries of 0 or 1 specifying if a given nucleotide is a true variant (1) or error (0).   **All positions in variant\_calls.csv are taken as variable by RegressHaplo and used for the analysis.  Including non-variable positions will lead to slower run times with no accuracy benefit.**
2. `read_table.csv`  A file containing a data.frame with columns count, pos1, pos2, ..., posn.  posi is the position of the ith variable position specified in `variant_calls.csv`.  Each row of the read\_table data.frame corresponds to a collection of reads that cover the same variable positions and have the same nucleotides at those positions.
3. `loci.csv` A file containing a data.frame with columns locus, pos.  The locus column gives the number of the locus, starting with 1 and incrementing by one.   The pos column gives the subset of the variable positions contained in each locus.  
4. `h.csv` Each row corresponds to a potential haplotype.   Each column has the position as a header.  Columns exactly correspond to variable positions.  
5. `y.csv, P.csv` The y vector and P matrix in the least squares \|y - Ph\|^2.
6. `solutions.csv` - A file containing all solutions generated for the penalized regression.   Each column represents a single solution to the penalized regression.   
7. `final_haplo.csv` - A csv file containing nucleotide values of reconstructed haplotypes, but only at the variable positions.  
8. `final_haplo.fasta` - A fasta file containing the reconstructed haplotypes over the whole reference.  

### Pipeline Accessor Functions

For convenience, the following functions allow the user to access data in files produced by the RegressHaplo pipeline:

1. `get_variant_calls.pipeline(out_dir)` - returns the data.frame written in `variant_calls.csv`.
2. `get_read_table.pipeline(out_dir)` - returns the data.frame written in `read_table.csv`.
3. `get_loci.pipeline(out_dir)` - returns a list of numeric vectors giving the positions, relative to the reference, composing each loci.  
4. `get_y.pipeline(out_dir), get_P.pipeline(out_dir)` - return a numeric vector, the y vector, and a numeric vector, the P matrix.  These functions access `y.csv` and `P.csv` respectively
5. `get_solutions.pipeline(out_dir)` - returns matrix containing solution information.  This function is offered for completeness, analysis of solutions should be done through the solution analysis functions described below.


