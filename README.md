# The RegressHaplo R Package

This package reconstructs haplotypes from a BAM file using a penalized regression approach.   The algorithm assumes that the BAM file reads are drawn from a low diversity population, roughly <2% diversity.  BAM files reflecting high diversity will lead to long run times and bad performance.

The penalized regression has the form

min_x ||y-Px||^2 + rho*J(x)

where J(x) is a quadratic penalty term.   The RegressHaplo algorithm is discussed in detail in the publication:   *Leviyang S., Griva I., Ita S, and Johnson W. A Penalized Regression Approach to Haplotype Reconstruction of Viral Populations Arising in Early HIV/SIV Infection. Bioinformatics (2017)*.

### If you have any questions or encounter a bug, please open an issue through github or e-mail Sivan Leviyang sr286@georgetown.edu.  This package can be modified and distributed under the Apache 2.0 License, see LICENSE file for details.  


## Installation


RegressHaplo depends on the following R packages: **igraph, plyr, dplyr, ggplot2, rmutil, Rsamtools, Biostrings, GenomicAlignments**.  These packages must be installed prior to using RegressHaplo.  All are CRAN packages, which can be installed using the `install.packages` command or through the Rstudio GUI, except for Biostrings, Rsamtools, and GenomicAlignments, which are Bioconductor packages and can be installed as follows:

```r
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
biocLite("Rsamtools")
biocLite("GenomicAlignments")
```

It's easiest to install RegressHaplo using the R package devtools.  After installing devtools, e.g. `install.packages("devtools")`, execute the following commands in R

```r
library(devtools)
install_github("SLeviyang/RegressHaplo")
```


## Quick Start: Running RegressHaplo

Once RegressHaplo is installed, you must load it into your R session:

```r
library("RegressHaplo")
```
The RegressHaplo analysis is split into components comprising the following pipeline:


**BAM file --> variant calls and pile up --> read\_table --> loci --> potential haplotypes --> parameters --> solutions --> haplotypes --> fasta file**

The intial input to the pipline is a BAM file.  **RegressHaplo assumes that the BAM file has been indexed and that the index file name is the same as bam\_file, but with a .bai appended**.  If your bam file does not have an index, use samtools with the following command

```
samtools index [bam_file]
```

The full pipeline can be run using the following command (**the directory out\_dir must exist, RegressHaplo will not create it**):

```r
full_pipeline(bam_file, out_dir, start_pos=NULL, end_pos=NULL, num_trials=700)
```
where 

* `bam_file` - path to the BAM file/
* `out_dir` - directory in which output datafiles are written. 
* `start_pos`,`end_pos` - variation will only be considered within the region between start\_pos and end\_pos, where position is relative to the reference.  All haplotypes returned will have the same length as the reference, but positions outside the specified region will be given consensus values.
* `num_trials` - number of times the penalized regression will be solved using different starting points and different values for rho.  Raising num\_trials improves reconstruction, but increases run time.


`full_pipeline` has several optional parameters that allow for greater control of the pipeline, use R help command, i.e. `help("full_pipeline")`, for further details.  RegressHaplo will overwrite any files in the output directory during execution of the pipeline. 

###After the pipeline is run, the final haplotype reconstruction is written to the file `final_haplo.fasta` in the specified output directory


The RegresHaplo haplotype reconstruction can be accessed directly through the `final_haplo.fasta` file.  In the file, each haplotype is named with the form **haplotype#\_freq**, e.g. haplotype1\_0.35, where freq gives the reconstructed frequency of the haplotype, eg. 0.35.

Alternatively, several functions provide an interface to pipeline results:

1. `get_variable_positions.pipeline(out_dir)` - returns the positions on the reference that were assumed variable during the reconstruction, all other positions take on consensus values.
2. `get_fasta.pipeline(out_dir)` - returns a list with two elements: haplotypes and freq.  The element haplotypes is a character.vector containing the reconstructed haplotypes, each equal in length to the reference.  The element freq is a numeric.vector containing the reconstructed frequencies.
3. `get_solutions_summary.pipeline(out_dir)` - returns a data.frame with columns rho (giving the rho value of the particular solution), K (the number of haplotypes reconstruted), and fit (the fit of haplotype reconstruction, lower being better).  

The example dataset, discussed directly below, shows how each of these functions is used.

## Quick Start Example

The **example/** folder has an example dataset (in data/) and example RegressHaplo output (in output/).  To generate your own version of the output, change your R working directory to the example/ directory, and then execute the following commands, which will create RegressHaplo output in the directory my\_output.  The reconstruction will be written to the file **my\_output/final_haplo.fasta**

```{r}
# load RegressHaplo
library("RegressHaplo")

# BAM file containing the dataset
bam_file <- "data/example.bam"

# directory in which RegressHaplo output will be written
out_dir <- "my_output/"
dir.create(out_dir)

# run the RegressHaplo pipeline
full_pipeline(bam_file, out_dir, start_pos=500, end_pos=1500, num_trials=700)
```

###See **example/example_readme.pdf** for further details on using the RegressHaplo pipeline.  
The code we used to produce the results in \output and further details regarding using the pipeline are contained in **example/example_readme.pdf** which was created using the R markdown file example/example_readme.Rmd.   


## The RegressHaplo Pipeline

Greater control of RegressHaplo can be acheived by running the pipeline components separately.  Each pipeline component takes the output of the previous pipeline components as its input, thereby allowing a user to alter parameters in a given pipeline component or replace any given pipeline component with their own analysis.  For example, variant calling can be modified by changing the significance level at which variants are called (*through the `sig` parameter in `variant_calls_to_read_table.pipeline`) or variant calling can be done with another tool.  

### Pipeline Functions

1.  `bam_to_variant_calls.pipeline(bam_file, out_dir,
                                          start_pos=NULL, end_pos=NULL,
                                          sig=.01, heavy_tail=T)` - outputs **variant\_calls.csv** and **pile\_up.csv**.
2. `variant_calls_to_read_table.pipeline(bam_file,
                                       out_dir,
                                       use_raw_read_table=F,
                                       sig=.01)` - outputs **raw\_read\_table.csv** and **read\_table.csv**.
3. `read_table_to_loci.pipeline(out_dir, max_num_haplotypes=1200)` - outputs **loci.csv**.
4. `loci_to_haplotypes.pipeline(out_dir,
                                        max_num_haplotypes=1000)` - outputs **h.csv**.
5. `haplotypes_to_parameters.pipeline(out_dir)` - outputs **y.csv, P.csv**
6. `parameters_to_solutions.pipeline(out_dir, num_trials=700,
                                             rho_vals=c(.1,1,5,10,20))` - outputs **solutions.csv**
7. `solutions_to_haplotypes.pipeline(out_dir, K_val=NULL)` - outputs **final_haplo.csv**
8. `haplotypes_to_fasta.pipeline(bam_file, out_dir)` - outputs **final_haplo.fasta**

**A few of the parameters to the pipeline functions are discussed below in the section. 'Examples of Controlling the Pipeline', use the R help command for a particualr function to see details for other parameters.**


### Pipeline Output Files

The following files are generated by the RegressHaplo pipeline.  The folder **example/output** contains an example of each of these files.


1. `variant_calls.csv`  A file containing a data.frame with columns: pos, A, C, G, T, -, i.   pos gives the position number relative to the reference.  The A,C,G,T,-,i columns have entries of 0 or 1 specifying if a given nucleotide or indel is a true variant (1) or error (0).   
2. `pile_up.csv`   A file containing a data.frame with columns pos, A, C, G, T, d, i, coverage.   The pos column gives the position on the reference, and all other columns give the frequency of the nucleotide and indels at the position over all reads covering the position.  **In contrast to `variant_calls.csv`, all positions on the reference are included, regardless of variability.**
3. `raw_read_table.csv`  A file containing a data.frame with columns count, pos1, pos2, ..., posn.  posi is the position of the ith variable position specified in `variant_calls.csv`.  Each row of the data.frame corresponds to a collection of reads that cover the same variable positions and have the same nucleotides at those positions.
4. `read_table.csv`  A file with the same format as `raw_read_table.csv`, but with error correction so that some read are removed.   Typically `read_table.csv` will be a much smaller file than `raw_read_table.csv`.   All subsequent pipeline components use `read_table.csv` for analysis and do not depend on `raw_read_table.csv`.
5. `loci.csv` A file containing a data.frame with columns locus, pos.  The locus column gives the number of the locus, starting with 1 and incrementing by one.   The pos column gives the subset of the variable positions contained in each locus.  
6. `h.csv` Each row corresponds to a potential haplotype.   Each column has the position as a header.  Columns exactly correspond to variable positions.  
7. `y.csv, P.csv` The y vector and P matrix in the least squares \|y - Ph\|^2.
8. `solutions.csv` - A file containing all solutions generated for the penalized regression.   Each column represents a single solution to the penalized regression.   
9. `final_haplo.csv` - A csv file containing nucleotide values of reconstructed haplotypes, but only at the variable positions.  
10. `final_haplo.fasta` - A fasta file containing the reconstructed haplotypes over the whole reference.  

### Pipeline Accessor Functions

For convenience, the following functions allow the user to access data in files produced by the RegressHaplo pipeline:

1. `get_variant_calls.pipeline(out_dir)` - returns the data.frame written in `variant_calls.csv`.
2. `get_pile_up.pipeline(out_dir)` - returns the data.frame written in `pile_up.csv`
3. `get_read_table.pipeline(out_dir)` - returns the data.frame written in `read_table.csv`.
4. `get_loci.pipeline(out_dir)` - returns a list of numeric vectors giving the positions, relative to the reference, composing each loci.  
5. `get_y.pipeline(out_dir), get_P.pipeline(out_dir)` - return a numeric vector, the y vector, and a numeric vector, the P matrix.  These functions access `y.csv` and `P.csv` respectively
6. `get_solutions.pipeline(out_dir)` - returns matrix containing solution information.  This function is offered for completeness, analysis of solutions should be done through the solution analysis functions described below.

## Evaluation of Solutions

Two functions allow the user to examine the solutions generated by the pipeline.  See **example/example_readme.pdf** and the help documentation for further details on using these functions.

1. `get_solutions_summary.pipeline(out_dir)` - returns a data.frame containing the fit, number of haplotypes, and penalty parameter for each solution.
2. `solution_accuracy.pipeline(out_dir, 
                                  solution, 
                                  nucs=c("A", "C"),
                                  plot=T)` - given a solution (through the solution number), the predicted frequencies and sampled frequency at each variable position are plotted and a corresponding data.frame returned.

## Examples of Controlling the Pipeline

### **Controlling variant calling**:

* Relevant file: `variant_calls.csv` 
* Relevant function: `bam_to_variant_calls.pipeline`:

Variant calling attempts to find errors at given positions:  for example, if at pos 100, the A nucleotide represents an error, then variant calling will remove all A's from the reads at pos 100 and replace them with the consensus nucleotide.  All positions in variant\_calls.csv are taken as variable by RegressHaplo and used for the analysis.  Including non-variable positions will lead to slower run times with no accuracy benefit.  If variability along the reference is very high, RegressHaplo performance may be poor because there are too many potential haplotypes.

With this in mind, sometimes, even if a position is variable, it is useful to remove it from consideration; for example, when the user is not interested in particular regions of the reference or when including all variable positions leads to poor results. To prevent RegressHaplo from considering a given position as variable, remove the corresponding row from the `variant_calls.csv` file.   Alternatively, the `sig` parameter to `bam_to_variant_calls.pipeline` can be set lower, leading to less variant calls, i.e. less rows in `variant_calls.csv`.

### **Controlling read table error correction**:

* Relevant files: `raw_read_table.csv`, `read_table.csv` 
* Relevant function:`variant_calls_to_read_table.pipeline`:

Variant calling cannot address errors that transform true variable nucleotides into other true variable nucleotides.  For example, suppose at position 100, there are true variants A and G, while at position 101 there are also true variants A and G.  Further, suppose that those variants occur only in the pairs AA (for positions 100,101) and GG.   Then there will be some errors that convert, say, the A at position 100 to a G.  These errors will not be caught be variant calling, but will create false GA reads.  The problem is particularly serious given the high coverage and relatively long read length (i.e. 200 base pair) in typical NGS dataset.  Failure to account for such errors often leads to an enourmous number of potential haplotypes and, consequently, poor RegressHaplo performance.

The file `raw_read_table.csv` contains all reads in the dataset.  A read error correction step is then performed, and reads that pass the error correction step are contained in `read_table.csv`.  Reads in `raw_read_table.csv` that do not pass the error correction step are thrown out.   The `sig` parameter to the function `variant_calls_to_read_table.pipeline` can be used to modulate the number of reads that are thrown out by the error correction (i.e. lowering `sig` leads to more reads being thrown out).   Choosing the right `sig` values amounts to a balance between overfitting and underfitting.  If too many reads are thrown out, RegressHaplo will run quickly but not be able to reconstruct the true haploytpes (underfitting).  If too few reads are thrown out, RegressHaplo will run slowly and likely reconstruct false haplotypes (overfitting).

The most time intensive step of `variant_calls_to_read_table.pipeline` is typically construction of `raw_read_table.csv`.  Once `raw_read_table.csv` is constructed, users can set the `use_raw_read_table=T` to use the existing `raw_read_table.csv` and just create a `read_table.csv` file through error correction.   This can be useful when interatively attempting to find the right `sig` value for read error correction.

### **Controlling locus construction**:

* Relevant files: `loci.csv`
* Relevant function:`read_table_to_loci.pipeline`

RegressHaplo groups the variable positions (as specified in `variant_calls.csv`) into loci.  On each locus, RegressHaplo attempts to construct all possible haplotypes reflecting the reads that cover the locus.  At first, RegressHaplo builds loci by selecting regions of the reference that have no linkage information across them.  Once this split is made, there are often too many haplotypes for RegressHaplo to consider in reasonable time, so loci are further split.   The parameter `max_num_haplotype` sets the maximum number of haplotypes RegressHaplo allows; reasonable run times on the order of several hours are acheived for `max_num_haplotype` in the 500-1200 range.  Higher values are not recommended and run times are uncertain.

The file `loci.csv` specifies the loci.  The file has two columns:  locus and pos.  The locus column just labels the locus number, i.e. 1-n where n is the total number of loci.  The pos column specifies the variable postions that are considered part of a single locus.  For example, an entry such as '1+10+23+40' represents a locus covering variable postions 1,10,23, and 40.  

Users can manually modify the `loci.csv` file to create specific loci of interest.  The positions included must exactly match the positions included in the `read_table.csv` file, with each position included in exactly one locus.  Note: Some positions in `variant_calls.csv` may not be included in `read_table.csv` because none of the reads covering those positions passed the read error correction step.  

### **Controlling potential haplotype construction**:

* Relevant files: `h.csv`
* Relevant function:`loci_to_haplotypes.pipeline`

Once loci are determined, the function `loci_to_haplotypes.pipeline` constructs global haplotypes, covering the full reference, by combining all possible combinations of locus haplotypes.  If the number of global haplotypes exceeds the value of `max_num_haplotypes` passed to `loci_to_haplotypes.pipeline`, then RegressHaplo does local regressions on each locus and selects a subset of the locus halpotypes.  This process is repeated until the number of global haplotypes is less than `max_num_haplotypes`.

Generally, constructing loci and choosing haplotype subsets is complex and the choices RegressHaplo makes may not be optimal, or the choices may not reflect user goals.  **Users can essentially skip the entire pipeline up to and including `loci_to_haplotypes.pipeline` by specifying their own haplotypes through the `h.csv` file**.  The columns of `h.csv` are the position numbers of the variable positions (i.e. postions in `variant_calls.csv`).  Each row gives the nucleotide and indel pattern of a global haplotype across the variable positions.

As an example application, users can use another haplotype reconstruction tool to construct haplotypes on two disjoint pieces of the reference, say locus 1 and locus 2.   The constructed haplotypes can be used to produce global haplotypes by pasting together all possible combinations of haplotypes on locus 1 and 2.   RegressHaplo can then be used to select the best global haplotypes, thereby inferring the linkage information that another tool is not able to infer.


