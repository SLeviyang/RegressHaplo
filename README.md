
# RegressHaplo R Package

This package infers haplotypes from a BAM file using a penalized regression approach. To install the package:

....

## RegressHaplo Pipeline

The package is split into components comprising the following pipeline:


**BAM file --> variant calls --> read\_table --> loci --> potential haplotypes --> parameters --> solutions --> solution haplotypes**

The function calls implementing each step of the pipeline are as follows. These calls should be executed sequentially.   

1.  `bam_to_variant_calls.pipeline(bam_file, out_dir, sig=.01)` 
2. `bam_to_read_table.pipeline(bam_file, out_dir)`
3. `read_table_to_loci.pipeline(out_dir)`
4. `loci_to_haplotypes.pipeline(out_dir)`
5. `haplotypes_to_parameters.pipeline(out_dir)`
6. `parameters_to_solutions.pipeline(out_dir, num_trials=100)`
7. `solutions_to_haplotypes.pipeline(out_dir)`

**Arguments**:

* `bam_file` - A BAM file of the NGS dataset.  The reference is not needed by RegressHaplo. 
* `out_dir` - directory in which output and input datafiles are written and read.
* `sig` - significance level at which to call variants.  A Bonferonni correction is automatically applied to the significance level accounting for multiplicity of positions.
* `num_trials` - number of times to solve the penalized regression in order to search for an optimal solution.  

**Generated Files**:

At each step in the pipeline, a file is produced in `out_dir`.  All input files are also expected to be in `out_dir`, except for the bam_file since its path is given as an argument in the first two steps of the pipeline (subsequent steps don't need the bam file).  For example, the function `bam_to_read_table.pipeline` expects `variant_calls.csv` to be in `out_dir` and writes `read_table.csv` to `out_dir`.

The pipeline format allows the user to edit the files generated at either step or to create their own files, thereby bypassing certain steps in the pipeline

1. `variant_calls.csv`  A file containing a data.frame with columns: pos, A, C, G, T, -.   pos gives the position number relative to the reference (which is not needed by RegressHaplo).  The A,C,G,T,- columns have entries of 0 or 1 specifying if a given nucleotide is a true variant (1) or error (0). 
2. `read_table.csv`  A file containing a data.frame with columns count, pos1, pos2, ..., posn.  posi is a number (e.g. the column pos1 will actually be names "56") giving the position of the ith variant call specified in `variant_calls.csv`.  Each row of the data.frame corresponds to a collection of reads that cover the same variable positions and have the same nucleotides at those positions....
3. `loci.csv` A file containing a data.frame with columns locus, pos.  The locus column gives the number of the locus, starting with 1 and incrementing by one.   The pos column gives the subset of variable positions (specified in `variant_calls.csv` contained in each locus.  The format of each entry is pos1+pos2+..+posn.  For example, 100+200+220 as an entry in the pos column represents a locus containing positions 100,200,220 relative to the reference.
4. `h.csv` Each row corresponds to a potential haplotype.   Each column has the position as a header.  Columns exactly correspond to variable positions.  
3. `y.csv, P.csv` The y vector and P matrix in the least squares \|y - Ph\|^2.
4. `solutions.csv`
5. `haplo.csv`

The final haplotype reconstruction is contained in HAPLO.csv.  

**Running the Whole Pipeline**:

The full pipleline can be executed through `full_pipeline(bam_file, out_dir, variant_calls, num_trials=100)`