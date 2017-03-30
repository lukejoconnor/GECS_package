# GE co-score regression code.
This is an initial version of GE co-score regression, as described at <url>http://biorxiv.org/content/early/2017/03/18/118018</url>. It is implemented in MATLAB. Code will be updated to be more user friendly in the future; if you have trouble using it, please contact Luke O'Connor (loconnor@g.harvard.edu) and I will be happy to assist you.

## Required input files
You will need to have MATLAB/Octave installed, and you will need to supply several input files. Required input files:

-GWAS summary statistics, in .sumstats format; see <url>https://github.com/bulik/ldsc/wiki/Summary-Statistics-File-Format</url>. Required order of columns is SNPs, Z scores, allele 1, allele 2, N, anything else.

-eQTL summary statistics, in fastQTL output format. These must be produced for the whole sample and for subsamples of the individuals in your dataset. If there are related individuals, you should prune them. If you have multiple samples from the same individuals, you should make sure that you are partitioning *individuals*, not samples.

-in-sample LD data from each subset of the eQTL cohort, in plink output format. This can be produced using e.g. plink --bfile genotype_filename --R --ld-window 999999 --ld-window-kb 1000 --remove $file_listing_samples_to_exclude --out output_filename

-Allele count data from each subset of the eQTL cohort in plink --hardy output format. This can be produced using e.g. plink --bfile genotype_filename --hardy --remove $file_listing_samples_to_exclude --out output_filename

-LD scores. These can be obtained at <url>https://data.broadinstitute.org/alkesgroup/LDSCORE/</url>, and they should be unzipped.

## Running GECS regression
After generating the required input files, download the repository to your computer or cluster. You need to process your eQTL files before running the main script. eQTL files should be in a single directory and should end with '.qtl'. You also need the allele count information for this step. Run the process_qtls function, passing in the path to the directory containing .qtl files, the name of the file containing allele counts, and the output file name. You need to do this three times: for the entire eQTL sample, and for each eQTL cohort subsample. You should name these output files some_path/some_filename.0.mat, some_path/some_filename.1.mat, and some_path/some_filename.2.mat.
, and open main_script.m in a text editor. 

Next, open main_scrip.m. Input the paths to your files (some instructions are contained within the script). Run the script. It will produce a .mat file containing estimates of rho^2 (and other things) for each phenotype. 
