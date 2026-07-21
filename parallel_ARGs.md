# Ancestral recombination graph workflow example
Author: Keaka Farleigh, Ph.D.
Date: July 21st, 2026
Email: keakafarleigh@virginia.edu

## Purpose

This workflow will allow us to estimate ancestral recombination graphs and generate summary statistics. We will leverage parallel computing to reduce computational time and also write custom functions to account for non-monophyletic relationships in our data.

Let's load our modules and activate our conda environment. 
```
module load gnu-parallel
conda activate argweaver
```

## Estimate ancestral recombination graphs
We will use [argweaver](https://anaconda.org/channels/genomedk/packages/argweaver/overview) to estimate ancestral recombination graphs. We will perform the analysis in 1 Mb windows across the genome. We get these windows from the `C_pyrrhus_1Mb.bed` file. We then use the while loop to loop through each line of the `C_pyrrhus_1Mb.bed` file to write an argweaver command into the `pyrrhus_b1_con.argweaver.commands.txt` file. The first column of the bed file becomes our `$chrom` variable, the second column becomes our `$start` variable, and the third column becomes our `$end` variable. We then output the results to a directory named out and with the chromosome, start, and end appended to the beginning of our output file (along with a string indicating the species and populations we are comparing). Then, we use [gnu-parallel](https://support.oakland.edu/TDClient/33/Support-Center/KB/PrintArticle?ID=847) to run analyses simultaneously and cut down our computational time.   
```
while read -r line; do chrom=$(echo $line | cut -f 1 -d ' '); start=$(echo $line | cut -f 2 -d ' '); end=$(echo $line | cut -f 3 -d ' '); echo "arg-sample --vcf pyrrhus_b1_con.${chrom}.vcf.gz --region ${chrom}:${start}-${end} -c 10 --maxtime 1e7 -N 1e5 -m 2.4e-9 --unphased -o ./out/pyrrhus_b1_con.${chrom}.${start}-${end}" >> pyrrhus_b1_con.argweaver.commands.txt; done <  ../C_pyrrhus_1Mb.bed

parallel -j 32 :::: pyrrhus_b1_con.argweaver.commands.txt

```

Now lets process the output in preperation for calculating our summary statistics. First, we move the arg files (.smc.gz) to a directory we create for processsing. Then 

```
# Make directory to process output
mkdir process_output

# Change into the directory with the raw ARG output
cd out

# Copy the smc.gz files to the processing directory. Files are generally small, but you could change it to mv if needed.
find ./ -type f -name "*.smc.gz" -exec cp {} ../process_output \; -print

# Change into directory for processing
cd ../process_output

# Make directory to hold the bed files
mkdir argBedFiles 

# Convert to bed files, only use the 1000th iteration for each block
for i in *.1000.smc.gz; do prefix=$(echo $i | sed 's/.1000.smc.gz//g');  echo "smc2bed $i --log-file ../out/${prefix}.log > ./argBedFiles/${prefix}.bed" >> smc2bed.commands.txt; done

# Run the command in parallel
parallel -j 8 :::: smc2bed.commands.txt

# Change into the directory of bedfile
cd argBedFiles

# Combine into one single bed file
cat *.bed > pyrrhus_b1_con_argweaver.bed
```

## Calculate summary statistics 
We will use the `sliding_argstats.R` script that contains our functions and `example_argstats.R` script to execute our analysis. The functions use a bed file that was generated following the `Estimating ancestral recombination graphs` workflow above. contains four columns. The columns are as follows: chromosome, start, end, and ARG of that particular window. We also need a file designating which populations individuals belong to. Note that ARGs are estimated per haplotype, so we will append `_1` and `_2` to individuals to match the argweaver output. Our data is diploid, which is why we add `_1` and `_2`.

```
Rscript example_argstats.R
```
