Need to make the (forward, reverse, and merged) and (trim/notrim) all separate libraries (i.e., 6 of them!)

Because they need to be deblurred all separately

# Requirements

- `cutadapt`
- `fastq-join`

# To do

## Errors

- Trimming and excluding files means that merging doesn't work?
- Instead, merge-notrim, and then trim all the files inside of Q2 using cutadapt
    - `--p-discard-untrimmed`
- Did I try to do this before and, for some reason, it didn't work?

## Later

- Open reference has some sequences that get clustered *de novo*
- In most cases, there is only one of those that is very big
- I look in the tax table for the big, unknown-taxonomy OTUs
- Then I export the open ref newref seqs and pull out that sequence
- Then I usearch for that sequence in 97 OTUs; needing to lower the threshold
