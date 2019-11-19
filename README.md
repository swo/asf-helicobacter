Need to make the (forward, reverse, and merged) and (trim/notrim) all separate libraries (i.e., 6 of them!)

Because they need to be deblurred all separately

# Requirements

- `cutadapt`
- `fastq-join`

I use `cutadapt` rather than Qiime 2's because I found that the option to
keep only trimmed sequences does not work in Qiime 2 version 2019.10.

I use `fastq-join` rather than Qiime 2's merging plugin, because that plugin
uses vsearch, which applies too strict a merging criteria to get most of the
reads merged, and those criteria cannot be relaxed using parameter options.

# To do

## Later

- Open reference has some sequences that get clustered *de novo*
- In most cases, there is only one of those that is very big
- I look in the tax table for the big, unknown-taxonomy OTUs
- Then I export the open ref newref seqs and pull out that sequence
- Then I usearch for that sequence in 97 OTUs; needing to lower the threshold
