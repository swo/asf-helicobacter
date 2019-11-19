Need to make the (forward, reverse, and merged) and (trim/notrim) all separate libraries (i.e., 6 of them!)

Because they need to be deblurred all separately

# Requirements

- `cutadapt`
- `fastq-join`

# To do

- move rdp to `qiime feature-classifier classify-sklearn`
- move closed ref to `qiime vsearch cluster-features-closed-reference`
- move open ref to `qiime vsearch cluster-features-open-reference`
