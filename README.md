# ASF-*Helicobacter*

## To do

## Overview of files

There are 3 inputs: forward, reverse, merged.

Each of these 3 are cleaned:

- Quality filtering
- Deblur

Then, OTUs are picked with one of 4 methods (all using Greengenes 13\_8):

- Qiime 1 closed reference
- Qiime 1 open reference
- Qiime 2 closed reference
- Qiime 2 open

## Outputs

Aside from the regular taxa plots, I also looked at where the *de novo*
clustered OTUs in the open-reference pickings is classified by `usearch` on the
Greengenes database, using the "default" Qiime 1 search (1 accepts and 8
rejects) versus an "exhaustive" search (unlimited accepts and rejects).

## Requirements

- `cutadapt`
- `fastq-join`

I use `cutadapt` rather than Qiime 2's because I found that the option to
keep only trimmed sequences does not work in Qiime 2 version 2019.10.

I use `fastq-join` rather than Qiime 2's merging plugin, because that plugin
uses vsearch, which applies too strict a merging criteria to get most of the
reads merged, and those criteria cannot be relaxed using parameter options.
