DIRECTIONS = ['forward', 'reverse']
SEARCHES = ['exhaustive', 'default']

# note that current directory is mounted as /data in the container
qiime = "docker run -t -i -v $(pwd):/data qiime2/core:2018.2 qiime"

rule all:
    input:
        expand("{direction}-derep-table.tsv", direction=DIRECTIONS),
        #expand("{direction}-derep-usearch-{search}.b6", direction=DIRECTIONS, search=SEARCHES)
        expand("{direction}-derep-usearch-default-taxonomy.tsv", direction=DIRECTIONS),
        expand("{direction}-derep-rdp.tsv", direction=DIRECTIONS)

rule rdp:
    output: "{x}-rdp.tsv"
    input: "{x}-seqs.fa"
    shell: "rdp_classify.py {input} {output}"

rule taxonomy:
    output: "{x}-taxonomy.tsv"
    input: "{x}.b6"
    params: taxonomy="db/97_otus_taxonomy.txt"
    script: "scripts/b6-to-tax.R"

rule usearch_default:
    output: "{x}-usearch-default.b6"
    input: "{x}-seqs.fa"
    params: db="db/97_otus.fasta"
    shell: "usearch -usearch_global {input} -db {params.db} -id 0.97 -strand both -blast6out {output} -maxaccepts 1 -maxrejects 8"

rule usearch_exhaustive:
    output: "{x}-usearch-exhaustive.b6"
    input: "{x}-seqs.fa"
    params: db = "97_otus.fasta"
    shell: "usearch -usearch_global {input} -db {params.db} -id 0.97 -strand both -blast6out {output} -maxaccepts 0 -maxrejects 0"

rule export_tsv:
    output: "{x}-table.tsv"
    input: "{x}-table.biom"
    shell: "biom convert -i {input} -o {output} --to-tsv"

rule export_fasta:
    output: "{x}-seqs.fa"
    input: "{x}-seqs.qza"
    shell:
        qiime + " tools export"
        " {input}"
        " --output-dir ."
        " && mv dna-sequences.fasta {output}"

rule export_table:
    output: "{x}-table.biom"
    input: "{x}-table.qza"
    shell:
        qiime + " tools export"
        " {input}"
        " --output-dir ."
        " && mv feature-table.biom {output}"

rule deblur:
    output:
        table = "{x}-deblur-table.qza",
        seqs = "{x}-deblur-seqs.qza",
        stats = "{x}-deblur-stats.qza"
    input: "{x}-filter.qza"
    shell:
        qiime + " deblur denoise-16S"
        " --i-demultiplexed-seqs {input}"
        " --p-trim-length 253"
        " --p-min-reads 1"
        " --p-jobs-to-start 2"
        " --p-sample-stats"
        " --o-table {output.table}"
        " --o-representative-sequences {output.seqs}"
        " --o-stats {output.stats}"

rule closed_ref:
    output:
        table = "{x}-closedref-table.qza",
        clusters = "{x}-closedref-clustered-seqs.qza",
        unmatched = "{x}-closedref-unmatched-seqs.qza"
    input:
        seqs = "{x}-derep-seqs.qza",
        table = "{x}-derep-table.qza",
        ref = "reference.qza"
    shell:
        qiime + " vsearch cluster-features-closed-reference"
        " --i-sequences {input.seqs}"
        " --i-table {input.table}"
        " --i-reference-sequences {input.ref}"
        " --p-perc-identity 0.97"
        " --p-threads 2"
        " --o-clustered-table {output.table}"
        " --o-clustered-sequences {output.clusters}"
        " --o-unmatched-sequences {output.unmatched}"

rule open_ref:
    output:
        table = "{x}-openref-table.qza",
        clusters = "{x}-openref-clustered-seqs.qza",
        seqs = "{x}-openref-newref-seqs.qza"
    input:
        seqs = "{x}-derep-seqs.qza",
        table = "{x}-derep-table.qza",
        ref = "reference.qza"
    shell:
        qiime + " vsearch cluster-features-open-reference"
        " --i-sequences {input.seqs}"
        " --i-table {input.table}"
        " --i-reference-sequences {input.ref}"
        " --p-perc-identity 0.97"
        " --p-threads 2"
        " --o-clustered-table {output.table}"
        " --o-clustered-sequences {output.clusters}"
        " --o-new-reference-sequences {output.seqs}"

rule reference:
    output: "reference.qza"
    input: "db/97_otus.fasta"
    shell:
        qiime + " tools import"
        " --input-path {input}"
        " --output-path {output}"
        " --type 'FeatureData[Sequence]'"

rule dereplicate:
    output:
        table = "{x}-derep-table.qza",
        seqs = "{x}-derep-seqs.qza"
    input: "{x}-filter.qza"
    shell:
        qiime + " vsearch dereplicate-sequences"
        " --i-sequences {input}"
        " --o-dereplicated-table {output.table}"
        " --o-dereplicated-sequences {output.seqs}"

rule filter:
    output:
        seqs = "{x}-filter.qza",
        stats = "{x}-filter-stats.qza"
    input: "{x}-demux.qza"
    shell:
        qiime + " quality-filter q-score"
        " --i-demux {input}"
        " --o-filtered-sequences {output.seqs}"
        " --o-filter-stats {output.stats}"

rule import_fastq:
    output: "{x}-demux.qza"
    input: "{x}-manifest.txt"
    shell:
        qiime + " tools import"
        " --type 'SampleData[SequencesWithQuality]'"
        " --input-path {input}"
        " --output-path {output}"
        " --source-format SingleEndFastqManifestPhred33"
