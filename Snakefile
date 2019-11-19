import glob, os.path

DIRECTIONS = ["forward", "reverse", "merge"]
TRIMS = ["notrim", "trim"]
GROUPS = expand("{direction}-{trim}", direction=DIRECTIONS, trim=TRIMS)

# library names from raw, input files
LIBRARIES = set([re.search("raw/(.*)_R1.fastq", x).group(1) for x in glob.glob("raw/*_R1.fastq")])

# groups of files expected as outputs
MANIFEST_FILES = expand("{group}/{library}_R1.fastq", group=GROUPS, library=LIBRARIES)

ANALYSES = ["deblur-table.tsv", "deblur-rdp.tsv", "derep-usearch-default-taxonomy.tsv", "derep-usearch-exhaustive-taxonomy.tsv"]

# note that current directory is mounted as /data in the container
qiime = "docker run -t -i -v $(pwd):/data qiime2/core:2019.10 qiime"

configfile: "config.yaml"

rule all:
    input:
        # expand("{group}-{analysis}", group=GROUPS, analysis=ANALYSES)
        "forward-notrim-deblur-rdp-tax.qza"
        , "forward-notrim-derep-seqs.qza"
        , "forward-notrim-derep-table.tsv"
        , "forward-notrim-closedref-table.tsv",

rule clean:
    shell:
        "rm -f"
        " *.txt *.qza *.biom *.tsv *.fa *.b6 *.log" + \
        " ".join(expand("{group}/*", group=GROUPS))

rule export_taxonomy:
    output: "{x}-tax.tsv"
    input: "{x}-tax.qza"
    shell:
        qiime + " tools export"
        " --i-input-path {input}"
        " --output-path ."
        " && mv taxonomy.tsv {output}"

rule rdp:
    output: "{x}-rdp-tax.qza"
    input: reads="{x}-seqs.qza", classifier="rdp-classifier.qza"
    shell:
        qiime + " feature-classifier classify-sklearn"
        " --i-classifier {input.classifier}"
        " --i-reads {input.reads}"
        " --o-classification {output}"

rule download_classifier:
    output: "rdp-classifier.qza"
    params: url=config["classifier-url"], md5=config["classifier-md5"]
    script: "scripts/download-and-check-md5.py"

# rule taxonomy:
#     output: "{x}-taxonomy.tsv"
#     input: "{x}.b6"
#     params: taxonomy="db/97_otu_taxonomy.txt"
#     script: "scripts/b6-to-tax.R"

# rule usearch_default:
#     output: "{x}-usearch-default.b6"
#     input: "{x}-seqs.fa"
#     params: db="db/97_otus.udb"
#     shell: "usearch -usearch_global {input} -db {params.db} -id 0.97 -strand both -blast6out {output} -maxaccepts 1 -maxrejects 8"

# rule usearch_exhaustive:
#     output: "{x}-usearch-exhaustive.b6"
#     input: "{x}-seqs.fa"
#     params: db="db/97_otus.udb"
#     shell: "usearch -usearch_global {input} -db {params.db} -id 0.97 -strand both -blast6out {output} -maxaccepts 0 -maxrejects 0"

rule deblur:
    output:
        table = "{x}-deblur-table.qza",
        seqs = "{x}-deblur-seqs.qza",
        stats = "{x}-deblur-stats.qza"
    input: "{x}-filter.qza"
    shell:
        qiime + " deblur denoise-16S"
        " --i-demultiplexed-seqs {input}"
        " --p-trim-length 100"
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
        ref = "97_otus.qza"
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
        ref = "97_otus.qza"
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

rule reference_taxonomy2:
    output: "97_otus_taxonomy.qza"
    input: "97_otus_taxonomy.tsv"
    shell:
        qiime + " tools import"
        " --input-path {input}"
        " --output-path {output}"
        " --type 'FeatureTable[Taxonomy]'"
        " --input-format TSVTaxonomyFormat"

rule reference_taxonomy1:
    output: "97_otus_taxonomy.tsv"
    input: "db/97_otu_taxonomy.txt"
    script: "scripts/convert-taxonomy.py"

rule reference_otus:
    output: "97_otus.qza"
    input: "db/97_otus.fasta"
    shell:
        qiime + " tools import"
        " --input-path {input}"
        " --output-path {output}"
        " --type 'FeatureData[Sequence]'"

rule export_tsv:
    output: "{x}-table.tsv"
    input: "{x}-table.biom"
    shell: "biom convert -i {input} -o {output} --to-tsv"

rule export_fasta:
    output: "{x}-seqs.fa"
    input: "{x}-seqs.qza"
    shell:
        qiime + " tools export"
        " --input-path {input}"
        " --output-path ."
        " && mv dna-sequences.fasta {output}"

rule export_table:
    output: "{x}-table.biom"
    input: "{x}-table.qza"
    shell:
        qiime + " tools export"
        " --input-path {input}"
        " --output-path ./"
        " && mv feature-table.biom {output}"

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
    output: seqs="{x}-filter.qza", stats="{x}-filter-stats.qza"
    input: "{x}-demux.qza"
    shell:
        qiime + " quality-filter q-score"
        " --i-demux {input}"
        " --o-filtered-sequences {output.seqs}"
        " --o-filter-stats {output.stats}"

rule demux:
    output: "{x}-demux.qza"
    input: "{x}-manifest.txt"
    shell:
        qiime + " tools import"
        " --type 'SampleData[SequencesWithQuality]'"
        " --input-path {input}"
        " --output-path {output}"
        " --input-format SingleEndFastqManifestPhred33"

rule write_manifest:
    output: expand("{group}-manifest.txt", group=GROUPS)
    input: MANIFEST_FILES
    params: libraries=LIBRARIES
    script: "scripts/write_manifest.py"

rule merge:
    output: "merge-{x}/{y}.fastq"
    input: forward="forward-{x}/{y}.fastq", reverse="reverse-{x}/{y}.fastq"
    script: "scripts/merge.py"

rule trim_forward:
    output: "forward-trim/{x}.fastq"
    input: "forward-notrim/{x}.fastq"
    params: primer=config["forward-primer"]
    script: "scripts/trim.py"

rule trim_reverse:
    output: "reverse-trim/{x}.fastq"
    input: "reverse-notrim/{x}.fastq"
    params: primer=config["reverse-primer"]
    script: "scripts/trim.py"

rule copy_reverse:
    output: "reverse-notrim/{x}_R1.fastq"
    input: "raw/{x}_R2.fastq"
    shell: "cp {input} {output}"

rule copy_forward:
    output: "forward-notrim/{x}_R1.fastq"
    input: "raw/{x}_R1.fastq"
    shell: "cp {input} {output}"
