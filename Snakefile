import glob, os.path

LIBS = set([re.search("notrim/(.*)_R1.fastq", x).group(1) for x in glob.glob("notrim/*_R1.fastq")])

ALL_FILES = expand("{folder}/{lib}_R{direct}.fastq", folder=["notrim", "trim"], lib=LIBS, direct=["1", "2"]) + \
    expand("{folder}/{lib}_R1.fastq", folder=["merge-notrim", "merge-trim"], lib=LIBS)

# DIRECTIONS = ['forward', 'reverse', "merge"]
# TRIMS = ["notrim", "trim"]
# SEARCHES = ['exhaustive', 'default']

# note that current directory is mounted as /data in the container
qiime = "docker run -t -i -v $(pwd):/data qiime2/core:2018.11 qiime"

configfile: "config.yaml"

rule all:
    input:
        "manifest.txt"
        # expand("{direction}-{trim}-derep-table.tsv", direction=DIRECTIONS, trim=TRIMS),
        #expand("{direction}-derep-usearch-{search}.b6", direction=DIRECTIONS, search=SEARCHES)
        # expand("{direction}-derep-usearch-default-taxonomy.tsv", direction=DIRECTIONS),
        # expand("{direction}-derep-rdp.tsv", direction=DIRECTIONS)

rule clean:
    shell: "rm trim/* merge-trim/* merge-notrim/* manifest.txt"

rule write_manifest:
    output: "manifest.txt"
    input: ALL_FILES
    script: "write_manifest.py"

rule merge:
    output: "merge-{x}/{y}_R1.fastq"
    input: forward="{x}/{y}_R1.fastq", reverse="{x}/{y}_R2.fastq"
    script: "merge.py"

rule trim_forward:
    output: "trim/{x}_R1.fastq"
    input: "notrim/{x}_R1.fastq"
    params: primer=config["forward-primer"]
    script: "trim.py"

rule trim_reverse:
    output: "trim/{x}_R2.fastq"
    input: "notrim/{x}_R2.fastq"
    params: primer=config["reverse-primer"]
    script: "trim.py"

#######################################################################

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
        " --input-path {input}"
        " --output-path ./"
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

rule filter_pair:
    output:
        seqs = "{x}-pair-filter.qza",
        stats = "{x}-pair-filter-stats.qza"
    input: "{x}-pair-demux.qza"
    shell:
        qiime + " quality-filter q-score-joined"
        " --i-demux {input}"
        " --o-filtered-sequences {output.seqs}"
        " --o-filter-stats {output.stats}"

rule filter_single:
    output:
        seqs = "{x}-single-filter.qza",
        stats = "{x}-single-filter-stats.qza"
    input: "{x}-single-demux.qza"
    shell:
        qiime + " quality-filter q-score"
        " --i-demux {input}"
        " --o-filtered-sequences {output.seqs}"
        " --o-filter-stats {output.stats}"

rule trim_pair_reverse:
    output: "merge-trim-demux.qza"
    input: "merge-trimforwardonly-demux.qza"
    params: primer=config["reverse-primer"]
    shell:
        qiime + " cutadapt trim-single"
        " --i-demultiplexed-sequences {input}"
        " --p-adapter {params.primer}"
        " --o-trimmed-sequences {output}"

rule trim_pair_forward:
    output: "merge-trimforwardonly-demux.qza"
    input: "merge-notrim-demux.qza"
    params: primer=config["forward-primer"]
    shell:
        qiime + " cutadapt trim-single"
        " --i-demultiplexed-sequences {input}"
        " --p-front {params.primer}"
        " --o-trimmed-sequences {output}"

rule import:
    output: "{x}-notrim-demux.qza"
    input: "{x}-manifest.txt"
    shell:
        qiime + " tools import"
        " --type 'SampleData[SequencesWithQuality]'"
        " --input-path {input}"
        " --output-path {output}"
        " --input-format SingleEndFastqManifestPhred33"

rule merge_manifest:
    input: forward="forward-manifest.txt", reverse="reverse-manifest.txt"
    output: "merge-manifest.txt"
    script: "merge-pairs.py"
