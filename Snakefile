import uuid

configfile: "config.yaml"

DIRECTIONS = ["forward", "reverse", "merge"]
PICKS = ["q2closed", "q2open", "q1closed", "q1open"]

# library names from raw, input files
LIBRARIES = config["libraries"].keys()

# note that current directory is mounted as /data in the container
qiime = "docker run -t -i -v $(pwd):/data qiime2/core:2019.10 qiime"

rule all:
    input: expand("{direction}-{pick}-{product}", pick=PICKS, direction=DIRECTIONS, product=["taxtable.tsv"])

# rule clean:
#     shell:
#         "rm -rf"
#         " *.txt *.qza *.biom *.tsv *.fasta *.b6 *.log *.udb *.tar.gz *.pdf" + \
#         " " + " ".join(expand("{group}/*", group=GROUPS))

rule taxa_plot:
    output:
        plot="taxaplot.pdf",
        data="taxaplot.tsv"
    input: expand("{direction}-{pick}-taxtable.tsv", direction=DIRECTIONS, pick=PICKS)
    script: "scripts/taxa-plot.R"

# rule usearch_openref:
#     output: "{x}-openref-usearch.tsv"
#     input:
#         table="{x}-openref-table.tsv",
#         seqs="{x}-openref-newref-seqs.fasta",
#         otus="97_otus.udb",
#         taxonomy="97-otu-tax.tsv"
#     script: "scripts/openref-usearch.R"

rule export_fasta:
    output: "{x}-seqs.fasta"
    input: "{x}-seqs.qza"
    params: tempdir=lambda wc: wc["x"] + "-" + uuid.uuid4().hex
    shell:
        qiime + " tools export"
        " --input-path {input}"
        " --output-path {params.tempdir}"
        " && mv {params.tempdir}/dna-sequences.fasta {output}"
        " && rmdir {params.tempdir}"

rule export_table:
    """Export to temporary directory to prevent race"""
    output: "{x}-table.tsv"
    input: "{x}-table.qza"
    params: tempdir=lambda wc: wc["x"] + "-" + uuid.uuid4().hex
    shell:
        qiime + " tools export"
        " --input-path {input}"
        " --output-path tmp/{params.tempdir}"
        " && biom convert -i tmp/{params.tempdir}/feature-table.biom -o {output} --to-tsv"
        " && rm -rf tmp/{params.tempdir}"


# Make Qiime 2 taxonomy tables -----------------------------------------

rule q1open_tax_table:
    output: "{x}-q1open-taxtable.tsv"
    input:
        otu_map="{x}-q1open-otus.txt",
        table="{x}-deblur-table.tsv",
        tax="{x}-q1open-shorttax.txt"
    shell:
        "scripts/q1-tax-table.R"
        " --otu_map {input.otu_map}"
        " --tax {input.tax}"
        " --table {input.table}"
        " --output {output}"

rule q1closed_tax_table:
    output: "{x}-q1closed-taxtable.tsv"
    input:
        otu_map="{x}-q1closed-otus.txt",
        table="{x}-deblur-table.tsv",
        tax="97-otu-shorttax.txt"
    shell:
        "scripts/q1-tax-table.R"
        " --otu_map {input.otu_map}"
        " --tax {input.tax}"
        " --table {input.table}"
        " --output {output}"

rule short_tax:
    output: "{x}-shorttax.txt"
    input: "{x}-tax.txt"
    shell: "scripts/clean-tax-map.R -i {input} -o {output}"

# Make Qiime 2 taxonomy tables -----------------------------------------

rule q2open_tax_table:
    output: "{x}-q2open-taxtable.tsv"
    input: table="{x}-q2open-table.tsv", tax="97-otu-shorttax.txt"
    shell: "scripts/q2-tax-table.R --tax {input.tax} --table {input.table} --output {output}"

rule q2closed_tax_table:
    output: "{x}-q2closed-taxtable.tsv"
    input: table="{x}-q2closed-table.tsv", tax="97-otu-shorttax.txt"
    shell: "scripts/q2-tax-table.R --tax {input.tax} --table {input.table} --output {output}"

# Qiime1 OTU picking --------------------------------------------------

rule q1:
    output: expand("{direction}-{pick}-{product}", direction=DIRECTIONS, pick=["q1closed", "q1open"], product=["otus.txt", "tax.txt"])
    shell: "echo 'Run Qiime1 scripts'"

# Qiime2 OTU picking --------------------------------------------------

rule closed_ref:
    output:
        table = "{x}-q1closed-table.qza",
        clusters = "{x}-q1closed-clustered-seqs.qza",
        unmatched = "{x}-q1closed-unmatched-seqs.qza"
    input:
        seqs = "{x}-seqs.qza",
        table = "{x}-table.qza",
        ref = "97_otus.qza"
    shell:
        qiime + " vsearch cluster-features-closed-reference"
        " --i-sequences {input.seqs}"
        " --i-table {input.table}"
        " --i-reference-sequences {input.ref}"
        " --p-strand both"
        " --p-perc-identity 0.97"
        " --o-clustered-table {output.table}"
        " --o-clustered-sequences {output.clusters}"
        " --o-unmatched-sequences {output.unmatched}"

rule open_ref:
    output:
        table = "{x}-q2open-table.qza",
        clusters = "{x}-q2open-clustered-seqs.qza",
        seqs = "{x}-q2open-newref-seqs.qza"
    input:
        seqs = "{x}-seqs.qza",
        table = "{x}-table.qza",
        ref = "97_otus.qza"
    shell:
        qiime + " vsearch cluster-features-open-reference"
        " --i-sequences {input.seqs}"
        " --i-table {input.table}"
        " --i-reference-sequences {input.ref}"
        " --p-strand both"
        " --p-perc-identity 0.97"
        " --o-clustered-table {output.table}"
        " --o-clustered-sequences {output.clusters}"
        " --o-new-reference-sequences {output.seqs}"

rule import_reference_otus:
    output: "97_otus.qza"
    input: "97_otus.fasta"
    shell:
        qiime + " tools import"
        " --input-path {input}"
        " --output-path {output}"
        " --type 'FeatureData[Sequence]'"

rule extract_reference_taxonomy:
    output: "97-otu-tax.txt"
    input: "gg_13_8_otus.tar.gz"
    params: archive="gg_13_8_otus/taxonomy/97_otu_taxonomy.txt"
    shell: "tar -O -f {input} -x {params.archive} > {output}"

rule udb:
    output: "{x}.udb"
    input: "{x}.fasta"
    shell: "usearch -makeudb_usearch {input} -output {output}"

rule extract_reference_otus:
    output: "97_otus.fasta"
    input: "gg_13_8_otus.tar.gz"
    params: archive="gg_13_8_otus/rep_set/97_otus.fasta"
    shell: "tar -O -f {input} -x {params.archive} > {output}"

rule download_reference:
    output: "gg_13_8_otus.tar.gz"
    params:
        url=config["taxonomy"]["url"],
        md5=config["taxonomy"]["md5"]
    shell:
        "wget {params.url} -O {output}"
        " && echo '{params.md5}  {output}' | md5sum --check"

rule deblur:
    output:
        table = "{x}-deblur-table.qza",
        seqs = "{x}-deblur-seqs.qza",
        stats = "{x}-deblur-stats.qza"
    input: "{x}-filter-seqs.qza"
    params: trim_length=lambda wc: config["deblur-trim-length"][wc.x]
    shell:
        qiime + " deblur denoise-16S"
        " --i-demultiplexed-seqs {input}"
        " --p-trim-length {params.trim_length}"
        " --p-min-reads 1"
        " --p-sample-stats"
        " --o-table {output.table}"
        " --o-representative-sequences {output.seqs}"
        " --o-stats {output.stats}"

rule filter:
    output:
        seqs="{x}-filter-seqs.qza",
        stats="{x}-filter-stats.qza"
    input: "{x}-demux-seqs.qza"
    shell:
        qiime + " quality-filter q-score"
        " --i-demux {input}"
        " --o-filtered-sequences {output.seqs}"
        " --o-filter-stats {output.stats}"

rule demux:
    output: "{x}-demux-seqs.qza"
    input: "{x}-manifest.txt"
    shell:
        qiime + " tools import"
        " --type 'SampleData[SequencesWithQuality]'"
        " --input-path {input}"
        " --output-path {output}"
        " --input-format SingleEndFastqManifestPhred33"

rule manifest:
    output: "{x}-manifest.txt"
    input: expand("{{x}}-trim/{library}_R1.fastq", library=LIBRARIES)
    params: libraries=config["libraries"]
    script: "scripts/write_manifest.py"

rule trim_merge:
    """Trim the forward and reverse primers (taking reverse complement)"""
    output: "merge-trim/{x}.fastq"
    input: "merge-notrim/{x}.fastq"
    params: forward=config["primers"]["forward"], reverse_rc=config["primers"]["reverse-rc"]
    shell: "cutadapt --front {params.forward}...{params.reverse_rc} --trimmed-only -o {output} {input}"

rule merge:
    output: "merge-notrim/{x}.fastq"
    input:
        forward="raw/{x}_R1.fastq",
        reverse="raw/{x}_R2.fastq"
    params:
        min_overlap=config["merge"]["min-overlap"],
        percent_max_diff=config["merge"]["percent-max-diff"]
    script: "scripts/merge.py"

rule trim_forward:
    output: "forward-trim/{x}.fastq"
    input: "raw/{x}_R1.fastq"
    params: primer=config["primers"]["forward"]
    shell: "cutadapt --front {params.primer} --trimmed-only -o {output} {input}"

rule trim_reverse:
    output: "reverse-trim/{x}.fastq"
    input: "raw/{x}_R2.fastq"
    params: primer=config["primers"]["reverse"]
    shell: "cutadapt --front {params.primer} --trimmed-only -o {output} {input}"
