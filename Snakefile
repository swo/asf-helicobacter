import uuid

configfile: "config.yaml"

DIRECTIONS = ["forward", "reverse", "merge"]
TRIMS = ["notrim", "trim"]
BLURS = ["derep", "deblur"]
GROUPS = expand("{direction}-{trim}", direction=DIRECTIONS, trim=TRIMS)
OTU_PICKINGS = ["closedref", "openref", "rdp"]

# library names from raw, input files
LIBRARIES = config["libraries"].keys()

ANALYSES = expand("{otu_picking}-{product}", otu_picking=OTU_PICKINGS, product=["taxtable.tsv", "taxaplot.pdf"]) + \
    ["openref-usearch.tsv"]

# note that current directory is mounted as /data in the container
qiime = "docker run -t -i -v $(pwd):/data qiime2/core:2019.10 qiime"

rule all:
    input:
        expand("{group}-{blur}-{analysis}", group=GROUPS, blur=BLURS, analysis=ANALYSES) + \
        # ["all-taxaplot.pdf"]
        []

rule clean:
    shell:
        "rm -rf"
        " *.txt *.qza *.biom *.tsv *.fasta *.b6 *.log *.udb *.tar.gz *.pdf" + \
        " gg_13_8_otus/"
        " " + " ".join(expand("{group}/*", group=GROUPS))

rule taxa_plot_all:
    output: plot="all-taxaplot.pdf", table="all-taxaplot-table.tsv"
    input: expand("{direction}-{trim}-{picking}-taxtable.tsv", direction=DIRECTIONS, trim=TRIMS, picking=OTU_PICKINGS)
    script: "scripts/taxa-plot-all.R"

rule taxa_plot:
    output: "{x}-taxaplot.pdf"
    input: "{x}-taxtable.tsv"
    script: "scripts/taxa-plot.R"

rule usearch_openref:
    output: "{x}-openref-usearch.tsv"
    input:
        table="{x}-openref-table.tsv",
        seqs="{x}-openref-newref-seqs.fasta",
        otus="97_otus.udb",
        taxonomy="97_otu_taxonomy.txt"
    script: "scripts/openref-usearch.R"

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
    """Shadow to prevent race condition with intermediate file"""
    output: "{x}-table.tsv"
    input: "{x}-table.qza"
    params: tempdir=lambda wc: wc["x"] + "-" + uuid.uuid4().hex
    shell:
        qiime + " tools export"
        " --input-path {input}"
        " --output-path tmp/{params.tempdir}"
        " && biom convert -i tmp/{params.tempdir}/feature-table.biom -o {output} --to-tsv"
        " && rm -rf tmp/{params.tempdir}"

rule export_taxonomy:
    """Taxonomy qza to tsv"""
    output: "{x}-tax.tsv"
    input: "{x}-tax.qza"
    params: tempdir=lambda wc: wc["x"] + "-" + uuid.uuid4().hex
    shell:
        qiime + " tools export"
        " --input-path {input}"
        " --output-path {params.tempdir}"
        " && mv {params.tempdir}/taxonomy.tsv {output}"
        " && rmdir {params.tempdir}"

rule rdp_tax_table:
    output: "{x}-rdp-taxtable.tsv"
    input: table="{x}-table.tsv", tax="{x}-rdp-tax.tsv"
    script: "scripts/tax-table.R"

rule ref_tax_table:
    output: "{x}ref-taxtable.tsv"
    input: table="{x}ref-table.tsv", tax="97_otu_taxonomy.txt"
    script: "scripts/tax-table.R"

rule rdp:
    output: "{x}-rdp-tax.qza"
    input: reads="{x}-seqs.qza", classifier="rdp-classifier.qza"
    shell:
        qiime + " feature-classifier classify-sklearn"
        " --i-classifier {input.classifier}"
        " --i-reads {input.reads}"
        " --verbose"
        " --p-pre-dispatch 1"
        " --o-classification {output}"

rule download_classifier:
    output: "rdp-classifier.qza"
    params:
        url=config["classifier"]["url"],
        md5=config["classifier"]["md5"]
    shell:
        "wget {params.url} -O {output}"
        " && echo '{params.md5}  {params.md5}' | md5sum "

rule closed_ref:
    output:
        table = "{x}-closedref-table.qza",
        clusters = "{x}-closedref-clustered-seqs.qza",
        unmatched = "{x}-closedref-unmatched-seqs.qza"
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
        table = "{x}-openref-table.qza",
        clusters = "{x}-openref-clustered-seqs.qza",
        seqs = "{x}-openref-newref-seqs.qza"
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
    output: "97_otu_taxonomy.txt"
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
        " && echo '{params.md5}  {params.md5}' | md5sum "

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

rule derep:
    output:
        table = "{x}-derep-table.qza",
        seqs = "{x}-derep-seqs.qza"
    input: "{x}-filter-seqs.qza"
    shell:
        qiime + " vsearch dereplicate-sequences"
        " --i-sequences {input}"
        " --o-dereplicated-table {output.table}"
        " --o-dereplicated-sequences {output.seqs}"

rule filter:
    output: seqs="{x}-filter-seqs.qza", stats="{x}-filter-stats.qza"
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
    input: expand("{{x}}/{library}_R1.fastq", library=LIBRARIES)
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
    input: forward="forward-notrim/{x}.fastq", reverse="reverse-notrim/{x}.fastq"
    params:
        min_overlap=config["merge"]["min-overlap"],
        percent_max_diff=config["merge"]["percent-max-diff"]
    script: "scripts/merge.py"

rule trim_forward:
    output: "forward-trim/{x}.fastq"
    input: "forward-notrim/{x}.fastq"
    params: primer=config["primers"]["forward"]
    shell: "cutadapt --front {params.primer} --trimmed-only -o {output} {input}"

rule trim_reverse:
    output: "reverse-trim/{x}.fastq"
    input: "reverse-notrim/{x}.fastq"
    params: primer=config["primers"]["reverse"]
    shell: "cutadapt --front {params.primer} --trimmed-only -o {output} {input}"

rule copy_reverse:
    output: "reverse-notrim/{x}_R1.fastq"
    input: "raw/{x}_R2.fastq"
    shell: "cp {input} {output}"

rule copy_forward:
    output: "forward-notrim/{x}_R1.fastq"
    input: "raw/{x}_R1.fastq"
    shell: "cp {input} {output}"
