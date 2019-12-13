import csv, re, os.path

group = re.match("(.*)-manifest.txt", snakemake.output[0]).group(1)

with open(snakemake.output[0], "w") as f:
    writer = csv.DictWriter(f, fieldnames=["sample-id", "absolute-filepath", "direction"])
    writer.writeheader()

    for filepath in snakemake.input:
        match = re.match("(.*)/(.*).fastq", filepath)
        assert group == match.group(1)

        library = match.group(2)
        sample_name = snakemake.params["libraries"][library]
        sample_id = f"{group}-{sample_name}"

        abs_filepath = "/data/" + filepath
        row = {"sample-id": sample_id, "absolute-filepath": abs_filepath, "direction": "forward"}

        writer.writerow(row)
