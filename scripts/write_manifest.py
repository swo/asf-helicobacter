import csv, re, os.path

def classify_file(fn):
    directory, basename = os.path.split(fn)

    m = re.match("^(.*)_R1.fastq", basename)

    if m is None:
        raise RuntimeError(f"bad filename {basename}")

    library = m.group(1)

    manifest = f"{directory}-manifest.txt"
    sample_id = f"{directory}-{library}"

    absolute_filepath = "/data/" + fn

    return {"manifest": manifest, "library": library, "row": {"sample-id": sample_id, "absolute-filepath": absolute_filepath, "direction": "forward"}}

data = [classify_file(x) for x in snakemake.input]

known_manifests = set([x["manifest"] for x in data])
assert known_manifests == set(snakemake.output)

known_libraries = set([x["library"] for x in data])
assert known_libraries == set(snakemake.params["libs"])

for output_fn in snakemake.output:
    rows = [x["row"] for x in data if x["manifest"] == output_fn]

    with open(output_fn, "w") as f:
        w = csv.DictWriter(f, fieldnames=["sample-id", "absolute-filepath", "direction"])
        w.writeheader()

        for row in rows:
            w.writerow(row)
