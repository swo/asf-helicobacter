import csv, re, os.path

def classify_file(fn):
    directory, basename = os.path.split(fn)

    m = re.search("^(.*)_R([12])", basename)
    library = m.group(1)
    r_direction = m.group(2)

    if r_direction == "1":
        direction = "forward"
    elif r_direction == "2":
        direction = "reverse"
    else:
        raise RuntimeError(f"Bad R-number: {r_direction}")

    manifest = f"{directory}-{direction}-manifest.txt"
    sample_id = f"{directory}-{direction}-{library}"

    absolute_filepath = "/data/" + fn

    return {"manifest": manifest, "row": {"sample-id": sample_id, "absolute-filepath": absolute_filepath, "direction": direction}}
    
data = [classify_file(x) for x in snakemake.input]
known_manifests = set([x["manifest"] for x in data])

assert known_manifests == set(snakemake.output)

for output_fn in snakemake.output:
    rows = [x["row"] for x in data if x["manifest"] == output_fn]

    with open(output_fn, "w") as f:
        w = csv.DictWriter(f, fieldnames=["sample-id", "absolute-filepath", "direction"])
        w.writeheader()

        for row in rows:
            w.writerow(row)
