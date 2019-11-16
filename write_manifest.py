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

    sample_id = f"{directory}-{direction}-{library}"

    absolute_filepath = "/data/" + fn

    return {"sample-id": sample_id, "absolute-filepath": absolute_filepath, "direction": direction}
    

rows = [classify_file(x) for x in snakemake.input]

with open(snakemake.output[0], "w") as f:
    w = csv.DictWriter(f, fieldnames=["sample-id", "absolute-filepath", "direction"])
    w.writeheader()

    for row in rows:
        w.writerow(row)
