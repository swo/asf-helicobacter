import csv, re, os.path, subprocess

def read_manifests(fns):
    manifest = []

    for fn in fns:
        with open(fn) as f:
            reader = csv.DictReader(f)
            rows = [x for x in reader]

        manifest += rows

    return manifest

def parse_manifest(manifest):
    data = {}

    for row in manifest:
        new_row = row.copy()
        # get filepaths without docker's "/data/" prefix
        new_row.update({"filepath": re.sub("^\/data\/", "", row["absolute-filepath"])})

        sample = row["sample-id"]
        direction = row["direction"]

        if sample in data:
            assert direction not in data[sample]
            data[sample][direction] = new_row
        else:
            data[sample] = {direction: new_row}

    return data

def merge(forward_fn, reverse_fn, merge_fn_template, min_overlap=4, percent_max_diff=25):
    command = ["fastq-join", forward_fn, reverse_fn, "-o", merge_fn_template, "-m", str(min_overlap), "-p", str(percent_max_diff)]
    print(" ".join(command))
    result = subprocess.run(command, capture_output=True)
    print(result.stdout)
    print(result.stderr)
    assert os.path.isfile(merge_fn_template + "join")

def merge_manifest(data):
    samples = data.keys()
    forward_fns = [data[sample]["forward"]["filepath"] for sample in samples]
    reverse_fns = [data[sample]["reverse"]["filepath"] for sample in samples]
    merge_fn_templates = ["merge/" + sample for sample in samples]
    merge_abs_fns = ["/data/" + x + "join" for x in merge_fn_templates]

    for forward_fn, reverse_fn, merge_fn_template in zip(forward_fns, reverse_fns, merge_fn_templates):
        merge(forward_fn, reverse_fn, merge_fn_template)

    merge_manifest = [{"sample-id": sample, "absolute-filepath": fn, "direction": "forward"} for sample, fn in zip(samples, merge_abs_fns)]

    return merge_manifest

def write_manifest(manifest, fn):
    with open(fn, "w") as f:
        writer = csv.DictWriter(f, fieldnames=["sample-id", "absolute-filepath", "direction"])
        writer.writeheader()

        for row in manifest:
            writer.writerow(row)

manifest = read_manifests([snakemake.input["forward"], snakemake.input["reverse"]])
data = parse_manifest(manifest)
new_manifest = merge_manifest(data)
write_manifest(new_manifest, snakemake.output[0])
