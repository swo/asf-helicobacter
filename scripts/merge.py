import csv, re, os.path, subprocess, os, tempfile

def merge(forward_fn, reverse_fn, merge_fn, min_overlap=None, percent_max_diff=None):
    with tempfile.TemporaryDirectory() as tmp_dir:
        command = ["fastq-join", forward_fn, reverse_fn, "-o", tmp_dir + "/"]

        if min_overlap is not None:
            command += ["-m", str(min_overlap)]

        if percent_max_diff is not None:
            command += ["-p", str(percent_max_diff)]

        print(" ".join(command))
        result = subprocess.run(command, capture_output=True)
        print(result.stdout.decode("utf-8"))
        print(result.stderr.decode("utf-8"))

        from_path = tmp_dir + "/join"
        os.replace(from_path, merge_fn)

    assert os.path.isfile(merge_fn)

merge(snakemake.input["forward"], snakemake.input["reverse"], snakemake.output[0], min_overlap=snakemake.params["min_overlap"], percent_max_diff=snakemake.params["percent_max_diff"])
