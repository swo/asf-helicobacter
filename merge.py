import csv, re, os.path, subprocess, os, tempfile

def merge(forward_fn, reverse_fn, merge_fn, min_overlap=4, percent_max_diff=25):
    with tempfile.TemporaryDirectory() as tmp_dir:
        command = ["fastq-join", forward_fn, reverse_fn, "-o", tmp_dir + "/", "-m", str(min_overlap), "-p", str(percent_max_diff)]

        print(" ".join(command))
        result = subprocess.run(command, capture_output=True)
        print(result.stdout)
        print(result.stderr)

        from_path = tmp_dir + "/join"
        os.replace(from_path, merge_fn)

    assert os.path.isfile(merge_fn)

merge(snakemake.input["forward"], snakemake.input["reverse"], snakemake.output[0])
