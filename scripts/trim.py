import subprocess, os.path

def trim(input_fn, output_fn, primer):
    command = ["cutadapt", "--front", primer, "-o", output_fn, "--trimmed-only", input_fn]
    print(" ".join(command))
    results = subprocess.run(command, capture_output=True)
    print(results.stdout.decode("utf-8"))
    print(results.stderr.decode("utf-8"))

    assert os.path.isfile(output_fn)

trim(snakemake.input[0], snakemake.output[0], snakemake.params["primer"])
