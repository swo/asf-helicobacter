import subprocess, os.path

def trim(input_fn, output_fn, primer):
    command = ["cutadapt", "--front", primer, "-o", output_fn, input_fn]
    print(" ".join(command))
    results = subprocess.run(command, capture_output=True)
    print(results.stdout)
    print(results.stderr)

    assert os.path.isfile(output_fn)

trim(snakemake.input[0], snakemake.output[0], snakemake.params["primer"])
