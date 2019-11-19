import requests, hashlib

def read_md5(fn):
    with open(fn, "rb") as f:
        contents = f.read()

    return hashlib.md5(contents).hexdigest()

r = requests.get(snakemake.params["url"])

with open(snakemake.output[0], "wb") as f:
    f.write(r.content)

assert read_md5(snakemake.output[0]) == snakemake.params["md5"]
