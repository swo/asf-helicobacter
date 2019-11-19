import re, itertools, csv

def clean_taxon(x):
    fields = x.split("; ")
    assert len(fields) == 7
    matches = [re.match("([kpcofgs])__(.*)", x) for x in fields]
    levels = [x.groups(1) for x in matches]
    values = [x.groups(2) for x in matches]
    assert levels == ["k", "p", "c", "o", "f", "g", "s"]

    new_values = list(itertools.takewhile(lambda x: x != "", values))
    return "; ".join([f"{level}__{value}" for level, value in zip(levels, new_values)])

with open(snakemake.input[0], "r") as fin, open(snakemake.output[0], "w") as fout:
    reader = csv.reader(fin, dialect="excel-tab")

    writer = csv.DictWriter(fout, dialect="excel-tab", fieldnames=["Feature ID", "Taxon", "Confidence"])
    writer.writeheader()

    for row in reader:
        fid = row[0]
        raw_taxon = row[1]
        taxon = clean_taxon(raw_taxon)
        writer.writerow({"Feature ID": fid, "Taxon": taxon, "Confidence": "1.0"})
