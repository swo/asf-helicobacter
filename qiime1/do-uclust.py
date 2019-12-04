#!/usr/bin/env python3

import subprocess, tempfile, os.path

def run_uclust(prefix, force=False):
    dest = "{}-otus.txt".format(prefix)

    if os.path.exists(dest) and not force:
        return None

    with tempfile.TemporaryDirectory() as tmp_dir:
        # -z enables reverse strand
        # -C suppresses new clusters (i.e., does reference-based)
        # -m selects the method
        command = ["pick_otus.py", "-z", "-C", "-m", "uclust_ref", "-i", "../{}-deblur-seqs.fasta".format(prefix), "-o", tmp_dir]
        subprocess.call(command)

        src = os.path.join(tmp_dir, "{}-deblur-seqs_otus.txt".format(prefix))
        command = ["mv", src, dest]
        subprocess.call(command)

    return dest

prefixes = ["{}-{}".format(direction, trim) for direction in ["forward", "reverse", "merge"] for trim in ["trim", "notrim"]]

for prefix in prefixes:
    result = run_uclust(prefix)
    print(prefix, result)
