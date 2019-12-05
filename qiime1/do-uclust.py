#!/usr/bin/env python3

import subprocess, tempfile, os.path

def run_uclust(prefix, force=False):
    dest = "{}-otus.txt".format(prefix)

    if os.path.exists(dest) and not force:
        return None

    with tempfile.TemporaryDirectory() as tmp_dir:
        # -f overwrites directories
        # -p points to the params file, which has an option for reverse strand
        # searching (for the reverse sequences)
        command = ["pick_closed_reference_otus.py", "-f", "-p", "params.config", "-i", "../{}-deblur-seqs.fasta".format(prefix), "-o", tmp_dir]

        # An alternative is to use pick_otus.py with -z (reverse strand search), -C (suppress new clusters, i.e., do reference-based search),
        # and -m uclust_ref (i.e., use the default closed reference method):
        # command = ["pick_otus.py", "-z", "-C", "-m", "uclust_ref", "-i", "../{}-deblur-seqs.fasta".format(prefix), "-o", tmp_dir]

        subprocess.call(command)

        src = os.path.join(tmp_dir, "uclust_ref_picked_otus", "{}-deblur-seqs_otus.txt".format(prefix))
        command = ["mv", src, dest]
        subprocess.call(command)

    return dest

prefixes = ["{}-{}".format(direction, trim) for direction in ["forward", "reverse", "merge"] for trim in ["trim", "notrim"]]

for prefix in prefixes:
    result = run_uclust(prefix)
    print(prefix, result)
