#!/usr/bin/env python3

import subprocess, tempfile, os.path

def closed_reference(prefix, force=False):
    dest = "{}-q1closedref-otus.txt".format(prefix)

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

def open_reference(prefix, force=False):
    dest1 = "{}-q1openref-otus.txt".format(prefix)
    dest2 = "{}-q1openref-tax.txt".format(prefix)

    if os.path.exists(dest1) and os.path.exists(dest2) and not force:
        return None

    with tempfile.TemporaryDirectory() as tmp_dir:
        command = ["pick_open_reference_otus.py", "-f", "-p", "params.config", "-i", "../{}-deblur-seqs.fasta".format(prefix), "-o", tmp_dir]

        subprocess.call(command)

        src1 = os.path.join(tmp_dir, "final_otu_map_mc2.txt")
        subprocess.call(["mv", src1, dest1])

        src2 = os.path.join(tmp_dir, "uclust_assigned_taxonomy", "rep_set_tax_assignments.txt")
        subprocess.call(["mv", src2, dest2])

    return (dest1, dest2)

prefixes = ["{}-{}".format(direction, trim) for direction in ["forward", "reverse", "merge"] for trim in ["trim", "notrim"]]

for prefix in prefixes:
    print("Closed reference")
    result = closed_reference(prefix)
    print(prefix, result)

    print("Open reference")
    result = open_reference(prefix)
    print(prefix, result)
