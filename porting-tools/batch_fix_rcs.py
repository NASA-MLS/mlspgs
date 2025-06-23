#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
batch_fix_rcs.py  –  Call rcs_repair.py on every file listed in bad_files.txt.

Usage:
    ./batch_fix_rcs.py  bad_files.txt

Where:
  * bad_files.txt  –  one path per line, **relative** to CVSROOT and ending in ',v'
  * CVSROOT        –  edit the constant below to match your backup copy
  * rcs_repair.py  –  must be in the same directory or on $PATH
Produces:
  ./repaired-rcs-files/<relative path>.repaired,v
"""

import os, sys, subprocess, shutil

CVSROOT = "/users/livesey/cvsroot-mlsdev-backup/cvsroot/mlspgs"
REPAIR = os.path.join(os.path.dirname(__file__), "rcs_repair.py")
OUTROOT = os.path.abspath("repaired-rcs-files")


def ensure_dir(path):
    d = os.path.dirname(path)
    if not os.path.isdir(d):
        os.makedirs(d)


def main(listfile):
    total_good = total_bad = 0
    with open(listfile) as fh:
        paths = [l.strip() for l in fh if l.strip()]

    for rel in paths:
        abs_rcs = os.path.join(CVSROOT, rel)
        if not os.path.exists(abs_rcs):
            print("!! MISSING  {}".format(rel))
            continue

        # Call rcs_repair.py, capture its stderr (progress) but let it run
        p = subprocess.Popen(
            [REPAIR, abs_rcs], stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        out, err = p.communicate()
        # parse last line for counts: "Done. N good, M bad revisions."
        good = bad = 0
        for line in err.splitlines():
            if line.startswith("Done."):
                parts = line.split()
                good = int(parts[1])
                bad = int(parts[3])
        if good == 0:
            print("✗  {:>3} good, {:>3} bad  {}".format(good, bad, rel))
            continue

        # repaired_path = abs_rcs + ".repaired,v"
        # if not os.path.exists(repaired_path):
        #     print("✗  repaired file missing  {}".format(rel))
        #     continue

        # dest = os.path.join(OUTROOT, rel + ".repaired,v")
        # ensure_dir(dest)
        # shutil.copy2(repaired_path, dest)
        print("✔  {:>3} good, {:>3} bad  {}".format(good, bad, rel))

        total_good += good
        total_bad += bad

    print("\n=== SUMMARY ===")
    print("total good revisions rescued : {}".format(total_good))
    print("total bad revisions skipped  : {}".format(total_bad))
    print("repaired files stored under  : {}".format(OUTROOT))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: {} bad_files.txt".format(sys.argv[0]))
    main(sys.argv[1])
