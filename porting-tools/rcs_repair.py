#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
rcs_repair.py  –  bulk-friendly RCS file rebuilder (Python 2.7)
----------------------------------------------------------------
• Extracts every trunk revision still readable from a corrupt CVS/RCS *,v* file.
• Builds a clean *,v* under `./repaired-rcs-files/<relative‑path>` so you can
  later `cp -a` them back into the real CVS repository.
• Skips branch revisions (>2 numeric components).
• Uses RCS CLI tools only (`co`, `ci`, `rlog`).
"""
from __future__ import print_function, unicode_literals
import os, re, sys, shutil, subprocess

INPUT_ROOT = "/users/livesey/cvsroot-mlsdev-backup/cvsroot/mlspgs"
OUTROOT = "./repaired-rcs-files"
WORK_TMP = "__tmp_workfile__"  # reused per‑revision

RLOG_RE = re.compile(
    r"^revision\s+(?P<rev>\S+)\s*\n"
    r"date:\s+(?P<date>\d{4}/\d{2}/\d{2})\s+"
    r"(?P<time>\d{2}:\d{2}:\d{2});\s+"
    r"author:\s+(?P<author>\w+);.*?state:\s+(?P<state>\w+);",
    re.M | re.S,
)


def is_branch_rev(rev):
    return rev.count(".") > 2


def run(cmd, **kw):
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, **kw)
    out, err = p.communicate()
    return p.returncode, out, err


def co_text(rcs_path, rev):
    rc, out, _ = run(["co", "-q", "-x,v", "-p" + rev, rcs_path])
    return out if rc == 0 else None


def ci_init(workfile, rev, author, datestr, logmsg, state):
    """Create brand-new RCS file (use -i)."""
    ci_cmd = [
        "ci",
        "-q",
        "-i",
        "-x,v",
        "-t-",
        "-w" + author,
        "-d" + datestr,
        "-m" + logmsg,
        "-s" + state,
        "-r" + rev,
        workfile,
    ]
    _run_ci(ci_cmd)


def ci_append(workfile, rev, author, datestr, logmsg, state):
    """Append further revisions (-f to skip locks)."""
    ci_cmd = [
        "ci",
        "-q",
        "-f",
        "-x,v",
        "-t-",
        "-w" + author,
        "-d" + datestr,
        "-m" + logmsg,
        "-s" + state,
        "-r" + rev,
        workfile,
    ]
    _run_ci(ci_cmd)


def _run_ci(cmd):
    p = subprocess.Popen(cmd, stdin=subprocess.PIPE)
    p.communicate(b"imported by repair script\n.\n")
    if p.returncode:
        raise subprocess.CalledProcessError(p.returncode, cmd)


def repair(rcs_path):
    # ── target paths ───────────────────────────────────────────────
    rel_path = os.path.relpath(rcs_path, INPUT_ROOT)
    repaired = os.path.join(OUTROOT, rel_path)
    workfile = os.path.join(
        os.path.dirname(repaired), os.path.basename(repaired)[:-2]
    )  # strip ',v'

    sys.stderr.write("Rebuilding %s -> %s\n" % (rcs_path, repaired))

    # fresh output dir
    outdir = os.path.dirname(repaired)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    if os.path.exists(repaired):
        os.remove(repaired)

    # ── collect trunk revisions (oldest → newest) ─────────────────
    rc, rlog_out, _ = run(["rlog", rcs_path])
    if rc:
        sys.exit("rlog failed on %s" % rcs_path)
    revs = RLOG_RE.findall(rlog_out.decode("utf-8", "replace"))
    revs.reverse()

    good = bad = 0
    first = True

    for rev, date, hms, author, state in revs:
        if is_branch_rev(rev):
            sys.stderr.write("   branch rev %s skipped\n" % rev)
            continue

        text = co_text(rcs_path, rev)
        if text is None:
            bad += 1
            sys.stderr.write(" n  %s\n" % rev)
            continue

        # write checkout to staging file
        open(WORK_TMP, "wb").write(text)
        datestr = date.replace("/", "-") + " " + hms
        logmsg = "rev " + rev

        if first:
            # clean any stale RCS left by earlier aborted runs
            for stale in (
                workfile + ",v",
                os.path.join(
                    os.path.dirname(workfile), "RCS", os.path.basename(workfile) + ",v"
                ),
            ):
                try:
                    os.remove(stale)
                except OSError:
                    pass

        if not first:
            # lock the head revision for the current user before appending
            run(["rcs", "-q", "-l", repaired])

        shutil.copy2(WORK_TMP, workfile)

        if first:
            ci_init(workfile, rev, author, datestr, logmsg, state)  # uses  ci -i
            first = False
        else:
            ci_append(workfile, rev, author, datestr, logmsg, state)  # uses ci -f

        try:
            os.remove(workfile)  # remove plain workfile
        except OSError:
            pass

        good += 1
        sys.stderr.write(" y  %s\n" % rev)

    # tidy staging files
    for f in (WORK_TMP, WORK_TMP + ",v"):
        try:
            os.remove(f)
        except OSError:
            pass

    sys.stderr.write("Done. %d good, %d bad revisions.\n" % (good, bad))
    sys.stderr.write("Repaired file written to: %s\n" % repaired)


if __name__ == "__main__":
    if len(sys.argv) != 2 or not sys.argv[1].endswith(",v"):
        sys.exit("Usage: %s path/to/bad/File,v" % sys.argv[0])
    repair(os.path.abspath(sys.argv[1]))
