#!/usr/bin/env python3
"""With more than one thread, which transcriptome alignment of a read is marked primary depends on the thread that maps
the read (ReadAlign_quantTranscriptome.cpp draws it from a per-thread RNG), so it varies from run to run even upstream.
bam_body_hash.py --sorted --anyPrimary must hash such BAMs alike, and still tell apart BAMs that differ otherwise."""
import gzip, os, struct, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))


def record(name, pos, flag):
    qname = name.encode() + b"\0"
    core = struct.pack("<iiBBHHHiiii", 0, pos, len(qname), 255, 0, 0, flag, 0, -1, -1, 0) + qname
    return struct.pack("<i", len(core)) + core


def bam(path, recs):
    ref = b"tr1\0"
    head = b"BAM\1" + struct.pack("<i", 0) + struct.pack("<i", 1) + struct.pack("<i", len(ref)) + ref + struct.pack("<i", 1000)
    with gzip.open(path, "wb") as f:
        f.write(head + b"".join(recs))


def digest(path, *opts):
    run = subprocess.run([sys.executable, os.path.join(HERE, "bam_body_hash.py"), *opts, path], capture_output=True, text=True)
    if run.returncode != 0:  # a failed run matches nothing, not even another failed run
        print("    bam_body_hash.py " + " ".join(opts) + " failed: " + run.stderr.strip().splitlines()[-1])
        return object()
    return run.stdout.strip()


SEC = 0x100
cases = {
    "ref":       [record("r1", 10, 0), record("r1", 50, SEC), record("r2", 7, 0)],
    "otherPrim": [record("r2", 7, 0), record("r1", 50, 0), record("r1", 10, SEC)],
    "noPrim":    [record("r1", 10, SEC), record("r1", 50, SEC), record("r2", 7, 0)],
    "twoPrim":   [record("r1", 10, 0), record("r1", 50, 0), record("r2", 7, 0)],
    "otherPos":  [record("r1", 10, 0), record("r1", 51, SEC), record("r2", 7, 0)],
}
fails = 0
with tempfile.TemporaryDirectory() as tmp:
    paths = {}
    for k, recs in cases.items():
        paths[k] = os.path.join(tmp, k + ".bam")
        bam(paths[k], recs)

    def expect(what, a, b, equal, *opts):
        global fails
        same = digest(paths[a], *opts) == digest(paths[b], *opts)
        ok = same == equal
        fails += not ok
        print(("PASS " if ok else "FAIL ") + what)

    expect("another primary alignment of the same read hashes alike", "ref", "otherPrim", True, "--sorted", "--anyPrimary")
    expect("without --anyPrimary the primary choice still counts", "ref", "otherPrim", False, "--sorted")
    expect("a read with no primary alignment is told apart", "ref", "noPrim", False, "--sorted", "--anyPrimary")
    expect("a read with two primary alignments is told apart", "ref", "twoPrim", False, "--sorted", "--anyPrimary")
    expect("a different alignment is told apart", "ref", "otherPos", False, "--sorted", "--anyPrimary")
sys.exit(1 if fails else 0)
