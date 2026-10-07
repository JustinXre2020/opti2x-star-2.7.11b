#!/usr/bin/env python3
"""Hash a BAM's alignment records only (skip text header, which embeds the command line).
With --sorted, the records are hashed in sorted order, for BAMs whose record order varies from run to run.
With --anyPrimary (needs --sorted), which alignment of a read is primary does not count, only how many records of each
read are primary: STAR picks the primary transcriptome alignment with a per-thread RNG, so with more than one thread
it varies from run to run."""
import gzip, hashlib, struct, sys
from collections import Counter
args = sys.argv[1:]
ordered = "--sorted" in args
anyPrimary = "--anyPrimary" in args
assert ordered or not anyPrimary, "--anyPrimary needs --sorted"
path = [a for a in args if a not in ("--sorted", "--anyPrimary")][0]
d = gzip.open(path).read()
assert d[:4] == b"BAM\1"
(l_text,) = struct.unpack_from("<i", d, 4)
p = 8 + l_text
(n_ref,) = struct.unpack_from("<i", d, p); p += 4
for _ in range(n_ref):
    (l_name,) = struct.unpack_from("<i", d, p); p += 4 + l_name + 4
body = d[p:]
SECONDARY = 0x100
if ordered:
    recs, q = [], 0
    primaries = Counter()
    while q < len(body):
        (n,) = struct.unpack_from("<i", body, q)
        rec = body[q:q + 4 + n]; q += 4 + n
        if anyPrimary:
            (flag,) = struct.unpack_from("<H", rec, 18)
            name = rec[36:36 + rec[12]]
            primaries[name] += 0 if flag & SECONDARY else 1
            rec = rec[:18] + struct.pack("<H", flag & ~SECONDARY) + rec[20:]
        recs.append(rec)
    body = b"".join(sorted(recs))
    if anyPrimary:
        body += b"".join(name + struct.pack("<i", k) for name, k in sorted(primaries.items()))
print(hashlib.sha256(body).hexdigest()[:16])
