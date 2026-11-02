#!/usr/bin/env python3
import pandas as pd
import argparse

p = argparse.ArgumentParser()
p.add_argument("--input", required=True)      # sample_bracken.kreport
p.add_argument("--threshold", type=float, default=1.0)  # percent, not fraction
p.add_argument("--output", required=True)
args = p.parse_args()

cols = ["percent", "reads_clade", "reads_direct", "rank_code", "ncbi_taxid", "sci_name"]
df = pd.read_csv(args.input, sep="\t", names=cols, header=None)

keep = df[df.percent >= args.threshold]
dropped = df[df.percent < args.threshold]

if len(dropped):
    other_row = pd.DataFrame([{
        "percent": dropped.percent.sum(),
        "reads_clade": dropped.reads_clade.sum(),
        "reads_direct": dropped.reads_direct.sum(),
        "rank_code": "S",
        "ncbi_taxid": 0,
        "sci_name": f"  Other (<{args.threshold}%, n={len(dropped)} taxa)"
    }])
    df = pd.concat([keep, other_row], ignore_index=True)
else:
    df = keep

df.to_csv(args.output, sep="\t", index=False, header=False)