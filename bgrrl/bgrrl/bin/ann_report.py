import os
import sys
import csv
import glob
import argparse

from collections import Counter

FEATURES = ["gene", "CDS", "tRNA", "rRNA", "tmRNA", "ncRNA"]

def countGFFFeatures(_in):
    return Counter(row[2] for row in csv.reader(_in, delimiter="\t") if row and not row[0].startswith("#"))


def main(args_in=sys.argv[1:]):
    ap = argparse.ArgumentParser()
    ap.add_argument("indir", type=str, default=".")
    ap.add_argument("--report-dir", type=str, default="")
    ap.add_argument("--ref-dir", type=str, default=".")

    args = ap.parse_args(args_in)
    print(args)


    gffs = list()
    owalk = os.walk(args.ref_dir)

    refcounts = dict()

    for d, dirs, files in os.walk(args.ref_dir):
        if os.path.basename(d) == "gff":
            for f in files:
                if f.endswith(".gff"):
                    with open(os.path.join(d, f)) as fin:
                        refcounts[f.strip(".gff")] = countGFFFeatures(fin)

    print(refcounts.keys())

    samples = next(os.walk(args.indir))[1]
    for sample in samples:
        sdir = os.path.join(args.indir, sample)
        with open(os.path.join(sdir, sample + ".ratt_report.tsv"), "w") as ratt_out:
            header, rows = list(), list()
            for d in sorted(next(os.walk(sdir))[1]):
                try:
                    f = glob.glob(os.path.join(sdir, d, "*.final.gff"))[0]
                except:
                    print("No .final.gff in {}.".format(os.path.join(sdir, d)), file=sys.stderr)
                    continue
                with open(f) as fin:
                    trfcounts = countGFFFeatures(fin)

                ref = d.replace(sample + ".", "")
                print("SAMPLE=", sample, "REF=", ref)
                row = [ref]
                if not header:
                    header = ["Reference"]
                    header.extend(map(lambda x:x+"[ref]", FEATURES))
                    header.extend(map(lambda x:x+"[transferred]", FEATURES))
                    header.extend(map(lambda x:x+"[%]", FEATURES))
                    header.extend(["Total[ref]", "Total[transferred]", "Total[%]"])
                    print(*header, sep="\t", file=ratt_out)
                     
                row.extend([refcounts.get(ref, Counter())[feature] for feature in FEATURES])
                row.extend([trfcounts[feature] for feature in FEATURES])
                rtotal, ttotal = 0, 0
                for feature in FEATURES:
                    rcount = refcounts.get(ref, Counter())[feature]
                    rtotal += rcount
                    ttotal += trfcounts[feature]
                    row.append(trfcounts[feature]/rcount if rcount > 0 else "NA")
                row.extend([rtotal, ttotal, ttotal/rtotal if rtotal > 0 else "NA"])
                rows.append(row)

            for row in sorted(rows, key=lambda x:x[-1] if x[-1] != "NA" else -1, reverse=True):            
                print(*row, sep="\t", file=ratt_out)

    return True  




if __name__ == "__main__":
    main()