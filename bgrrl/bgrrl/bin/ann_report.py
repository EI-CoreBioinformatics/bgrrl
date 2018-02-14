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
        # sample = os.path.basename(sdir)
        with open(os.path.join(sdir, sample + ".ratt_report.tsv"), "w") as ratt_out:
            header = list()
            for d in sorted(next(os.walk(sdir))[1]):
                # print("D=", d)
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
                # keys = list(sorted(refcounts.get(ref, Counter()).keys()))
                # print("KEYS=", keys)
                keys = FEATURES
                if not header:
                    header = ["Reference"] + keys + keys + keys + ["Total[ref]", "Total[transferred]", "Total[%]"]
                    print(*header, sep="\t", file=ratt_out)
                     
                row.extend([refcounts.get(ref, Counter())[key] for key in keys])
                row.extend([trfcounts[key] for key in keys])
                rtotal, ttotal = 0, 0
                for key in keys:
                    rcount = refcounts.get(ref, Counter())[key]
                    rtotal += rcount
                    ttotal += trfcounts[key]
                    row.append(trfcounts[key]/rcount if rcount > 0 else "NA")
                row.extend([rtotal, ttotal, ttotal/rtotal if rtotal > 0 else "NA"])

                # row.extend([trfcounts[key]/refcounts.get(d, Counter())[key] for key in keys])
                # row.extend([sum(refcounts.get(d, Counter()).values()), sum(trfcounts.values()), sum(trfcounts.values())/sum(refcounts.get(d, Counter()).values())])
                
                print(*row, sep="\t", file=ratt_out)
            



            




if __name__ == "__main__":
    main()
