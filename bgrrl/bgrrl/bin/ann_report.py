import os
import sys
import csv
import glob
import argparse

from collections import Counter, namedtuple, OrderedDict

# FEATURES = ["gene", "CDS", "tRNA", "rRNA", "tmRNA", "ncRNA"]
FEATURES = ["gene", "tRNA", "rRNA", "tmRNA", "ncRNA"]
ANN_REPORT_HEADER = ["Reference"]
ANN_REPORT_HEADER.extend(["Transfer coordinate correction?"])
ANN_REPORT_HEADER.extend(["#Features[ref]", "#Features[transferred]", "Transferred[%]"])
ANN_REPORT_HEADER.extend(map(lambda x:x+"[ref]", FEATURES))
ANN_REPORT_HEADER.extend(map(lambda x:x+"[transferred]", FEATURES))
ANN_REPORT_HEADER.extend(map(lambda x:x+"[%]", FEATURES))
ANN_REPORT_HEADER.extend(["Total[ref]", "Total[transferred]", "Total[%]"])

GFFeature = namedtuple("GFFeature", "seqname source feature start end score strand frame attribute".split(" "))




def countGFFFeatures(_in):
    fcounter = dict((feature, set()) for feature in FEATURES)
    genecounter = {None}
    # genecounter = {None: set()}

    # return Counter(row[2] for row in csv.reader(_in, delimiter="\t") if row and not row[0].startswith("#"))
    for row in csv.reader(_in, delimiter="\t"):
        if row and not row[0].startswith("#"):
            attr = OrderedDict((item.split("=")[0], item.split("=")[1]) for item in row[8].split(";"))
            ltag = attr.get("locus_tag", None)
            if ltag is not None and row[2] != "CDS":
                fcounter.setdefault(row[2], set()).add(ltag)
                # genecounter.add(attr.get("ID", None))
                genecounter.add(ltag)
    
    for k in fcounter:
        if k != "gene":
            fcounter["gene"].difference_update(fcounter[k])
    genecounter.remove(None)

    return dict((k, len(fcounter[k])) for k in fcounter), genecounter
            
        
 


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
    refgenecounts = dict()

    for d, dirs, files in os.walk(args.ref_dir):
        if os.path.basename(d) == "gff":
            for f in files:
                if f.endswith(".gff"):
                    with open(os.path.join(d, f)) as fin:
                        refcounts[f.strip(".gff")], refgenecounts[f.strip(".gff")] = countGFFFeatures(fin)

    print(refcounts.keys())

    samples = next(os.walk(args.indir))[1]
    for sample in samples:
        sdir = os.path.join(args.indir, sample)
        with open(os.path.join(sdir, sample + ".ratt_report.tsv"), "w") as ratt_out:
            header, rows = list(), list()
            for d in sorted(next(os.walk(sdir))[1]):
                try:
                    f = glob.glob(os.path.join(sdir, d, "*.final.gff.delta"))[0]
                except:
                    print("No delta file in {}.".format(os.path.join(sdir, d)), file=sys.stderr)
                    continue
                has_transfer_issue = os.stat(f).st_size > 0
                
                try:
                    f = glob.glob(os.path.join(sdir, d, "*.final.gff"))[0]
                except:
                    print("No .final.gff in {}.".format(os.path.join(sdir, d)), file=sys.stderr)
                    continue
                with open(f) as fin:
                    trfcounts, trfgenecounts = countGFFFeatures(fin)

                ref = d.replace(sample + ".", "")
                print("SAMPLE=", sample, "REF=", ref)
                row = [ref, has_transfer_issue]
                if not header:
                    header = ANN_REPORT_HEADER
                    #header = ["Reference"]
                    #header.extend(map(lambda x:x+"[ref]", FEATURES))
                    #header.extend(map(lambda x:x+"[transferred]", FEATURES))
                    #header.extend(map(lambda x:x+"[%]", FEATURES))
                    #header.extend(["Total[ref]", "Total[transferred]", "Total[%]"])
                    print(*header, sep="\t", file=ratt_out)
                     
                row.append(len(refgenecounts.get(ref, set())))
                row.append(len(trfgenecounts))
                row.append(row[-1]/row[-2] if row[-2] else "NA")

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

            for row in sorted(rows, key=lambda x:(x[-2], x[-1]) if (x[-2] != "NA" and x[-1] != "NA") else -1, reverse=True):            
                print(*row, sep="\t", file=ratt_out)

    return True  




if __name__ == "__main__":
    main()
