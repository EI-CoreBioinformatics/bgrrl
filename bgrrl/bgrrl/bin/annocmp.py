import sys
import csv
import os
from os.path import join, dirname, basename
import subprocess
import pathlib

from collections import namedtuple, OrderedDict, Counter
from argparse import ArgumentParser

from bgrrl.bin.ann_report import ANN_REPORT_HEADER

GFFeature = namedtuple("GFFeature", "seqname source feature start end score strand frame attribute".split(" "))

PROKKA_FEATURES = ["CDS", "rRNA", "tRNA", "tmRNA"]
BEDTOOLS_CMD = "source bedtools-2.26.0_github20170207; bedtools subtract -a {} -b {}"

HEADER = ["Sample", "Best ratt reference"] 
HEADER.extend(ANN_REPORT_HEADER[1:])
HEADER.extend(map(lambda s:s + "(Prokka)", PROKKA_FEATURES))
HEADER.extend(["Novel features", "Novel features (hypothetical)", "Novel features (other)", "Partial features", "Partial features (hypothetical)", "Partial features (other)"])
HEADER.extend(["Unpredicted features"]) #, "Unpredicted features (hypothetical)", "Unpredicted features (other)"])

def readProkkaGFF(_in):
    data = OrderedDict() 
    for row in csv.reader(_in, delimiter="\t"):
        if row and not row[0].startswith("#"):
            attr = OrderedDict((item.split("=")[0], item.split("=")[1]) for item in row[8].split(";"))
            key = attr.get("ID", None)
            if key is not None:
                # print(row, file=sys.stderr)
                data[key] = GFFeature(row[0], row[1], row[2], int(row[3]), int(row[4]), row[5], row[6], row[7], attr)
    return data


def _process_prokka_ratt_delta(prokka_gff, ratt_gff, prokka_full, output_gff):         
    feat_count = Counter()

    pr = subprocess.Popen(BEDTOOLS_CMD.format(prokka_gff, ratt_gff), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)    
    po, pe = pr.communicate()
    prokka_sub = readProkkaGFF(po.decode().strip().split("\n"))

    with open(output_gff, "w") as _out:
        for k in prokka_sub:
            sub_feat = prokka_sub[k]
            sub_len = sub_feat.end - sub_feat.start + 1
            
            try:
                full_feat = prokka_full[k]
            except:
                print("ERROR: MISSING FEATURE FROM PROKKA-ANNOTATION: " + k)
                sys.exit(1)
    
            full_len = full_feat.end - full_feat.start + 1
    
            if sub_len == full_len:
                ftype, fcov = "novel", "0.00%"
            else:
                ftype, fcov = "partial", "{:.2f}%".format((full_len-sub_len)/full_len*100)
            feat_count[ftype] += 1
            feat_count[(sub_feat.attribute.get("product", None), ftype)] += 1
            sub_feat.attribute["feature_type"] = ftype
            sub_feat.attribute["feature_cov"] = fcov
            attribute = ";".join("{}={}".format(k, sub_feat.attribute[k]) for k in sub_feat.attribute)        
    
            print(*(sub_feat[:-1]), attribute, sep="\t", file=_out)

    return feat_count

def _process_ratt_prokka_delta(prokka_gff, ratt_gff, ratt_full, output_gff):
    feat_count = Counter()

    pr = subprocess.Popen(BEDTOOLS_CMD.format(ratt_gff, prokka_gff), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    po, pe = pr.communicate()
    ratt_sub = readProkkaGFF(po.decode().strip().split("\n"))

    with open(output_gff, "w") as _out:
        for k in ratt_sub:
            sub_feat = ratt_sub[k]
            sub_len = sub_feat.end - sub_feat.start + 1
 
            try:
                full_feat = ratt_full[k]
            except:
                print("ERROR: MISSING FEATURE FROM RATT ANNOTATION: " + k)
                sys.exit(1)
     
            full_len = full_feat.end - full_feat.start + 1

            if sub_len == full_len:
                ftype, fcov = "unpredicted", "0.00%"
            else:
                ftype, fcov = "partial", "{:2f}%".format((full_len-sub_len)/full_len*100)
            feat_count[ftype] += 1
            # feat_count[(sub_feat.attribute.get("pseudo
            sub_feat.attribute["feature_type"] = ftype
            sub_feat.attribute["feature_cov"] = fcov
            attribute = ";".join("{}={}".format(k, sub_feat.attribute[k]) for k in sub_feat.attribute)

            print(*(sub_feat[:-1]), attribute, sep="\t", file=_out)

    return feat_count


def compareGFFs(prokka_gff, ratt_tsv, prokka_ratt_delta_gff, ratt_prokka_delta_gff):
    def _getBestRATTTransfer(ratt_tsv):
        with open(ratt_tsv) as _in:
            for i, r in enumerate(csv.reader(_in, delimiter="\t")):
                if i == 1:
                    return r
        return None

    with open(prokka_gff) as _in:
        prokka_full = readProkkaGFF(_in)

    best_ratt_ref = _getBestRATTTransfer(ratt_tsv)
    try:
        refname = best_ratt_ref[0]
    except:
        print("ERROR: RATT ANNOTATION IS EMPTY " + ratt_tsv, file=sys.stderr)
        sys.exit(1)

    prokka_fcount = Counter(feature.feature for feature in prokka_full.values())    
    ratt_gff = join(dirname(ratt_tsv), refname, basename(ratt_tsv).replace(".ratt_report.tsv", "_" + refname + ".final.gff"))

    with open(ratt_gff) as _in:
        ratt_full = readProkkaGFF(_in)

    prokka_ratt_delta = _process_prokka_ratt_delta(prokka_gff, ratt_gff, prokka_full, prokka_ratt_delta_gff)
    ratt_prokka_delta = _process_ratt_prokka_delta(prokka_gff, ratt_gff, ratt_full, ratt_prokka_delta_gff)

    return prokka_ratt_delta, ratt_prokka_delta, best_ratt_ref, prokka_fcount

def main(args=sys.argv[1:]):
    ap = ArgumentParser()
    ap.add_argument("prokka_dir", type=str)
    ap.add_argument("ratt_dir", type=str)
    ap.add_argument("report_dir", type=str)
    args = ap.parse_args(args)

    delta_dir = join(args.prokka_dir, "..", "delta")
    pathlib.Path(delta_dir).mkdir(parents=True, exist_ok=True) 

    samples = next(os.walk(args.prokka_dir))[1]

    try:
        with open(join(report_dir, "blobtools_report.tsv")) as _in:
            tax_report = dict((row[0], row) for i, row in enumerate(csv.reader(_in, delimiter="\t")))
            HEADER.extend(tax_report.get("Sample", list())[1:])
            # "Sample", "#contigs", "Predominant genus", "#contigs(Predominant genus)", "%(Predominant genus)", "span(Predominant genus)[bp]", "Subdominant genus", "#contigs(Subdominant genus)", "%(Subdominant genus)", "span(Subdominant genus)[bp]"
    except:
            tax_report = dict()

    

    with open(join(args.report_dir, "annotation_report.tsv"), "w") as delta_out:
        print(*HEADER, sep="\t", file=delta_out)
        for sample in samples:
            prokka_gff = join(args.prokka_dir, sample, sample + ".prokka.gff")
            ratt_tsv = join(args.ratt_dir, sample, sample + ".ratt_report.tsv")
            prokka_ratt_delta_gff = join(delta_dir, sample + ".prokka_ratt_delta.gff")
            ratt_prokka_delta_gff = join(delta_dir, sample + ".ratt_prokka_delta.gff")
            
            prokka_ratt_delta, ratt_prokka_delta, best_ratt_ref, prokka_fcount = compareGFFs(prokka_gff, ratt_tsv, prokka_ratt_delta_gff, ratt_prokka_delta_gff)

            row = [sample]
            row.extend(best_ratt_ref)
            row.extend(prokka_fcount[feature] for feature in PROKKA_FEATURES)            
            row.extend([prokka_ratt_delta["novel"], 
                        prokka_ratt_delta[("hypothetical protein", "novel")], 
                        prokka_ratt_delta["novel"] - prokka_ratt_delta[("hypothetical protein", "novel")],
                        prokka_ratt_delta["partial"], 
                        prokka_ratt_delta[("hypothetical protein", "partial")], 
                        prokka_ratt_delta["partial"] - prokka_ratt_delta[("hypothetical protein", "partial")],
                        ratt_prokka_delta["unpredicted"]])
            row.extend(tax_report.get(sample, list())[1:])            

            print(*row, sep="\t", file=delta_out)

    return True

if __name__ == "__main__":
    main()
