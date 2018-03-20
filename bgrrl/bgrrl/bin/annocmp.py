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


    

def readProkkaGFF(_in):
    data = dict() 
    for row in csv.reader(_in, delimiter="\t"):
        if row and not row[0].startswith("#"):
            attr = OrderedDict((item.split("=")[0], item.split("=")[1]) for item in row[8].split(";"))
            key = attr.get("ID", None)
            if key is not None:
                # print(row, file=sys.stderr)
                data[key] = GFFeature(row[0], row[1], row[2], int(row[3]), int(row[4]), row[5], row[6], row[7], attr)
    return data
         
def compareGFFs(prokka_gff, ratt_tsv, output_gff):
    with open(prokka_gff) as _in:
        prokka_full = readProkkaGFF(_in)

    with open(ratt_tsv) as _in:
        for i, r in enumerate(csv.reader(_in, delimiter="\t")):
            if i > 1:
                break
            if i == 1:
                refname = r[0]
                best_ratt_ref = r
    # NC_012563.1/PRO1880_Plate3_F1_TACGGCGTT-CTTGTA_NC_012563.1.final.gff
    ratt_gff = os.path.join(os.path.dirname(ratt_tsv), refname, os.path.basename(ratt_tsv).replace(".ratt_report.tsv", "") + "_" + refname + ".final.gff")

    cmd = "source bedtools-2.26.0_github20170207; bedtools subtract -a {} -b {}"
    pr = subprocess.Popen(cmd.format(prokka_gff, ratt_gff), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)    
    po, pe = pr.communicate()
    # print(po)
    # print(pe)
    
    prokka_sub = readProkkaGFF(po.decode().strip().split("\n"))
    # print(prokka_sub)
    feat_count = Counter()

    with open(output_gff, "w") as _out:
        for k in sorted(prokka_sub):
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

    return feat_count, best_ratt_ref

def main(args=sys.argv[1:]):
    ap = ArgumentParser()
    ap.add_argument("prokka_dir", type=str)
    ap.add_argument("ratt_dir", type=str)
    ap.add_argument("report_dir", type=str)
    args = ap.parse_args(args)

    delta_dir = join(args.prokka_dir, "..", "delta")
    pathlib.Path(delta_dir).mkdir(parents=True, exist_ok=True) 

    samples = next(os.walk(args.prokka_dir))[1]

    with open(join(args.report_dir, "annotation_report.tsv"), "w") as delta_out:
        print("Sample", "Best ratt reference", *ANN_REPORT_HEADER[1:], "Novel features", "Novel features (hypothetical)", "Novel features (other)", "Partial features", "Partial features (hypothetical)", "Partial features (other)", sep="\t", file=delta_out)
        for sample in samples:
            prokka_gff = join(args.prokka_dir, sample, sample + ".prokka.gff")
            ratt_tsv = join(args.ratt_dir, sample, sample + ".ratt_report.tsv")
            output_gff = join(delta_dir, sample + ".delta.gff")
            
            feat_count, best_ratt_ref = compareGFFs(prokka_gff, ratt_tsv, output_gff)
            row = [sample]
            row.extend(best_ratt_ref)
            row.extend([feat_count["novel"], 
                        feat_count[("hypothetical protein", "novel")], 
                        feat_count["novel"] - feat_count[("hypothetical protein", "novel")],
                        feat_count["partial"], 
                        feat_count[("hypothetical protein", "partial")], 
                        feat_count["partial"] - feat_count[("hypothetical protein", "partial")],])

            print(*row, sep="\t", file=delta_out)

    return True

if __name__ == "__main__":
    main()
