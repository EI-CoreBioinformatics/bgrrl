#!/usr/bin/env python
import sys
import os
import glob
import re
import csv
import argparse
import pathlib

from collections import Counter, namedtuple

from bgrrl.bgrrl import readSamplesheet, Sample

ECriteria = namedtuple("ECriteria", "minsize maxsize n50 ncontigs ncount spcount".split(" "))
# https://bitbucket.org/enterobase/enterobase-web/wiki/EnteroBase%20Backend%20Pipeline%3A%20QA%20evaluation
ENTERO_CRITERIA = { "Salmonella": ECriteria(4000000, 5800000, 20000, 600, 0.03, 0.7),
                    "Escherichia": ECriteria(3700000, 6400000, 20000, 600, 0.03, 0.7),
                    "Shigella": ECriteria(3700000, 6400000, 20000, 600, 0.03, 0.7),
                    "Yersinia": ECriteria(3700000, 5500000, 15000, 600, 0.03, 0.65),
                    "Moraxella": ECriteria(1800000, 2600000, 20000, 600, 0.03, 0.65)
 }


def compileQUASTReport(quast_dir, out=sys.stdout):
    header = ""
    for cdir, dirs, files in os.walk(quast_dir):
        if "transposed_report.tsv" in files:
            with open(os.path.join(cdir, "transposed_report.tsv")) as qin:
                r = csv.reader(qin, delimiter="\t")
                first = next(r)
                if not header:
                    header = first
                    print(*header, file=out, sep="\t")
                    yield header
                for row in r:
                    print(*row, sep="\t", file=out)
                    yield row

def ENTERO_FILTER(_in, organism="Salmonella", out=sys.stdout):
    header = ""
    crit = ENTERO_CRITERIA.get(organism, None)
    print("Running ENTERO_FILTER for organism: {}".format(organism), crit)
    # for row in csv.reader(_in, delimiter="\t"):
    for row in _in:
        if not header:
            header = row
            print(*header, sep="\t", file=out)
        elif crit is None:
            print(*row, sep="\t", file=out)
            yield row[0]
        else:
            # print(row)
            if crit.minsize <= int(row[8]) <= crit.maxsize:
                if int(row[2]) < crit.ncontigs:
                    if int(row[17]) > crit.n50:
                        if float(row[21]) < crit.ncount:
                            print(*row, sep="\t", file=out)
                            yield row[0]


"""
parse quast reports and check Enterobase criteria
1. obtain header from one sample
head -n1 Analysis/qa/quast/FD01543398_L006/transposed_report.tsv > quast_report.tsv
2. get data from all samples' transposed reports
find Analysis/qa/quast -name 'transposed_report.tsv' -exec tail -n +2 {} \; >> quast_report.tsv
3. Filter by S.enterica criteria
awk -v FS="\t" -v OFS="\t" '/^Assembly/ {print $0; next;} {if (4000000 <= $9 && $9 <= 5800000 && $3 < 600 && $18 > 20000 && $22 < 0.03) print $0;}' quast_report.tsv > quast_report.enterobase.tsv

check blobtools tables for taxonomic classification
"""

def TAX_FILTER(blob_dir, organism="Salmonella", out=sys.stdout):
    crit = ENTERO_CRITERIA.get(organism, None)
    print("Running TAX_FILTER for organism: {}".format(organism), crit)
    lastCol = "MeetsEnterobaseCriteria? (fO >= {:.3f})".format(crit.spcount) if crit is not None else ""
    print("Sample", "Organism", "nC=nContigs", "nCO=nContigs[Organism]", "fO=nCO/nC", lastCol, sep="\t", file=out)
    for cdir, dirs, files in os.walk(blob_dir):
        blobtable = list(filter(lambda s:s.endswith(".blobDB.table.txt"), files))
        if blobtable:
            taxcounter = Counter()
            with open(os.path.join(cdir, blobtable[0])) as tin:
                for row in csv.reader(tin, delimiter="\t"):
                    if not row[0].startswith("#"):
                        taxcounter[row[5].split(" ")[0]] += 1
            sample = blobtable[0].split(".")[0]
            orgcount = sum(taxcounter[org] for org in taxcounter if org.startswith(organism))
            orgfrac = orgcount/sum(taxcounter.values())
            meets_enterobase = crit is None or orgfrac >= crit.spcount
            dominant_org = [org for org in taxcounter if taxcounter[org] == max(taxcounter.values())][0] 
            if not organism or dominant_org.startswith(organism):
                lastCol = int(meets_enterobase) if crit is not None else ""
                print(sample, dominant_org, sum(taxcounter.values()), orgcount, "{:.3f}".format(orgfrac), lastCol, sep="\t", file=out)
            if meets_enterobase:
                yield sample


def main(args_in=sys.argv[1:]):
    ap = argparse.ArgumentParser()
    ap.add_argument("indir", type=str, default=".")
    ap.add_argument("enterobase_groups", type=str, default="")
    ap.add_argument("--report-dir", type=str, default="")
    args = ap.parse_args(args_in)
 
    
    print("Running asm:report generation...") #, end="", flush=True)
    if not args.report_dir:
        report_dir = os.path.join(args.indir, "reports")
    else:
        report_dir = os.path.join(args.report_dir)
    pathlib.Path(report_dir).mkdir(parents=True, exist_ok=True)

    print("Reading global QUAST report...", end="", flush=True)
    try:
        with open(os.path.join(report_dir, "quast_report.tsv")) as quast_in:
            quast_report = list(row for row in csv.reader(quast_in, delimiter="\t"))
    except:
        print()
        print("Error: No QUAST report found at {}. Exiting.".format(report_dir))
        sys.exit(1)
        pass
    print(" Done")        

    print("Determining downstream processing mode...", end="", flush=True)
    if args.enterobase_groups:
        if args.enterobase_groups == "all":
            print(" EB:all")
            check_sample_groups = list(ENTERO_CRITERIA.keys())
        else:
            check_sample_groups = args.enterobase_groups.split(",")
        valid_sample_groups = list(filter(lambda g:g in ENTERO_CRITERIA, check_sample_groups))        
        invalid_sample_groups = set(check_sample_groups).difference(valid_sample_groups)
        if not valid_sample_groups:
            print(" DEFAULT:FALLBACK")
            print("No valid Enterobase groups found. Proceeding without checking Enterobase criteria.")
            pass
        elif valid_sample_groups:
            if not args.enterobase_groups == "all":
                print(" EB")
            print("Moving forward with Enterobase groups: {}.".format(",".join(sorted(valid_sample_groups))))
        if invalid_sample_groups:
            print("Found invalid Enterobase groups: {}.".format(",".join(sorted(invalid_sample_groups)))) 
    else:
        print(" DEFAULT")
        valid_sample_groups = list()
    
    blob_path = os.path.join(args.indir, "qa", "blobtools", "blob")

    if not valid_sample_groups:
        print("No Enterobase groups specified. Proceeding without checking Enterobase criteria.")
        with open(os.path.join(report_dir, "all_samples.txt"), "w") as samples_out:
            print(quast_report)
            print(*sorted(line[0] for line in quast_report[1:]), sep="\n", file=samples_out)
        with open(os.path.join(report_dir, "all_taxonomy_report.tsv"), "w") as tax_out:
            list(TAX_FILTER(blob_path, organism="", out=tax_out))
            """try:
            except:
                print("Error: No BLOBTOOLS data found at {}. Exiting.".format(blob_path))
                sys.exit(1)
                pass
            """
    else:
        print("Checking Enterobase groups...") 
        for sgroup in valid_sample_groups:
            print(sgroup)
            with open(os.path.join(report_dir, "eb_{}_quast_report.tsv".format(sgroup)), "w") as eb_asm_out, open(os.path.join(report_dir, "eb_{}_taxonomy_report.tsv".format(sgroup)), "w") as eb_tax_out, open(os.path.join(report_dir, "eb_{}_samples.txt".format(sgroup)), "w") as eb_samples_out:
                eb_pass_assembly = set(ENTERO_FILTER(quast_report, organism=sgroup, out=eb_asm_out))
                try:
                    eb_pass_taxonomy = set(TAX_FILTER(blob_path, organism=sgroup, out=eb_tax_out))
                except:
                    print("Error: No BLOBTOOLS data found at {}. Exiting.".format(blob_path))
                    sys.exit(1)
                    pass
                print(*sorted(eb_pass_assembly.intersection(eb_pass_taxonomy)), sep="\n", file=eb_samples_out)

    print(" Done.\n Generated asm reports in {}.".format(report_dir))



if __name__ == "__main__":
    main()
