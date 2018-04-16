#!/usr/bin/env python
import sys
import os
from os.path import join, basename
import glob
import re
import csv
import argparse
import pathlib

from collections import Counter, namedtuple, OrderedDict

from bgrrl.samplesheet import readSamplesheet, Sample
from bgrrl.enterobase_helpers import validateEnterobaseInput, ECriteria, loadEnterobaseCriteria

ENTERO_CRITERIA = dict()

BlobSample = namedtuple("BlobSample", "sample ncontigs dom_org dom_org_ncontigs dom_org_perc dom_org_span subdom_org subdom_org_ncontigs subdom_org_perc subdom_org_span".split(" "))

def ENTERO_FILTER(_in, organism="Salmonella", out=sys.stdout, full_out=None):
    header = ""
    header_written = False
    crit = ENTERO_CRITERIA.get(organism, None)
    print("Running ENTERO_FILTER for organism: {}".format(organism), crit)
    for row in _in:
        if not header:
            header = list(row)
            if crit is not None:
                header.extend(["Enterobase category", "Meets Enterobase criteria?", "Genome Size check", "Contig amount", "N50 check", "N content"])
        elif crit is None:
            if not header_written:
                header_written = True
                print(*header, sep="\t", file=out)
                if full_out is not None:
                    print(*header, sep="\t", file=full_out)
  
            print(*row, sep="\t", file=out)
            if full_out is not None:
                print(*row, sep="\t", file=full_out)
            yield row[0]
        else:
            genome_size_ok = crit.minsize <= int(row[8]) <= crit.maxsize
            contig_amt_ok = int(row[2]) < crit.ncontigs
            n50_ok = int(row[17]) > crit.n50
            ncontent_ok = float(row[21]) < crit.ncount
            meets_criteria = all((genome_size_ok, contig_amt_ok, n50_ok, ncontent_ok))
            if not header_written:
                header_written = True
                print(*header, sep="\t", file=out)
                if full_out is not None:
                    print(*header, sep="\t", file=full_out)
            print(*row, organism if meets_criteria else "", meets_criteria, genome_size_ok, contig_amt_ok, n50_ok, ncontent_ok, sep="\t", file=out)
            if full_out is not None:
                print(*row, organism if meets_criteria else "", meets_criteria, genome_size_ok, contig_amt_ok, n50_ok, ncontent_ok, sep="\t", file=full_out)
            if meets_criteria:
                yield row[0]

def TAX_FILTER(blob_data, organism="Salmonella", out=sys.stdout, full_out=None):
    crit = ENTERO_CRITERIA.get(organism, None)
    print("Running TAX_FILTER for organism: {}".format(organism), crit)
    lastColH = "MeetsEnterobaseCriteria? (%predom. >= {:.3f})".format(crit.spcount) if crit is not None else ""
    header, full_header = False, False
    for sample in blob_data:
        if sample.dom_org != organism:
            continue
        meets_enterobase = crit is None or (sample.dom_org_perc is not None and float(sample.dom_org_perc) >= crit.spcount)
        lastCol = int(meets_enterobase) if crit is not None else ""
        if full_out is not None:                
            if not full_header:
                full_header = True
                print("Sample", "#contigs", "Predominant taxon", "#contigs(Predominant taxon)", "%(Predominant taxon)", "span(Predominant taxon)", "Subdominant taxon", "#contigs(Subdominant taxon)", "%(Subdominant taxon)", "span(Subdominant taxon)", lastColH, sep="\t", file=full_out)
            print(*sample, lastCol, sep="\t", file=full_out)

        if meets_enterobase:
            if not header:
                header = True
                print("Sample", "#contigs", "Predominant taxon", "#contigs(Predominant taxon)", "%(Predominant taxon)", "span(Predominant taxon)", "Subdominant taxon", "#contigs(Subdominant taxon)", "%(Subdominant taxon)", "span(Subdominant taxon)", lastColH, sep="\t", file=out)
                
            print(*sample, lastCol, sep="\t", file=out)
            yield sample.sample


def main(args_in=sys.argv[1:]):
    ap = argparse.ArgumentParser()
    ap.add_argument("indir", type=str, default=".")
    ap.add_argument("enterobase_groups", type=str, default="")
    ap.add_argument("entero_criteria", type=str, default="")
    ap.add_argument("--report-dir", type=str, default="")
    ap.add_argument("--mode", "-m", type=str, default="asm") # asm/survey
    
    args = ap.parse_args(args_in)

    global ENTERO_CRITERIA
    ENTERO_CRITERIA = loadEnterobaseCriteria(args.entero_criteria)
    
    print("Running asm:report generation...") #, end="", flush=True)
    if not args.report_dir:
        report_dir = join(args.indir, "reports")
    else:
        report_dir = join(args.report_dir)
    pathlib.Path(report_dir).mkdir(parents=True, exist_ok=True)

    print("Reading global QUAST report...", end="", flush=True)
    try:
        with open(join(report_dir, "quast_report.tsv")) as quast_in:
            quast_report = list(row for row in csv.reader(quast_in, delimiter="\t"))
    except:
        print()
        print("Error: No QUAST report found at {}. Exiting.".format(report_dir))
        sys.exit(1)
        pass
    print(" Done")       

    print("Reading global Blobtools report...", end="", flush=True)
    try:
        with open(join(report_dir, "blobtools_report.tsv")) as blob_in:
            blob_report = list(BlobSample(*row) for row in csv.reader(blob_in, delimiter="\t"))
    except:
        print()
        print("Error: No Blobtools report found at {}. Exiting.".format(report_dir))
        sys.exit(1)
    print(" Done")
    valid_sample_groups = validateEnterobaseInput(args.enterobase_groups, ENTERO_CRITERIA)

    if not valid_sample_groups:
        print("No Enterobase groups specified. Proceeding without checking Enterobase criteria.")
        with open(join(report_dir, "all_samples.txt"), "w") as samples_out:
            print(quast_report)
            print(*sorted(line[0] for line in quast_report[1:]), sep="\n", file=samples_out)
    else:
        print("Checking Enterobase groups...") 
        
        with open(join(report_dir, "all_quast_taxonomy_report.tsv"), "w") as full_out, open(join(report_dir, "eb_taxonomy_report.tsv"), "w") as full_tax_out:
            for sgroup in valid_sample_groups:
                print(sgroup)
                eb_asm_out_fn = join(report_dir, "eb_{}_quast_report.tsv".format(sgroup))
                eb_samples_out_fn = join(report_dir, "eb_{}_samples.txt".format(sgroup))
                eb_tax_out_fn = join(report_dir, "eb_{}_taxonomy_report.tsv".format(sgroup))
                with open(eb_asm_out_fn, "w") as eb_asm_out, open(eb_samples_out_fn, "w") as eb_samples_out, open(eb_tax_out_fn, "w") as eb_tax_out:
                    eb_pass_taxonomy = set(TAX_FILTER(blob_report, organism=sgroup, out=eb_tax_out, full_out=full_tax_out))
                    quast_stream = list(filter(lambda x:x[0] in eb_pass_taxonomy or x[0].startswith("Assembly"), quast_report))
                    eb_pass_assembly = set(ENTERO_FILTER(quast_stream, organism=sgroup, out=eb_asm_out, full_out=full_out))
                    all_passed = eb_pass_assembly.intersection(eb_pass_taxonomy)
                    if all_passed:
                        print(*sorted(all_passed), sep="\n", file=eb_samples_out)

                    blob_report = list(filter(lambda x:x.sample not in eb_pass_taxonomy, blob_report))

    print(" Done.\n Generated asm reports in {}.".format(report_dir))
    return True


if __name__ == "__main__":
    main()
