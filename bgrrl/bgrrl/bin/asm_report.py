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

from bgrrl.bgrrl import readSamplesheet, Sample, validateEnterobaseInput

ECriteria = namedtuple("ECriteria", "minsize maxsize n50 ncontigs ncount spcount".split(" "))
BlobSample = namedtuple("BlobSample", "sample ncontigs dom_org dom_org_ncontigs dom_org_perc dom_org_span subdom_org subdom_org_ncontigs subdom_org_perc subdom_org_span".split(" "))

# https://bitbucket.org/enterobase/enterobase-web/wiki/EnteroBase%20Backend%20Pipeline%3A%20QA%20evaluation
ENTERO_CRITERIA = { "Salmonella": ECriteria(4000000, 5800000, 20000, 600, 0.03, 0.7),
                    "Escherichia": ECriteria(3700000, 6400000, 20000, 600, 0.03, 0.7),
                    "Shigella": ECriteria(3700000, 6400000, 20000, 600, 0.03, 0.7),
                    "Yersinia": ECriteria(3700000, 5500000, 15000, 600, 0.03, 0.65),
                    "Moraxella": ECriteria(1800000, 2600000, 20000, 600, 0.03, 0.65)
 }

ASSEMBLY_STAGES = OrderedDict([("asm_main_ucn", "Main,Unicycler,normalized"),
                               ("asm_fb1_uct", "Fallback1,Unicycler,trimmed"),
                               ("asm_fb2_spn", "Fallback2,Spades,normalized"),
                               ("asm_fb3_spt", "Fallback3,Spades,trimmed"),
                               ("asm_main_ven", "Main,Velvet,normalized"),
                               ("NA", "not_assembled")])

"""
def compileASMInfo(asm_dir, out=sys.stdout, asm_stat_out=sys.stdout):
    asm_tag_ctr = Counter()
    for cdir, dirs, files in os.walk(asm_dir):
        sample = basename(cdir)
        if sample != "log" and os.path.dirname(cdir) == asm_dir:
            asm_tag = ([f for f in files if f.startswith("asm_")] + ["NA"])[0]
            asm_tag_ctr[asm_tag] += 1 # .setdefault(asm_tag, list()).append(cdir)
            print(sample, asm_tag, sep="\t", file=out)
    print(asm_tag_ctr)
    for asm_tag in ASSEMBLY_STAGES:
        if asm_tag_ctr[asm_tag] > 0:
            print(ASSEMBLY_STAGES[asm_tag], asm_tag_ctr[asm_tag], asm_tag_ctr[asm_tag]/sum(asm_tag_ctr.values()), sep="\t", file=asm_stat_out)
    print("Total", "", "", sum(asm_tag_ctr.values()), sep="\t", file=asm_stat_out)
"""

"""
def compileQUASTReport(quast_dir, out=sys.stdout):
    header = ""
    for cdir, dirs, files in os.walk(quast_dir):
        if "transposed_report.tsv" in files:
            with open(join(cdir, "transposed_report.tsv")) as qin:
                r = csv.reader(qin, delimiter="\t")
                first = next(r)
                pif not header:
                    header = first
                    print(*header, file=out, sep="\t")
                    yield header
                for row in r:
                    print(*row, sep="\t", file=out)
                    yield row
"""

def ENTERO_FILTER(_in, organism="Salmonella", out=sys.stdout, full_out=None):
    header = ""
    header_written = False
    crit = ENTERO_CRITERIA.get(organism, None)
    print("Running ENTERO_FILTER for organism: {}".format(organism), crit)
    # for row in csv.reader(_in, delimiter="\t"):
    for row in _in:
        if not header:
            header = list(row)
            if crit is not None:
                header.extend(["Enterobase category", "Meets Enterobase criteria?", "Genome Size check", "Contig amount", "N50 check", "N content"])
            # print(*header, sep="\t", file=out)
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

            #checks_ok = genome_size_ok << 3 | contig_amt_ok << 2 | n50_ok << 1 | ncontent_ok
            # print(row)
            # if crit.minsize <= int(row[8]) <= crit.maxsize:
            #    if int(row[2]) < crit.ncontigs:
            #        if int(row[17]) > crit.n50:
            #            if float(row[21]) < crit.ncount:
            #                print(*row, sep="\t", file=out)
            #                yield row[0]


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

# def BLOBREPORT(blob_dir, out=sys.stdout):    
#     print("Running BLOBREPORT... " + blob_dir)
#     print("Sample", "#contigs", "Predominant genus", "#contigs(Predominant genus)", "%(Predominant genus)", "span(Predominant genus)[bp]", "Subdominant genus", "#contigs(Subdominant genus)", "%(Subdominant genus)", "span(Subdominant genus)[bp]", sep="\t", file=out)
#     for cdir, dirs, files in os.walk(blob_dir):
#         blobtable = list(filter(lambda s:s.endswith(".blobDB.table.txt"), files))
#         if blobtable:
#             sample = blobtable[0].split(".")[0]
#             taxcounter, spancounter, taxmap = Counter(), Counter(), dict()
#             with open(join(cdir, blobtable[0])) as tin:
#                 for row in csv.reader(tin, delimiter="\t"):
#                     if not row[0].startswith("#"):
#                         genus = row[5].split(" ")[0]
#                         taxcounter[genus] += 1
#                         spancounter[genus] += int(row[1])
#                         taxmap.setdefault(genus, Counter())[row[5]] += 1
#             orgs = sorted(taxcounter.items(), key=lambda x:x[1], reverse=True)
#             ncontigs = sum(taxcounter.values())
#             dom_org, vdom_org = orgs.pop(0), (None, 0)
#             if orgs:
#                 vdom_org = orgs.pop(0)
#             blob_data = BlobSample(sample, ncontigs, dom_org[0], dom_org[1], dom_org[1]/ncontigs if ncontigs else None, spancounter[dom_org[0]], vdom_org[0], vdom_org[1], vdom_org[1]/ncontigs if ncontigs else None, spancounter[vdom_org[0]])
#             print(*blob_data, sep="\t", file=out)
#             yield blob_data

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
    ap.add_argument("--report-dir", type=str, default="")
    ap.add_argument("--mode", "-m", type=str, default="asm") # asm/survey
    args = ap.parse_args(args_in)
 
    
    print("Running asm:report generation...") #, end="", flush=True)
    if not args.report_dir:
        report_dir = join(args.indir, "reports")
    else:
        report_dir = join(args.report_dir)
    pathlib.Path(report_dir).mkdir(parents=True, exist_ok=True)

    """
    print("Gathering assembly stage information...")
    with open(join(report_dir, "assembly_stage_summary.tsv"), "w") as asm_stage_summary, open(join(report_dir, "assembly_stages.tsv"), "w") as asm_stages:
        compileASMInfo(join(args.indir, "assembly"), out=asm_stages, asm_stat_out=asm_stage_summary)
    """

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

    # print("Determining downstream processing mode...", end="", flush=True)
    # if args.enterobase_groups:
    #     if args.enterobase_groups == "all":
    #        print(" EB:all")
    #        check_sample_groups = list(ENTERO_CRITERIA.keys())
    #    else:
    #        check_sample_groups = args.enterobase_groups.split(",")
    #    valid_sample_groups = list(filter(lambda g:g in ENTERO_CRITERIA, check_sample_groups))        
    #    invalid_sample_groups = set(check_sample_groups).difference(valid_sample_groups)
    #    if not valid_sample_groups:
    #        print(" DEFAULT:FALLBACK")
    #        print("No valid Enterobase groups found. Proceeding without checking Enterobase criteria.")
    #        pass
    #    elif valid_sample_groups:
    #        if not args.enterobase_groups == "all":
    #            print(" EB")
    #        print("Moving forward with Enterobase groups: {}.".format(",".join(sorted(valid_sample_groups))))
    #    if invalid_sample_groups:
    #        print("Found invalid Enterobase groups: {}.".format(",".join(sorted(invalid_sample_groups)))) 
    #else:
    #    print(" DEFAULT")
    #    valid_sample_groups = list()
    
    # blob_path = join(args.indir, "qa", args.mode, "blobtools", "blob")
    valid_sample_groups = validateEnterobaseInput(args.enterobase_groups, ENTERO_CRITERIA)

    if not valid_sample_groups:
        print("No Enterobase groups specified. Proceeding without checking Enterobase criteria.")
        with open(join(report_dir, "all_samples.txt"), "w") as samples_out:
            print(quast_report)
            print(*sorted(line[0] for line in quast_report[1:]), sep="\n", file=samples_out)
        # with open(join(report_dir, "all_taxonomy_report.tsv"), "w") as tax_out:
        #    # list(TAX_FILTER(blob_path, organism="", out=tax_out))
        #    list(BLOBREPORT(blob_path, out=tax_out))
        #    """try:
        #    except:
        #        print("Error: No BLOBTOOLS data found at {}. Exiting.".format(blob_path))
        #        sys.exit(1)
        #        pass
        #    """
    else:
        print("Checking Enterobase groups...") 
        #with open(join(report_dir, "all_taxonomy_report.tsv"), "w") as tax_out:
        #    blob_report = list(BLOBREPORT(blob_path, out=tax_out))
        

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
                    # try:
                    #    eb_pass_taxonomy = set(TAX_FILTER(blob_path, organism=sgroup, out=eb_tax_out))
                    #except:
                    #    print("Error: No BLOBTOOLS data found at {}. Exiting.".format(blob_path))
                    #    sys.exit(1)
                    #    pass
                    all_passed = eb_pass_assembly.intersection(eb_pass_taxonomy)
                    if all_passed:
                        print(*sorted(all_passed), sep="\n", file=eb_samples_out)

                    blob_report = list(filter(lambda x:x.sample not in eb_pass_taxonomy, blob_report))

    print(" Done.\n Generated asm reports in {}.".format(report_dir))
    return True


if __name__ == "__main__":
    main()
