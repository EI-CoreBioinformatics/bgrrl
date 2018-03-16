#!/usr/bin/env python
import sys
import os
import glob
import csv
import argparse
from collections import namedtuple, Counter
from os.path import exists
import yaml
# https://stackoverflow.com/a/600612
import pathlib

from bgrrl import run_snakemake

Sample = namedtuple("Sample", "sampleID customerSampleID R1 R2 S taxonomyID taxonomyTxt fastqcR1 fastqcR2 fastqcS".split(" "))

def validateEnterobaseInput(eb_input_str, eb_criteria):
    print("Determining downstream processing mode...", end="", flush=True)
    if not eb_input_str:
        print(" DEFAULT")
        valid_sample_groups = list()
    else:
        if eb_input_str == "all":
            print(" EB:all")
            check_sample_groups = list(eb_criteria.keys())
        else:
            check_sample_groups = eb_input_str.strip().split(",")
       
        valid_sample_groups = list(filter(lambda g:g in eb_criteria, check_sample_groups))
        invalid_sample_groups = set(check_sample_groups).difference(valid_sample_groups)

        if not valid_sample_groups:
            print(" DEFAULT:FALLBACK")
            print("No valid Enterobase groups found. Proceeding without checking Enterobase criteria.")
            pass
        elif not eb_input_str == "all":
            print(" EB")
            print("Moving forward with Enterobase groups: {}.".format(",".join(sorted(valid_sample_groups))))
        if invalid_sample_groups:
            print("Found invalid Enterobase groups: {}.".format(",".join(sorted(invalid_sample_groups))))

    return valid_sample_groups



def compileBUSCO(busco_dir, out=sys.stdout):
    print("Sample", "Mode", "Complete BUSCOs (C)", "Complete and single-copy BUSCOs (S)", "Complete and duplicated BUSCOs (D)", "Fragmented BUSCOs (F)", "Missing BUSCOs (M)", "Total BUSCO groups searched", "C%", "S%", "D%", "F%", "M%", sep="\t", file=out)
    for cdir, dirs, files in os.walk(busco_dir):
        summary = list(filter(lambda s:"short_summary" in s, files))
        if summary:
            with open(os.path.join(cdir, summary[0])) as _in:
                for line in _in:
                    if line.startswith('# BUSCO was run in mode:'):
                        mode = re.search(' (transcriptome|proteins|genome)', line).group().strip()
                        break
                [next(_in), next(_in), next(_in)]
                numbers = [int(line.strip().split()[0]) for line in _in]
                percentages = list(map(lambda x:"{:.1f}".format(int(x)/int(numbers[-1] * 100)), numbers[:-1]))
                print("_".join(summary[0].split("_")[1:3]), mode, *numbers, *percentages, sep="\t", file=out)

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
    print("Sample", "Organism", "nC=nContigs", "nCO=nContigs[Organism]", "fO=nCO/nC", "MeetsEnterobaseCriteria? (fO >= {:.3f})".format(crit.spcount), sep="\t", file=out)
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
            print(sample, organism, sum(taxcounter.values()), orgcount, "{:.3f}".format(orgfrac), int(meets_enterobase), sep="\t", file=out)
            if meets_enterobase:
                yield sample

def compileQUASTReport(quast_dir, out=sys.stdout):
    header = ""
    for cdir, dirs, files in os.walk(quast_dir):
        if "transposed_report.tsv" in files:
            with open(os.path.join(cdir, "transposed_report.tsv")) as qin:
                first = next(qin)
                if not header:
                    header = first
                    print(header, file=out, end="")
                    yield header
                for row in qin:
                    print(row, file=out, end="")
                    yield row

ECriteria = namedtuple("ECriteria", "minsize maxsize n50 ncontigs ncount spcount".split(" "))
# https://bitbucket.org/enterobase/enterobase-web/wiki/EnteroBase%20Backend%20Pipeline%3A%20QA%20evaluation
ENTERO_CRITERIA = { "Salmonella": ECriteria(4000000, 5800000, 20000, 600, 0.03, 0.7),
                    "Escherichia": ECriteria(3700000, 6400000, 20000, 600, 0.03, 0.7),
                    "Shigella": ECriteria(3700000, 6400000, 20000, 600, 0.03, 0.7),
                    "Yersinia": ECriteria(3700000, 5500000, 15000, 600, 0.03, 0.65),
                    "Moraxella": ECriteria(1800000, 2600000, 20000, 600, 0.03, 0.65)
 }

def ENTERO_FILTER(_in, organism="Salmonella", out=sys.stdout):
    header = ""
    crit = ENTERO_CRITERIA.get(organism, None)
    print("Running ENTERO_FILTER for organism: {}".format(organism), crit)
    for row in csv.reader(_in, delimiter="\t"):
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
            
            

        



def readSamplesheet(fn, delimiter=","):
    with open(fn) as fi:
        for row in csv.reader(fi, delimiter=delimiter):
            # print(row)
            sample = Sample(*row)
            # print(sample)
            # print((sample.R1 and exists(sample.R1) and sample.R2 and exists(sample.R2)))
            assert sample.sampleID
            # assert (sample.R1 and exists(sample.R1) and sample.R2 and exists(sample.R2)) or (sample.S and exists(sample.S))
            if not sample.customerSampleID:
                row[1] = sample.sampleID
            yield (sample.sampleID, Sample(*row))

def verifySamplesheet(fn, delimiter=","):
    for sample_id, sample in readSamplesheet(fn, delimiter=delimiter):
        r1_exists, r2_exists, s_exists = map(os.path.exists, (sample.R1, sample.R2, sample.S))
        if not (r1_exists and r2_exists):
            raise ValueError("Cannot find R1/R2 data at R1={}, R2={}.".format(sample.R1, sample.R2))
        return True

def run_fin(samplesheet, out_dir, args, exe_env, bgrrl_config=dict()):
    config = bgrrl_config
    if verifySamplesheet(samplesheet):
        config["samplesheet"] = samplesheet
    config["out_dir"] = out_dir
    config["package_dir"] = os.path.join(os.path.dirname(out_dir), "Data_Package")
    config["project_prefix"] = args.project_prefix

    config["enterobase_groups"] = validateEnterobaseInput(args.enterobase_groups, ENTERO_CRITERIA) if "enterobase_groups" in args else list()

    if "fin_mode" in args:
        config["finalize_mode"] = args.fin_mode


    config_file = os.path.join(out_dir, "bg-fin.conf.xml")
    with open(config_file, "w") as outfile:
        yaml.dump(config, outfile, default_flow_style=False)

    print("Running BG-FIN")
    res = run_snakemake(os.path.join(os.path.dirname(__file__), "zzz", "bgrrl-fin.smk.py"), out_dir, config_file, exe_env, dryrun=False, unlock=args.unlock)
    return res


def run_qc(samplesheet, out_dir, args, exe_env, bgrrl_config=dict()):
    print(bgrrl_config)
    config = bgrrl_config # dict()
    if verifySamplesheet(samplesheet):
        config["samplesheet"] = samplesheet
    config["out_dir"] = out_dir
    config["no_normalization"] = args.no_normalization
    
    config_file = os.path.join(out_dir, "bg-qc.conf.xml")
    with open(config_file, "w") as outfile:
        yaml.dump(config, outfile, default_flow_style=False)

    print("Running BG-QC") 
    res = run_snakemake(os.path.join(os.path.dirname(__file__), "zzz", "bgrrl-qc.smk.py"), out_dir, config_file, exe_env, dryrun=False, unlock=args.unlock) 

  
    # write samplesheet into configfile
    # test if sample data exists
    # run bg-qc on samples
    # upon return, parse bg-qc's output and annotate samples
    # generate samplesheet for bg-asm
    return res

def run_asm(samplesheet, out_dir, args, exe_env, bgrrl_config=dict()):
    print(bgrrl_config)
    config = bgrrl_config
    if verifySamplesheet(samplesheet):
        config["samplesheet"] = samplesheet
    config["out_dir"] = out_dir
    config["etc"] = os.path.join(os.path.dirname(__file__), "..", "etc")
    config["cwd"] = os.getcwd()
    config["assembler"] = args.assembler
    config["reapr_correction"] = False
    # print(args)
    config["no_normalization"] = args.no_normalization 

    if args.contig_minlen:
        config["use_asm_lengthfilter"] = True
        config["asm_lengthfilter_contig_minlen"] = args.contig_minlen
    else:
        config["use_asm_lengthfilter"] = False
        config["asm_lengthfilter_contig_minlen"] = 0
        
        # print(config)

    config_file = os.path.join(out_dir, "bg-asm.conf.xml")
    with open(config_file, "w") as outfile:
        yaml.dump(config, outfile, default_flow_style=False)

    print("Running BG-ASM")
    # print("THIS:", __file__)
    res = run_snakemake(os.path.join(os.path.dirname(__file__), "zzz", "bgrrl-asm.smk.py"), out_dir, config_file, exe_env, dryrun=False, unlock=args.unlock)
    return res

def run_ann(samplesheet, out_dir, args, exe_env, bgrrl_config=dict()):
    print(bgrrl_config)
    config = bgrrl_config
    if verifySamplesheet(samplesheet):
        config["samplesheet"] = samplesheet
    config["out_dir"] = out_dir
    config["etc"] = os.path.join(os.path.dirname(__file__), "..", "etc")
    config["cwd"] = os.getcwd()

    config["run_ratt"] = args.annotation in ("both", "ratt")
    config["run_prokka"] = args.annotation in ("both", "prokka")
    config["ratt_reference"] = args.ratt_reference_dir

    assert not config["run_ratt"] or os.path.exists(config["ratt_reference"]), "Missing reference data for ratt. Please make sure to use the --ratt-reference-dir parameter."

    config_file = os.path.join(out_dir, "bg-ann.conf.xml")
    with open(config_file, "w") as outfile:
        yaml.dump(config, outfile, default_flow_style=False)

    print("Running BG-ANN")
    res = run_snakemake(os.path.join(os.path.dirname(__file__), "zzz", "bgrrl-ann.smk.py"), out_dir, config_file, exe_env, dryrun=False, unlock=args.unlock)
    return res
    





if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('samplesheet', type=str)

    args = ap.parse_args()

    assert args.samplesheet and os.path.exists(args.samplesheet)
    samples = dict(readSamplesheet(args.samplesheet))

    for sample in samples:
        print(*samples[sample])
