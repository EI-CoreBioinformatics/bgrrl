import os
from os.path import join 
import sys
import argparse
import csv

def readDatasum(dfile):
    import xlrd
    header, data = list(), dict()
    try:
        book = xlrd.open_workbook(dfile)
    except:
        print("Cannot open datasum file " + dfile, file=sys.stderr)
        exit(1)
    try:    
        sheet = book.sheet_by_name("Sheet1")
    except:
        sheet = book.sheet_by_index(0)
    header = sheet.row_values(0)
    for i in range(1, sheet.nrows):
        data[sheet.row_values(i)[0].replace("_R1", "")] = dict(zip(header, sheet.row_values(i)))
    return data, header


def writeDatasum(data, header, dfile, qc_pass=set(), asm_report=dict(), tax_report=dict(), tax_eb_report=dict(), asm_eb_report=dict()):
    def write_row(sheet, row_index, data, header):
        for j, col in enumerate(header):
            sheet.write(row_index, j, data[col])

    in_pool_cheader = header[-1]
    ip_header, nip_header = list(header[:-1]), list(header[:-1])
    ip_header.append("status")

    asm_header = asm_report.get("Assembly", list())
    asm_eb_header = asm_eb_report.get("Assembly", list())
    tax_header = tax_report.get("Sample", list())
    eb_header = tax_eb_report.get("Sample", "")
    # qc_header = qc_pass.get("Sample", list())
    qc_header = ["Sample", "qc:Status", "fqc:Status", "fqc:Description", "#R1_reads", "#R2_reads", "Read_Length", "kat:Status", "kat:Description", "#K-mer_peaks", "Max_Peak_Volume", "Max_Peak_Volume/total", "GC_peaks", "Est_GenomeSize", "Est_GenomeSize_unit", "K-mer_freq", "mean_GC", "tadpole:Status", "tadpole:Description", "tadpole:Assembly_size"]

    ip_header.extend(qc_header[1:])
    ip_header.extend(asm_header[1:] + asm_eb_header[-6:] + tax_header[2:] + [eb_header])


    import xlsxwriter
    book = xlsxwriter.Workbook(dfile)
    sheet_ip = book.add_worksheet("In Pool")
    sheet_nip = book.add_worksheet("Not in Pool")
    for j, col in enumerate(ip_header):
        sheet_ip.write(0, j, col)
    for j, col in enumerate(nip_header):
        sheet_nip.write(0, j, col)

    row_index_ip, row_index_nip = 1, 1
    for i, k in enumerate(sorted(data), 1):
        if data[k][in_pool_cheader]:
            data[k].update(dict(zip(asm_header[1:], asm_report.get(k, [""]*len(asm_header))[1:])))
            data[k].update(dict(zip(asm_eb_header[-6:], asm_eb_report.get(k, [""]*len(asm_eb_header))[-6:])))
            data[k].update(dict(zip(tax_header[2:], tax_report.get(k, [""]*len(tax_header))[2:])))
            data[k].update(dict(zip([eb_header], [tax_eb_report.get(k, "")])))
            data[k].update(dict(zip(qc_header[1:], qc_pass.get(k, [""]*len(qc_header))[1:])))

            # data[k]["status"] = "assembled" if k in qc_pass else "failed Bioinf_QC"
            data[k]["status"] = "assembled" if qc_pass.get(k, [None, "FAIL"])[1] == "PASS" else "failed Bioinf_QC"
            write_row(sheet_ip, row_index_ip, data[k], ip_header)
            row_index_ip += 1
        else:
            write_row(sheet_nip, row_index_nip, data[k], nip_header)
            row_index_nip += 1

    book.close()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("indir", type=str, help="Path to input directory.")
    ap.add_argument("datasum_file", type=str, help="Path to datasummary file. (xlsx)")
    ap.add_argument("--mode", "-m", type=str, choices=["asm", "survey"], default="asm")
    args = ap.parse_args()

    data, header = readDatasum(args.datasum_file)
    report_dir = join(args.indir, "reports")
    
    # with open(join(report_dir, "samplesheets", "samplesheet.qc_pass.tsv")) as qc_in:
    #     qc_pass = set(row[0] for i, row in enumerate(csv.reader(qc_in, delimiter=",")))
    with open(join(report_dir, "qc_eval.tsv")) as qc_in:
        qc_pass = dict((row[0], row) for row in csv.reader(qc_in, delimiter="\t"))

    if args.mode == "asm":
        quast_report, tax_report = "quast_report.tsv", "blobtools_report.tsv"
        try:
            with open(join(report_dir, "eb_taxonomy_report.tsv")) as eb_in:
                tax_eb_report = dict((row[0], row[-1]) for row in csv.reader(eb_in, delimiter="\t"))
        except:
                tax_eb_report = dict() 
        try:
            with open(join(report_dir, "all_quast_taxonomy_report.tsv")) as asm_eb_in:
                asm_eb_report = dict((row[0], row[-6:]) for row in csv.reader(asm_eb_in, delimiter="\t"))
        except:
                asm_eb_report = dict()
    else:
        quast_report, tax_report = "quast_survey_report", "blobtools_survey_report.tsv"
        tax_eb_report, asm_eb_report = dict(), dict()

    with open(join(report_dir, quast_report)) as quast_in:
        quast_data = dict((row[0], row) for i, row in enumerate(csv.reader(quast_in, delimiter="\t")))
    with open(join(report_dir, tax_report)) as tax_in:
        tax_report = dict((row[0], row) for i, row in enumerate(csv.reader(tax_in, delimiter="\t")))

    writeDatasum(data, header, args.datasum_file.replace(".xlsx", ".bgsum.xlsx"), qc_pass=qc_pass, asm_report=quast_data, tax_report=tax_report, tax_eb_report=tax_eb_report, asm_eb_report=asm_eb_report)

    pass
