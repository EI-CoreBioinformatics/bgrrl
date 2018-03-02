import os
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


def writeDatasum(data, header, dfile, qc_pass=set(), asm_report=dict()):
    def write_row(sheet, row_index, data, header):
        for j, col in enumerate(header):
            sheet.write(row_index, j, data[col])

    in_pool_cheader = header[-1]
    ip_header, nip_header = list(header[:-1]), list(header[:-1])
    ip_header.append("status")

    asm_header = asm_report.get("Assembly", list())
    # print("ASM", asm_header)
    # print(asm_report)
    # print(asm_report.keys())
    ip_header.extend(asm_header[1:])

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
            data[k]["status"] = "" if k in qc_pass else "failed Bioinf_QC"
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
    args = ap.parse_args()

    data, header = readDatasum(args.datasum_file)

    with open(os.path.join(args.indir, "reports", "samplesheets", "samplesheet.qc_pass.tsv")) as qc_in:
        qc_pass = set(row[0] for i, row in enumerate(csv.reader(qc_in, delimiter=",")))
    with open(os.path.join(args.indir, "reports", "quast_report.tsv")) as quast_in:
        quast_data = dict((row[0], row) for i, row in enumerate(csv.reader(quast_in, delimiter="\t")))


    # for k in data:
    #    print(data[k])

    writeDatasum(data, header, args.datasum_file.replace(".xlsx", ".bgsum.xlsx"), qc_pass=qc_pass, asm_report=quast_data)

    pass
