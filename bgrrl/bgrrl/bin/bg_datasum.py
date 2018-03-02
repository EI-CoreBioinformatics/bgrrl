import sys
import argparse

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
        data[sheet.row_values(i)[0]] = dict(zip(header, sheet.row_values(i)))
    return data, header


def writeDatasum(data, header, dfile):
    """
    import xlwt
    book = xlwt.Workbook()
    sheet1 = book.add_sheet("Sheet1")
    # sheet1 = book.get_sheet(0)
    row = sheet1.row(0)
    for j, col in enumerate(header):
        row.write(j, col)

    for i, k in enumerate(sorted(data), 1):
        row = sheet1.row(i)
        for j, col in enumerate(header):
            row.write(j, data[k][col])

    book.save(dfile)
    """
    def write_row(sheet, row_index, data, header):
        for j, col in enumerate(header):
            sheet.write(row_index, j, data[col])


    import xlsxwriter
    book = xlsxwriter.Workbook(dfile)
    sheet_ip = book.add_worksheet("In Pool")
    sheet_nip = book.add_worksheet("Not in Pool")
    for j, col in enumerate(header):
        sheet_ip.write(0, j, col)
        sheet_nip.write(0, j, col)
    row_index_ip, row_index_nip = 1, 1
    for i, k in enumerate(sorted(data), 1):
        if data[k][header[-1]]:
            write_row(sheet_ip, row_index_ip, data[k], header)
            row_index_ip += 1
        else:
            write_row(sheet_nip, row_index_nip, data[k], header)
            row_index_nip += 1

    book.close()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("indir", type=str, help="Path to input directory.")
    ap.add_argument("datasum_file", type=str, help="Path to datasummary file. (xlsx)")
    args = ap.parse_args()

    data, header = readDatasum(args.datasum_file)
    for k in data:
        print(data[k])

    writeDatasum(data, header, args.datasum_file.replace(".xlsx", ".bgsum.xlsx"))

    pass