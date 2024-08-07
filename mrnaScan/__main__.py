import os, sys, re, getopt, functools, pysam, tables
import pandas as pd
import numpy as np
from collections import defaultdict
from multiprocess import Pool
from scipy import sparse

from mrnaScan import __version__

class ArgumentList:
    file_bam = ""
    file_ref = ""
    file_name = ""
    file_bc = ""
    chr_filter = ""
    chr_list = ""
    cn_len = 10
    num_thread = 1
    def __init__(self):
        self.file_bam = ""
        self.file_ref = ""
        self.file_name = ""
        self.file_bc = ""
        self.chr_filter = ""
        self.chr_list = ""
        self.cn_len = 10
        self.num_thread = 1

def chr_cmp(a, b):
    sa = str(a)
    sb = str(b)
    la = len(sa)
    lb = len(sb)
    lm = min(la, lb)
    for i in range(0, lm):
        if sa[i] != sb[i]:
            oa = ord(sa[i]) if sa[i] != "M" and sa[i] != "m" else 0x7A
            ob = ord(sb[i]) if sb[i] != "M" and sb[i] != "m" else 0x7A
            if oa < 0x3A and oa > 0x2F and ob < 0x3A and ob > 0x2F and la != lb:
                return la - lb
            cd = oa - ob
            return cd
    return la - lb

def getTrans(name, exons):
    exons = sorted(exons, key=lambda x: str(x[0])+"_"+str(x[1]))
    res = ""
    length = 0
    start, end, bps = exons[0]
    for region in exons[1:]:
        if (region[0] >= end) or (region[0] - start > 8 and end - region[0] < min(region[1] - region[0], end - start) * 2 / 3):
            length += end- start
            res += "_" + str(start) + "-" + str(end)
            start, end, bps = region
        elif region[2] > bps:
            start, end, bps = region
    length += end- start
    res += "_" + str(start) + "-" + str(end)
    return name+"_"+str(length)+res

def chr_handle_exon(chr, ref_raw, cells, samfile):
    ref = ref_raw[ref_raw["seq_id"] == chr]
    ref.index = range(1, ref.shape[0] + 1)
    samfile = pysam.AlignmentFile(samfile)
    trans_info = defaultdict(lambda: defaultdict(list))
    feature_count = 0
    for index, row in ref.iterrows():
        if index % 10000 == 0:
            print(chr+": "+str(index)+" terms were processed...")
        chr, start, end, gene, length = list(row)
        if length < 4:
            continue
        for reads in samfile.fetch(chr, start, end):
            cb = reads.get_tag("CB")
            if not cb in cells:
                continue
            ref_pos, cigar_list = reads.reference_start+1, reads.cigar
            if any(map(lambda x : x[0] > 4, cigar_list)):
                continue
            bps = 0
            block_count = 0
            block_end = 0 if reads.flag & 16 == 0 else [x[0] for x in cigar_list].count(0) - 1
            end_test = False
            for cigar in cigar_list: # S=4;M=0;D=2;I=1:N=3
                if ref_pos > end:
                    break
                if cigar[0] == 0:
                    if bps < length and ref_pos + cigar[1] > start:
                        if block_count == block_end:
                            end_test = True
                        bps += cigar[1]
                        if ref_pos < start:
                            bps -= start - ref_pos
                    ref_pos += cigar[1]
                elif cigar[0] == 2 or cigar[0] == 3:
                    ref_pos += cigar[1]
            if end_test or bps / length > 0.8:
                trans_info[reads.qname+"_"+cb][gene].append([start, end, bps])
        feature_count += 1
    res = defaultdict(lambda: defaultdict(int))
    for qname, info in trans_info.items():
        id = ["", 0, 0]
        for k, v in info.items():
            l = sum([x[2] for x in v])
            if len(v) > id[1] or (len(v) == id[1] and l > id[2]):
                id = [k, len(v), l]
        res[getTrans(id[0], info[id[0]])][cells[qname.split("_")[-1]]] += 1
    samfile.close()
    return res

def chr_handle_utr(chr, ref_raw, cells, samfile, pad):
    ref = ref_raw[ref_raw["seq_id"] == chr]
    ref.index = range(1, ref.shape[0] + 1)
    samfile = pysam.AlignmentFile("tagged.bam")
    utr_dict = defaultdict(lambda: defaultdict(list))
    for index, row in ref.iterrows():
        if index % 10000 == 0:
            print(chr+": "+str(index)+" terms were processed...")
        chr, start, end, strand, gene = list(row)
        for reads in samfile.fetch(chr, start, end):
            cb = reads.get_tag("CB")
            if not cb in cells:
                continue
            ref_pos, cigar_list = reads.reference_start+1, reads.cigar
            if any(map(lambda x : x[0] > 4, cigar_list)):
                continue
            pos = 0
            for cigar in cigar_list: # S=4;M=0;D=2;I=1:N=3
                if ref_pos > start:
                    break
                if cigar[0] == 0:
                    if ref_pos <= start and ref_pos + cigar[1] > start:
                        pos += start - ref_pos
                        utr_len = 0
                        if reads.flag & 16 == 0:
                            t = 0 if cigar_list[-1][0] != 4 else cigar_list[-1][1]
                            utr_len = reads.rlen - t - pos - 3 # reads.seq[pos+3:reads.rlen-t]
                        else:
                            t = 0 if cigar_list[0][0] != 4 else cigar_list[0][1]
                            utr_len = pos - t # reads.seq[t:pos]
                        if utr_len > 0:
                            utr_dict[gene][cb].append(utr_len)
                        break
                    ref_pos += cigar[1]
                    pos += cigar[1]
                elif cigar[0] == 2 or cigar[0] == 3:
                    ref_pos += cigar[1]
                else:
                    pos += cigar[1]
    samfile.close()
    return utr_dict

def bedScan(args):
    (pathname, extension) = os.path.splitext(args.file_bam)
    (filepath, filename) = os.path.split(pathname)
    if args.file_name != "":
        filename = args.file_name
    pathname = os.path.join(filepath, filename)
    if not os.path.isfile(args.file_bam+".bai"):
        print("Creating the index file of bam...")
        if os.system("which samtools") == 0:
            os.system("samtools index "+args.file_bam)
        else:
            print("Failed, please check your environment variables...")
            return
    print("Processing the reference...")
    ref_raw = pd.read_table(args.file_ref, comment="#", header=None, dtype={0:str})
    ref_raw.columns = ["seq_id", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    type_list = list(set(ref_raw["type"]))
    if any(map(lambda x: x not in type_list, ["exon", "stop_codon"])):
        print("There is no suitable term in the gtf for mRNA scaning...")
        return
    chr_list = list(set(ref_raw["seq_id"]))
    if args.chr_list != "":
        chr_list = list(set(chr_list).intersection(set(args.chr_list.split(","))))
    if args.cn_len > 0:
        chr_list = [x for x in chr_list if len(x) < min(args.cn_len, min(list(map(lambda x: len(x), chr_list)))*4)]
    if len(chr_list) == 0:
        print("There is no chromosome would be calculated...")
        return
    fs = pysam.AlignmentFile(args.file_bam, "rb")
    chr_detect = []
    fs_header = fs.header.to_dict()
    if "SQ" in fs_header.keys():
        for term in fs_header["SQ"]:
            if "SN" in term.keys():
                chr_detect.append(str(term["SN"]))
    if len(chr_detect) > 0:
        chr_list = list(set(chr_list).intersection(set(chr_detect)))
    if len(chr_list) == 0:
        print("There is no chromosome would be calculated...")
        return
    chr_list = sorted(list(set(chr_list).difference(set(args.chr_filter.split(",")))), key=functools.cmp_to_key(chr_cmp))
    chr_list = [x for x in chr_list if re.match(r".*[Mm]+,*", x) == None]
    ref_exon = ref_raw.copy(deep=True)
    ref_exon = ref_exon[(ref_exon["type"] == "exon") & (ref_exon["seq_id"].isin(chr_list))]
    ref_exon["attributes"] = list(map(lambda x: re.findall(r'gene_id "(.*?)";', x)[0], ref_exon["attributes"]))
    ref_exon = ref_exon[["seq_id", "start", "end", "attributes"]].sort_values(by=["seq_id", "start", "end", "attributes"])
    ref_exon.drop_duplicates(subset=["seq_id", "start", "end", "attributes"], keep="first", inplace=True)
    ref_exon["len"] = ref_exon["end"] - ref_exon["start"]
    ref_exon["start"] = ref_exon["start"] - 1
    ref_utr = ref_raw.copy(deep=True)
    ref_utr = ref_utr[(ref_utr["type"] == "stop_codon") & (ref_utr["seq_id"].isin(chr_list))]
    ref_utr["attributes"] = list(map(lambda x: re.findall(r'gene_id "(.*?)";', x)[0], ref_utr["attributes"]))
    ref_utr = ref_utr[["seq_id", "start", "end", "strand", "attributes"]].sort_values(by=["seq_id", "start", "end", "attributes"])
    ref_utr.drop_duplicates(subset=["seq_id", "start", "attributes"], keep="last", inplace=True)
    cells = {}
    cell_count = 0
    if os.path.isfile(args.file_bc):
        bcs = pd.read_csv(args.file_bc, header=0, index_col=0)
        for term in bcs.index:
            cells[term[:-2]] = cell_count
            cell_count += 1
    else:
        cells["NNNNNNNN"] = 0
        cell_count += 1
    print("Scaning exons and 3'utf of reads...")
    pool = Pool(args.num_thread)
    res_list = {"exon": [], "utr": []}
    for chr in chr_list:
        res_list["exon"].append(pool.apply_async(chr_handle_exon, (chr, ref_exon, cells, args.file_bam,)))
        res_list["utr"].append(pool.apply_async(chr_handle_utr, (chr, ref_utr, cells, args.file_bam, 3,)))
    pool.close()
    pool.join()
    print("Saving the results...")
    feature_id = 0
    feature_counts = defaultdict(int)
    features = defaultdict(list)
    count_row = []
    count_col = []
    count_data = []
    for term in res_list["exon"]:
        res = term.get()
        for trans, counts in res.items():
            infos = trans.split("_")
            features["name"].append(infos[0] + "_" + str(feature_counts[infos[0]]))
            features["id"].append(feature_id)
            features["gene"].append(infos[0])
            features["pos"].append("_".join(infos[2:]))
            features["len"].append(infos[1])
            feature_counts[infos[0]] += 1
            for cid, c in counts.items():
                count_row.append(feature_id)
                count_col.append(cid)
                count_data.append(c)
            feature_id += 1
    count_matrix = sparse.csc_matrix((np.array(count_data), (np.array(count_row), np.array(count_col))), shape=(len(features["id"]), len(cells)))
    with tables.open_file(pathname+"_exon_info.h5", "w", title="") as fs:
        group = fs.create_group("/", "matrix", "Datasets of count")
        feature_group = fs.create_group(group, "features", "Genes and other features measured")
        for k, v in features.items():
            r = fs.create_array(feature_group, k, np.array(v))
        r = fs.create_array(group, "genes", np.array(features["name"]))
        r = fs.create_array(group, "barcodes", np.array(list(cells.keys())))
        r = fs.create_array(group, "data", count_matrix.data)
        r = fs.create_array(group, "indices", count_matrix.indices)
        r = fs.create_array(group, "indptr", count_matrix.indptr)
        r = fs.create_array(group, "shape", count_matrix.shape)
    with open(pathname+"_utr_info.tsv", "w") as fs:
        msg = "Gene\tCell\tUTRs\n" if len(cells) > 1 else "Gene\tUTRs\n"
        r = fs.write(msg)
        for term in res_list["utr"]:
            res = term.get()
            for gene, utrs in res.items():
                for k, v in utrs.items():
                    msg = gene+"\t"+k+"\t"+",".join([str(x) for x in v])+"\n" if len(cells) > 1 else gene+"\t"+",".join([str(x) for x in v])+"\n"
                    r = fs.write(msg)

def main():
    opts, args = getopt.getopt(sys.argv[1:], 
        "hi:r:s:b:f:c:n:t:", 
        ["help", "input=", "reference=", "save=", "barcode=", "filter=", "chr=", "nl=","nthread="])
    arguments = ArgumentList()
    help_flag = False
    help_info = "Usage:\nmrnaScan [options] -i <input.bam> -r <reference.gtf>\nArguments:\n"\
        +"-h, --help\t\tShow this help information\n"\
        +"-i, --input <file>\tA aligned & deduped BAM file\n"\
        +"-r, --reference <file>\tGTF genome annotation\n"\
        +"-s, --save <name>\tThe name of results (default: same as the bam)\n"\
        +"-b, --barcode <name>\tA file with barcodes as list (default: none, process as bulk)\n"\
        +"-f, --filter [aaa,bbb]\tThe list of chromosomes which should be filtered (default: none)\n"\
        +"-c, --chr [aaa,bbb]\tThe list of chromosomes would be used (default: all)\n"\
        +"-n, --nl [0-X]\tThe length limit of chromosome names (default: 10, 0 is no limit)\n"\
        +"-t, --nthread [1-X]\tThe number of processing thread (default: 1)\n"
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            help_flag = True
        elif opt in ("-i", "--input"):
            arguments.file_bam = arg
        elif opt in ("-r", "--reference"):
            arguments.file_ref = arg
        elif opt in ("-s", "--save"):
            arguments.file_name = arg
        elif opt in ("-b", "--barcode"):
            arguments.file_bc = arg
        elif opt in ("-f", "--filter"):
            arguments.chr_filter = arg
        elif opt in ("-c", "--chr"):
            arguments.chr_filter = arg
        elif opt in ("-n", "--nl"):
            arguments.cn_len = int(arg)
        elif opt in ("-t", "--nl"):
            if int(arg) > 0:
                arguments.num_thread = int(arg)
    print("mrnaScan - Version: "+__version__)
    if help_flag or arguments.file_bam == "" or arguments.file_ref == "":
        print(help_info)
        sys.exit()
    bedScan(arguments)

if __name__ == "__main__":
    main()
