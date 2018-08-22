from Bio import SeqIO
from Bio.Seq import Seq
import csv

def parse_var_file(varscan_file):
    """
    parses varscan file
    returns list of tuples with snp details
    """
    var_list = []
    with open(varscan_file, 'r') as f:
        for line in f:
            if line[:5] == 'Chrom':
                pass
            else:
                snp_line = line.split('\t', 5)
                snp_pos = int(snp_line[1])
                snp_ref = snp_line[2]
                snp_var = snp_line[3]
                freq_var = snp_line[4].split(':')[4]
                result = (snp_pos, snp_ref, snp_var, freq_var)
                var_list.append(result)
    f.close()
    return var_list

def get_genomic_seq(gbk):
    seq_record = SeqIO.read(gbk, "genbank")
    seq_obj = seq_record.seq
    return str(seq_obj)

def parse_gbk_cds(gbk):
    cds = []
    with open(gbk, 'r') as f:
        start = False
        for line in f:
            line = line.strip().split()

            if len(line) > 0 and line[0] == "CDS":
                results = []
                start = True
                if line[1][-1] != ")":
                    if "join" in line[1]:
                        for join in line[1].split("join(")[1].split(","):
                            if join != "" and join[0] != "<" and join[-1] != ">":
                                results.append(join)
                    else:
                        cds.append(line[1].strip())
                        start=False
                else:
                    start = False
                    if "<" not in line[1] and ">" not in line[1]:
                        cds.append(line[1].split("join(")[1][:-1].split(","))
            elif start:
                if line[0][-1] == ")":
                    start = False
                    for join in line[0][:-1].split(","):
                        if join != "" and join[0] != "<" and join[-1] != ">":
                            results.append(join)
                    cds.append(results)
                else:
                    for join in line[0].split(","):
                        if join != "" and join[0] != "<" and join[-1] != ">":
                            results.append(join)
    f.close()
    return cds


def calc_aa_changes(cdss, genomic_seq, var_list):
    all_changes = {}
    count = 1
    for cds in cdss:
        nucl_changes = get_changes(cds, genomic_seq, var_list)
        mod_var_list = []
        for i in range(len(var_list)):
            mod_var_list.append(var_list[i]+nucl_changes[i])
        all_changes[count] = mod_var_list
        count += 1
    return all_changes


def get_changes(cds, ori_genomic_seq, var_list):
    nucl_changes=len(var_list) * [("Intron","-","-")]
    var_list_pos = 0
    for a in var_list:
        ori_cds_seq = Seq("")
        mod_cds_seq = Seq("")
        mod_genomic_seq = ori_genomic_seq[:a[0]-1]+a[2]+ori_genomic_seq[a[0]:]
        pos = 0
        var_in_exon = False
        for exon in cds:
            exon_pos = exon.split("..")
            if int(exon_pos[0]) <= a[0] <= int(exon_pos[1]):
                var_in_exon = True                
            ori_cds_seq+=ori_genomic_seq[int(exon_pos[0])-1:int(exon_pos[1])]
            mod_cds_seq+=mod_genomic_seq[int(exon_pos[0])-1:int(exon_pos[1])]
        if var_in_exon:
            found = False
            ori_str = str(ori_cds_seq.translate())
            mod_str = str(mod_cds_seq.translate())
            while pos < len(ori_str):
                if ori_str[pos] != mod_str[pos]:
                    nucl_changes[var_list_pos] = (pos+1, ori_str[pos], mod_str[pos])
                    found = True
                pos+=1
            if not found:
                nucl_changes[var_list_pos] = ("No change","-","-")
        var_list_pos+=1
    return nucl_changes


def mod_gb(all_changes, gb_in_file, gb_out_file, csv_out_file):
    """
    parses gb_in_file
    inserts varscan SNP infor
    writes to new gb_out_file and simplified csv_out_file
    """
    insert_pos = False
    with open(gb_in_file, 'r') as f:
        with open(gb_out_file, 'w') as g:
            for line in f:
                if insert_pos:
                    for var in all_changes[1]:
                        g.write('     misc_feature     {0}..{0}\n'.format(str(var[0])))
                        g.write('                     /vntifkey="21"\n')
                        g.write('                     /label={0}-->{1}_{2}\n'.format(var[1], var[2], var[3]))
                        if var[-1] != "-":
                            g.write('     misc_feature     {0}..{0}\n'.format(str(var[0])))
                            g.write('                     /vntifkey="21"\n')
                            g.write('                     /label={0}-->{1}\n'.format(var[5], var[6]))                            
                    insert_pos = False
                    g.write(line)
                elif line[:8] == "FEATURES":
                    insert_pos = True
                    g.write(line)
                else:
                    g.write(line)
    g.close()
    f.close()
    with open(csv_out_file, 'w', newline='') as csvfile:
        spamwriter = csv.writer(csvfile)
        spamwriter.writerow(["CDS","Var nucl pos","Ref nucl", "Var nucl", "Var percent", "Var prot pos", "Ref aa", "Var aa"])
        for i in range(len(all_changes)):
            for j in all_changes[i+1]:
                row = list((i+1,)+(j))
                spamwriter.writerow(row)
    csvfile.close()
    
def annotate(in_gb, var_csv, out_gb, out_csv):
    var_list = parse_var_file(var_csv)
    genomic_seq = get_genomic_seq(in_gb)
    cdss = parse_gbk_cds(in_gb)
    all_changes = calc_aa_changes(cdss, genomic_seq, var_list)
    mod_gb(all_changes, in_gb, out_gb, out_csv)


IDS, = glob_wildcards("raw/{smp}_1.fq")
gbs, = glob_wildcards("gb/{gb}.gb")

rule alignments:
    input:
        expand("annot_csv/{smp}.{fa}.csv", smp=IDS, fa=gbs)

rule folders:
    output: "seq"
    shell: "mkdir {output}"

rule trimming:
    input:
        fwd = "raw/{smp}_1.fq",
        rev = "raw/{smp}_2.fq"
    output:
        fwd = "seq/{smp}_1_val_1.fq",
        rvs = "seq/{smp}_2_val_2.fq"
    shell:
        "trim_galore --paired --retain_unpaired --phred33 --length 36 -q 5 -o seq --stringency 1 -e 0.1 {input.fwd} {input.rev}"

rule converting_gb:
    input:
        ingb = "gb/{fa}.gb"

    output:
        outfa = "fa/{fa}.fa"
    run:
        SeqIO.convert(input.ingb, "gb", output.outfa, "fasta"),
        shell("bowtie2-build -f {output.outfa} {output.outfa}")

rule align:
    input:
        index = "fa/{fa}.fa",
        fwd = "seq/{smp}_1_val_1.fq",
        rvs = "seq/{smp}_2_val_2.fq",
    output:
        out_align = "align/{smp}.{fa}.bam"
    threads: 40
    shell:
        "bowtie2 -x {input.index} -p {threads} -1 {input.fwd} -2 {input.rvs} | samtools view -@ {threads} -bS - | samtools sort -@ {threads} - > {output.out_align}"

rule snp_id:
    input:
        fasta = "fa/{fa}.fa",
        alignment = "align/{smp}.{fa}.bam"
    output:
        csv_file = "snp/{smp}.{fa}.csv"
    shell:
        "samtools mpileup -f {input.fasta} {input.alignment}  | varscan mpileup2snp --min-var-freq 0.0075 --min-coverage 200 --min-reads2 30 > {output.csv_file}"

rule gb_annot:
    input:
        csv_file = "snp/{smp}.{fa}.csv",
        ingb = "gb/{fa}.gb"
    output:
        outgb = "annot_gb/{smp}.{fa}.gb",
        outcsv = "annot_csv/{smp}.{fa}.csv"
    run:
        annotate(input.ingb, input.csv_file, output.outgb, output.outcsv)
