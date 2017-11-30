from Bio import SeqIO


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
    print(var_list)
    return var_list

def mod_gb(a,var_list, gb_in_file, gb_out_file, csv_out_file):
    """
    parses gb_in_file
    inserts carscan SNP infor
    writes to new gb_out_file and simplified csv_out_file
    """
    compiled_list = []
    print(a)
    insert_pos = False
    with open(gb_in_file, 'r') as f:
        with open(gb_out_file, 'w') as g:
            for line in f:
                if insert_pos:
                    print(var_list)
                    for var in var_list[0]:
                        print(var)
                        g.write('     misc_feature     {0}..{0}\n'.format(str(var[0])))
                        g.write('                     /vntifkey="21"\n')
                        g.write('                     /label={0}-->{1}_{2}\n'.format(var[1], var[2], var[3]))
                        compiled_list.append(
                            (gb_in_file, var[0], var[1], var[2], var[3]))
                    insert_pos = False
                    g.write(line)
                elif line[:8] == "FEATURES":
                    insert_pos = True
                    g.write(line)
                else:
                    g.write(line)
    g.close()
    f.close()
    with open(csv_out_file, 'w') as h:
        h.write("Gene, Position, Reference, Variant, Frequency\n")
        for var in compiled_list:
            h.write("{0},{1},{2},{3},{4}\n".format(var[0].split('/')[1].split('.')[0],
                                                   var[1],
                                                   var[2],
                                                   var[3],
                                                   var[4]))
    h.close()


IDS, = glob_wildcards("fastq/{smp}_1.fastq")
gbs, = glob_wildcards("gb/{gb}.gb")

rule alignments:
    input:
        expand("annot_csv/{smp}.{fa}.csv", smp=IDS, fa=gbs)

rule folders:
    output: "seq"
    shell: "mkdir {output}"
rule trimming:
    input:
        fwd = "fastq/{smp}_1.fastq",
        rev = "fastq/{smp}_2.fastq"
    output:
        fwd = "seq/{smp}_1_val_1.fq",
        rvs = "seq/{smp}_2_val_2.fq"
    threads: 12
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
    threads: 12
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
        var_list = parse_var_file(input.csv_file),
        a=8,
        compiled_list = mod_gb(a,var_list, input.ingb,
                               output.outgb, output.outcsv)
