from Bio import SeqIO
hla_all = open("./data/hla_prot.fasta", 'r')
hla_class_1 = open("./output/hla_class_1.fasta", 'w+')
record_dict = SeqIO.to_dict(SeqIO.parse(hla_all, "fasta"))

for i in record_dict.keys():
    if " A*" in record_dict[i].description or " B*" in record_dict[i].description or " C*" in record_dict[i].description:
        hla_class_1.write(str('>' + record_dict[i].description + '\n' + record_dict[i].seq) + '\n')
        
hla_all.close()
hla_class_1.close()