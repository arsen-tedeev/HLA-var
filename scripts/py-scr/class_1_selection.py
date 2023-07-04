from Bio import SeqIO
hla_class_1 = open("./output/hla_class_1.fasta", 'w+')

record_dict = SeqIO.to_dict(SeqIO.parse("./data/hla_prot.fasta", "fasta"))
for i in record_dict.keys():
    if " A*" in record_dict[i].description or " B*" in record_dict[i].description or " C*" in record_dict[i].description:
        hla_class_1.write(record_dict[i].description + "\n")

hla_class_1.close()