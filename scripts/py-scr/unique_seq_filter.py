# running this script separates only unique sequences from hla_class_1.fasta

from Bio import SeqIO
hla_class_1 = open("./output/hla_class_1.fasta", 'r')
unique_seq = open("./output/unique_seq.fasta", 'w+') 
record_dict = SeqIO.to_dict(SeqIO.parse(hla_class_1, "fasta"))
dict_inverse = {}
for i, j in record_dict.items():
    dict_inverse[record_dict[i].description] = str(record_dict[i].seq)
dict_inverse = {j:i for i, j in dict_inverse.items()}

for seq, desc in dict_inverse.items():
    if "N" not in desc: 
        unique_seq.write(">" + desc + "\n" + seq + "\n")

hla_class_1.close()
unique_seq.close()