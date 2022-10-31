from Bio import SeqIO
names = open("./data/allele_names.txt", "r")
selected_alleles = open("./output/freq_all.fasta", "a+")

record_dict = SeqIO.to_dict(SeqIO.parse("./data/hla_prot.fasta", "fasta"))
for name in names:
    for i in record_dict.keys():
        if str(" " + name.strip() + ":") in record_dict[i].description or str(" " + name.strip() + " ") in record_dict[i].description:
            selected_alleles.write(str(">" + record_dict[i].description + "\n" + record_dict[i].seq) + "\n")
            break
        
names.close()
selected_alleles.close()