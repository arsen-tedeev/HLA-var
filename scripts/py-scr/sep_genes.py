from Bio import SeqIO
hla_a = open("./output/hla_a.fasta", 'a+')
hla_b = open("./output/hla_b.fasta", 'a+')
hla_c = open("./output/hla_c.fasta", 'a+')
freq_alleles = open("./output/freq_all.fasta", 'r')
names = open("./data/allele_names.txt", 'r')
unique_list = []
a_names = []
b_names = []
c_names = []

for i in names:
    unique_list.append(i[:-4].strip())

unique_list =  list(dict.fromkeys(unique_list)) #thus, alleles would differ by first set of digits

for name in unique_list:
    if str("A") in name:
        a_names.append(name)
    elif str("B") in name:
        b_names.append(name)
    else:
        c_names.append(name)

record_dict = SeqIO.to_dict(SeqIO.parse(freq_alleles,"fasta"))
for a_name in a_names:
    for z in record_dict.keys():
        if str(a_name) in record_dict[z].description:
            hla_a.write(str('>' + record_dict[z].description + '\n' + record_dict[z].seq) + '\n')
            break

for b_name in b_names:
    for k in record_dict.keys():
        if str(b_name) in record_dict[k].description:
            hla_b.write(str('>' + record_dict[k].description + '\n' + record_dict[k].seq) + '\n')
            break

for c_name in c_names:
    for p in record_dict.keys():
        if str(c_name) in record_dict[p].description:
            hla_c.write(str('>' + record_dict[p].description + '\n' + record_dict[p].seq) + '\n')
            break
        
hla_a.close()
hla_b.close()
hla_c.close()
freq_alleles.close()
names.close()
