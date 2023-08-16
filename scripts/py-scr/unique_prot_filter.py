# this script reads the file with unique sequences, and then writes the longest sequence for each protein to the new file
# (as, for example, each of the sequences A*01:01:99 (365 bp) and A*01:01:69 (181 bp) is unique, but these are "pieces" of the same protein)

from Bio import SeqIO
unique_seq = open("./output/unique_seq.fasta", 'r')
record_dict = SeqIO.to_dict(SeqIO.parse(unique_seq,"fasta"))
unique_prot = open("./output/unique_prot.fasta", 'w+')
my_dict = {}

for key in record_dict.keys():
    my_dict[record_dict[key].description.split(" ")[0]] = [record_dict[key].description.split(" ")[1], record_dict[key].description.split(" ")[2], record_dict[key].seq]

for g, h in list(my_dict.items()):
    my_dict[g] = [str(h[0].split(":")[0] + ":" + h[0].split(":")[1]), h[1], h[2]] # leaving only first two sets of digits 

fl_list = []
curr_prot = None

for i, j in my_dict.items():
    if curr_prot != None and curr_prot == j[0]: 
        if int(fl_list[-1][2]) < int(j[1]):
            fl_list.pop()
            fl_list.append((i, j[0], j[1], j[2]))
    else:
        if curr_prot == None:
            curr_prot = j[0]
            fl_list.append((i, j[0], j[1], j[2]))
        else:
            curr_prot = j[0]
            unique_prot.write(">" + fl_list[-1][0] + " " + fl_list[-1][1] + " " + fl_list[-1][2] + "\n" + str(fl_list[-1][3]) + "\n")
            fl_list.pop()
            fl_list.append((i, j[0], j[1], j[2]))

if fl_list != []:
    unique_prot.write(">" + fl_list[-1][0] + " " + fl_list[-1][1] + " " + fl_list[-1][2] + "\n" + str(fl_list[-1][3]) + "\n")
unique_seq.close()
unique_prot.close()