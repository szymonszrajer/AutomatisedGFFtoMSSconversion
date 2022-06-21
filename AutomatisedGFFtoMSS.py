#!/usr/bin/python
# coding: UTF-8

import sys
from Bio import SeqIO
from Bio import Seq
from BCBio import GFF
import fileinput
import re
import csv
import pandas as pd


# function handling transformation of CDS/exon records in gff file
def CDSmRNA(COUNT, CDS_f, POSITION, out_JOINT, out_JOINT_CLOSE):
                COUNT += 1
                transl_table = CDS_f.qualifiers.get("transl_table", ["1"])
                if COUNT==1: #First CDSs/exons in the corresponding mRNA
                    CDS_START = CDS_f.location.start + 1
                    CDS_END = CDS_f.location.end
                    POSITION = POSITION + str(CDS_START) + ".." + str(CDS_END)
                else: #The second and subsequent CDSs/exons in the corresponding mRNA
                    CDS_START = CDS_f.location.start +1
                    CDS_END = CDS_f.location.end
                    POSITION = POSITION + ","+ str(CDS_START) + ".." + str(CDS_END)
                    out_JOINT = "join(" # more than one CDS/exon
                    out_JOINT_CLOSE = ")"
                return COUNT, transl_table, POSITION, out_JOINT, out_JOINT_CLOSE
                
# function printing resulting CDS 
def printCDS(out_STRAND, out_JOINT, POSITION, out_JOINT_CLOSE, out_STRAND_CLOSE, locus_tag_ID, mRNA_ID, product_name, transl_table):
                print(f"\tCDS\t"+ out_STRAND + out_JOINT + POSITION + out_JOINT_CLOSE + out_STRAND_CLOSE + "\tcodon_start\t1" )
                print("\t\t\t" + "locus_tag\t" + mRNA_ID[0])
                #print("\t\t\t" + "note\t" + mRNA_ID[0])
                print("\t\t\t" + "product\t" + product_name[0])
                print("\t\t\t" + "transl_table\t" + transl_table[0])

# function printing resulting mRNA
def printmRNA(out_STRAND, out_JOINT, POSITION, out_JOINT_CLOSE, out_STRAND_CLOSE, locus_tag_ID, mRNA_ID, product_name, transl_table):
                print(f"\tmRNA\t"+ out_STRAND + out_JOINT + POSITION + out_JOINT_CLOSE + out_STRAND_CLOSE + "\tlocus_tag\t" + mRNA_ID[0])            #+ "\tcodon_start\t1")
                #print("\t\t\t" + "locus_tag\t" + locus_tag_ID[0])
                #print("\t\t\t" + "note\t" + mRNA_ID[0])
                #print("\t\t\t" + "product\t" + product_name[0])
                #print("\t\t\t" + "transl_table\t" + transl_table[0])

args = sys.argv
gff_file = args[1] # gff input file
fasta_file = args[2] # fasta input file 



a = {}
fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
lengths = {} # dictionary of lenghts of Scaffolds - needed for a header of each Scaffold
scaffolds = {}
precedence = {}
order = 0 
for fasta in fasta_sequences:
    name, length = fasta.id, str(len(fasta.seq))  
    lengths[name] = length 
    scaffolds[order] = name 
    precedence[name] = order 
    gaps = ([[m.start(0), m.end(0)] for m in re.finditer(r"N{10,}", str(fasta.seq))]) 
    a[name] = gaps
    order += 1 


PreContig = "" # variable used to check which Scaffold the program is on in order to add headers 
Contig_Count = 0     
i = 0
cap = len(precedence)

filename = "scafforder"
fileW = open(filename, "w")
  

"""
for n in range (0,cap):
    limit_info = scaffolds[n]
    fileW.write(str(n) +  "\n")
    in_handle = open(gff_file)
"""    
in_handle = open(gff_file)
for rec in GFF.parse(in_handle): # for every scaffold   
    fileW.write(rec.id + "\n")
    NowContig = rec.id # contig name from a line in a gff file
    NowContingBegin = 1 
    NowContigEnd = lengths[NowContig] # contig end from a lengths dictionary
    # header of Scaffolds using the lenghts dictionary 
    if PreContig != NowContig: # new header on new contig encounter
        while (precedence[NowContig]>=precedence[scaffolds[i]]):
                scaffold = scaffolds[i]
                print(scaffold + "\t" + "source" + "\t" + str(1) + ".." + lengths[scaffold] + "\t" + "ff_definition" + "\t" + "@@[organism]@@ DNA, @@[submitter_seqid]@@" + "\n" + \
                "\t" + "\t" + "\t" + "mol_type" + "\t" + "genomic DNA" + "\n" + \
                "\t" + "\t" + "\t" + "organism" + "\t" + 	"Gryllus bimaculatus" + "\n" + \
                "\t" + "\t" + "\t" + "strain" + "\t" + "white eyes" + "\n" + \
                "\t" + "\t" + "\t" + "submitter_seqid" + "\t" + "@@[entry]@@" + "\n" + \
                "\t" + "\t" + "\t" + "country" + "\t" + "Japan:Tokushima")       

                for list in a[scaffold]:
                    print("\t" + "assembly_gap" + "\t" + str(list[0]+1) + ".." + str(list[1]) + "\t" + "estimated_length" + "\t" + "known" )
                    print("\t" + "\t" + "\t" + "gap_type" + "\t" + "within scaffold" )
                    print("\t" + "\t" + "\t" + "linkage_evidence" + "\t" + "paired-ends" )   

                if i < cap - 1 : i += 1                                 


    PreContig = rec.id # remembering current contig as an old contig in next iteration

    for gene_f in rec.features: # for every gene

        for mRNA_f in gene_f.sub_features: # for every CDS-mRNA subfeature
            # gene informations
            mRNA_ID = mRNA_f.qualifiers["ID"]
            locus_tag_ID = mRNA_f.qualifiers.get("Alias", ["Unknown_gene"])
            product_name = mRNA_f.qualifiers.get("Note", ["Unknown_product"])

            strand = mRNA_f.strand # orientation of the gene
            if strand == -1:
                out_STRAND = "complement("
                out_STRAND_CLOSE = ")"
            else :
                out_STRAND=""
                out_STRAND_CLOSE=""

            COUNT = 0 # set count of CDS to 0 when entering new RNA
            POSITION="" # output
            COUNT2 = 0 # set count of mRNA to 0 when entering new RNA
            POSITION2="" # output

            # more than one CDS/exon within a feature calls for join in MSS format
            out_JOINT = ""
            out_JOINT_CLOSE=""
            out_JOINT2 = ""
            out_JOINT_CLOSE2=""

            for CDS_f in mRNA_f.sub_features: # CDS/mRNA within a subfeature

                if str(CDS_f.type) == "CDS": # handling the CDS records
                    COUNT, transl_table, POSITION, out_JOINT, out_JOINT_CLOSE = CDSmRNA(COUNT, CDS_f, POSITION, out_JOINT, out_JOINT_CLOSE)      

                elif str(CDS_f.type) == "exon": # handling the exon records forming mRNA
                    COUNT2, transl_table, POSITION2, out_JOINT2, out_JOINT_CLOSE2 = CDSmRNA(COUNT2, CDS_f, POSITION2, out_JOINT2, out_JOINT_CLOSE2)

                else :
                    transl_table = CDS_f.qualifiers.get("transl_table", ["1"])

            #printing CDS and mRNA 
            if COUNT != 0 and COUNT2 != 0: 
            
                printCDS(out_STRAND, out_JOINT, POSITION, out_JOINT_CLOSE, out_STRAND_CLOSE, locus_tag_ID, mRNA_ID, product_name, transl_table)
                printmRNA(out_STRAND, out_JOINT2, POSITION2, out_JOINT_CLOSE2, out_STRAND_CLOSE, locus_tag_ID, mRNA_ID, product_name, transl_table)


for j in range (i, cap):   
        scaffold = scaffolds[j]
        print(scaffold + "\t" + "source" + "\t" + str(1) + ".." + lengths[scaffold] + "\t" + "ff_definition" + "\t" + "@@[organism]@@ DNA, @@[submitter_seqid]@@" + "\n" + \
        "\t" + "\t" + "\t" + "mol_type" + "\t" + "genomic DNA" + "\n" + \
        "\t" + "\t" + "\t" + "organism" + "\t" + 	"Gryllus bimaculatus" + "\n" + \
        "\t" + "\t" + "\t" + "strain" + "\t" + "white eyes" + "\n" + \
        "\t" + "\t" + "\t" + "submitter_seqid" + "\t" + "@@[entry]@@" + "\n" + \
        "\t" + "\t" + "\t" + "country" + "\t" + "Japan:Tokushima")       

        for list in a[scaffold]:
            print("\t" + "assembly_gap" + "\t" + str(list[0]+1) + ".." + str(list[1]) + "\t" + "estimated_length" + "\t" + "known" )
            print("\t" + "\t" + "\t" + "gap_type" + "\t" + "within scaffold" )
            print("\t" + "\t" + "\t" + "linkage_evidence" + "\t" + "paired-ends" )   



in_handle.close()
 

# fixes the "First codon is not a start codon error"
def first(location):
    rewrite = location.split('(')
    if (rewrite[0] == "complement"): # complement(....) and complement(join(....)) case , fix: (xxx..xxx,xxx..>xxx) (including complement/complemntjoin)
        adding_sign = location.rsplit("..",1) 
        location = adding_sign[0] + "..>" + adding_sign[1]

    else : # (if rewrite[0] != complement)
        if (rewrite[0] == "join"):  # join(....) case , fix: (<xxx..xxx,xxx..xxx) (including join)
            location = "join(<" + rewrite[1]
        else :
            location = "<" + location # (....) case , fix: (<xxx..xxx,xxx..xxx)
    return location

# fixes the "Final codon is not a stop codon error"
def final(location):
    rewrite = location.split('(')
    if (rewrite[0] == "complement"): 
        if (rewrite[1] == "join"): # complement(join(....)) case , fix: (<xxx..xxx,xxx..xxx) (including complementjoin)
            location = "complement(join(<" + rewrite[2]
        else : # complement(....) case , fix: (<xxx..xxx,xxx..xxx) (including complement)
            location = "complement(<" + rewrite[1]

    else : # (if rewrite[0] != complement) # join(....) and (....) case , fix: (xxx..xxx,xxx..>xxx) (including join or not)
        adding_sign = location.rsplit("..",1) 
        location = adding_sign[0] + "..>" + adding_sign[1]

    return location


filename = "errors_dup_bef20"    # errors file
fileR = open(filename, "r")

filename = "what_to_fix"    # generated file with information how to solve things
fileW = open(filename, "w")




for line in fileR:
  if line[0] == '>': # new location of error
    location = line.split()[1]
    oldLocation = location
    while line[0:2] != "//": # until errors for this location dont end
        line = fileR.readline()
        if (line.split()[0] == "TC0017:ER2:"): # first codon error
            location = first(location)
        if (line.split()[0] == "TC0018:ER2:"): # final codon error
            location = final(location)
    uninterrupted = oldLocation + ";" + location + "\n" # write to fix file in "bad_location;location_fix" manner
    fileW.write(uninterrupted)


fileR.close()
fileW.close()

        


with fileinput.FileInput("mss_duplicates_20", inplace=True) as file:
    x = 0 # position in replacement file
    for line in file: # for line in mss file
            splitted = re.split(r'\t+', line) # splitted line by tabs
            flag = 0 #flag to check if replacement in the line is being made
            if splitted[1] == "CDS": # if given line is CDS of mRNA header
                with open("what_to_fix") as f: # with open file with things that need to be fixed
                    for linee in f.readlines()[x:]: # opening the file from a new position as to not overwrite the fixed errors
                        replacement = linee.split(";") # version giving an error and the new version are separated my a ';' in a file
                        r1 = replacement[0] # old version
                        r2 = replacement[1] # new version
                        if r1 == splitted[2]: # if location is approptiate
                            x += 1 # new file opening position
                            flag = 1 # replacement needs to be done                           ###
                            savedR1 = r1.rstrip('\n') # stripping the line of the newline     ###
                            savedR2 = r2.rstrip('\n')                                         ###
                            break;
            if flag == 1 : print(line.replace(savedR1,savedR2), end='') # replacement to be done ###
            else :  print(line, end = '') # no replacement


array = []
gff_file = "gff_tosort.gff"
in_handle = open(gff_file)

for line in in_handle : # for every line
    if line[0] != '#' :
        whole_line = line
        line = line.split('\t')
        rec_id = int(re.findall(r"\d+",line[0])[0])
        start = int(line[3])
        if  None != re.findall(r"Parent=",line[8]): parent = 1 
        else : parent = 0
        array.append([rec_id,start,parent,whole_line])


gff_tobesorted = pd.DataFrame(array, columns=("rec_id","start","parent","line"))
sorted_gff = gff_tobesorted.sort_values(by=['rec_id','start','parent'],ascending=[True,True,False])
"""
print(sorted_gff)
sorted_a = []
sorted_a = sorted_gff.loc[:,"line"]
"""

sorted_gff.to_csv('gff_sorted',columns=['line'],index=False,header=None)


"""
with open("gff_sorted", "w") as out_handle:
    out_handle.write("##gff-version 3" + "\n")
    for i in range (len(sorted_a)):
        print(sorted_gff.loc[i,"line"])
        out_handle.write(sorted_gff.loc[i,"line"])
"""




filename = "mss_duplicates_20_re"
fileW = open(filename, "w")
#res = []

with fileinput.FileInput("mss_duplicates_20") as file:
    for line in file:
        line = re.sub(r"\(\w+\s\w+.+\)","",line)
        line = re.sub(r"\\","",line)
        match  = re.search(r"(Scaffold)([0-9]+)",line)
        if match:
            replacement = "BOPP0" + str(int(match.group(2)) + 1000000)
            line = re.sub(match.group(0),replacement,line)
        fileW.write(line)



fileW.close()




