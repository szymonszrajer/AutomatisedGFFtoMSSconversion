#!/usr/bin/python
# coding: UTF-8

import sys
from Bio import SeqIO
from BCBio import GFF
import re
import copy

#monkey patching a function from a GFF module in order to make it not sort the GFF file

def monkey_parse_in_parts(self, gff_files, base_dict=None, limit_info=None, target_lines=None):
        """Parse a region of a GFF file specified, returning info as generated.

        target_lines -- The number of lines in the file which should be used
        for each partial parse. This should be determined based on available
        memory.
        """
        for results in self.parse_simple(gff_files, limit_info, target_lines):
            if base_dict is None:
                cur_dict = dict()
            else:
                cur_dict = copy.deepcopy(base_dict)
            cur_dict = self._results_to_features(cur_dict, results)
            all_ids = list(cur_dict.keys())
            # all_ids.sort()               <---------------- the patched part
            for cur_id in all_ids:
                yield cur_dict[cur_id]

GFF.parse_in_parts = monkey_parse_in_parts


#-----------------------------------------------
#Part 1 Sorting


# input
args = sys.argv
gff_file = args[1] # gff input file
fasta_file = args[2] # fasta input file 
mol_type = args[3]  # type of molecule 
organism = args[4] # type of organism
strain = args[5] # organism strain
country = args[6] # country

"""
array = []

with fileinput.FileInput(gff_file) as in_file:
    for line in in_file: # for every line
        if line[0] != '#' :
            whole_line = line
            whole_line = whole_line.strip('\n')
            line = line.split('\t')
            rec_id = int(re.findall(r"\d+",line[0])[0])
            start = int(line[3])
            if  None != re.findall(r"Parent=",line[8]): parent = 1 
            else : parent = 0
            array.append([rec_id,start,parent,whole_line])


gff_tobesorted = pd.DataFrame(array, columns=("rec_id","start","parent","line"))
sorted_gff = gff_tobesorted.sort_values(by=['rec_id','start','parent'],ascending=[True,True,False])


sorted_gff.to_csv('gff_sorted',columns=['line'],index=False,header=None)



with fileinput.FileInput("gff_sorted", inplace=True) as in_file:
    for line in in_file:
        line = line.replace('"','').strip("\n")
        print(line)

gff_file = "gff_sorted"
"""

#-----------------------------------------------
#Part 2 GFF to MSS file conversion


# function handling transformation of CDS/exon records in gff file
def CDSmRNA(count, CDS_f, position, out_joint, out_joint_close):
                count += 1
                transl_table = CDS_f.qualifiers.get("transl_table", ["1"])
                if count==1: #First CDSs/exons in the corresponding mRNA
                    CDS_start = CDS_f.location.start + 1
                    CDS_end = CDS_f.location.end
                    position = position + str(CDS_start) + ".." + str(CDS_end)
                else: #The second and subsequent CDSs/exons in the corresponding mRNA
                    CDS_start = CDS_f.location.start +1
                    CDS_end = CDS_f.location.end
                    position = position + ","+ str(CDS_start) + ".." + str(CDS_end)
                    out_joint = "join(" # more than one CDS/exon
                    out_joint_close = ")"
                return count, transl_table, position, out_joint, out_joint_close
                
# function printing resulting CDS 
def printCDS(out_strand, out_joint, position, out_joint_close, out_strand_close, locus_tag_ID, mRNA_ID, product_name, transl_table):
                
                noteprint = 0
                case = re.search(r"(Similar to [Uu]ncharacterized)", product_name[0])
                if case:
                    product = "uncharacterized protein" 
                else:
                    case = re.search(r"(Similar to)", product_name[0])
                    if case:
                        case2 = re.search(r"(Similar to )([\'\:\(\)\[\]\/\\A-Za-z0-9\._-]+)(:)(\s)(.+)",product_name[0])
                        if case2:
                            product = case2.group(5)
                            noteprint = 1
                            note =  case2.group(1) + case2.group(2) + " NCBI"

                        else:
                            case2 = re.search(r"(Similar to )(.+)",product_name[0])
                            product = case2.group(2)
                        
                    else:
                        product = product_name[0]


                print(f"\tCDS\t"+ out_strand + out_joint + position + out_joint_close + out_strand_close + "\tcodon_start\t1" )
                print("\t\t\t" + "locus_tag\t" + mRNA_ID[0])
                print("\t\t\t" + "product\t" + product)
                if (noteprint): print("\t\t\t" + "note\t" + note)
                print("\t\t\t" + "transl_table\t" + transl_table[0])
                


# function printing resulting mRNA
def printmRNA(out_strand, out_joint, position, out_joint_close, out_strand_close, locus_tag_ID, mRNA_ID, product_name, transl_table):
                print(f"\tmRNA\t"+ out_strand + out_joint + position + out_joint_close + out_strand_close + "\tlocus_tag\t" + mRNA_ID[0])     




# fetching information about assembly gaps and scaffolds from the files
fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
assembly_gaps = {} # scaffold_name:gap list dictionary
lengths = {} # dictionary of lenghts of Scaffolds - needed for a header of each Scaffold
scaffolds = {} # dictionary with names of scaffolds based on order
precedence = {} # dictionary of scaffold precedence based on the name

order = 0 # iterating through order 
for fasta in fasta_sequences:
    name, length = fasta.id, str(len(fasta.seq))  
    lengths[name] = length 
    scaffolds[order] = name 
    precedence[name] = order 
    gaps = ([[m.start(0), m.end(0)] for m in re.finditer(r"N{10,}", str(fasta.seq))]) # annotating assembly gaps
    assembly_gaps[name] = gaps
    order += 1 


pre_contig = "" # variable used to check which Scaffold the program is on in order to add headers 
i = 0 # iterating through scaffolds
cap = len(precedence) # last scaffold precedence




in_handle = open(gff_file)
for rec in GFF.parse(in_handle): # for every scaffold   
    now_contig = rec.id # contig name from a line in a gff file
    now_conting_begin = 1 
    now_contig_end = lengths[now_contig] # contig end from a lengths dictionary
    # header of Scaffolds using the lenghts dictionary 
    if pre_contig != now_contig: # new header on new contig encounter
        # going through scaffolds until the last scaffold with features other than assembly gap is annotated
        while (precedence[now_contig]>=precedence[scaffolds[i]]):
                scaffold = scaffolds[i]
                print(scaffold + "\t" + "source" + "\t" + str(1) + ".." + lengths[scaffold] + "\t" + "ff_definition" + "\t" + "@@[organism]@@ DNA, @@[submitter_seqid]@@" + "\n" + \
                "\t" + "\t" + "\t" + "mol_type" + "\t" + mol_type + "\n" + \
                "\t" + "\t" + "\t" + "organism" + "\t" + 	organism + "\n" + \
                "\t" + "\t" + "\t" + "strain" + "\t" + strain + "\n" + \
                "\t" + "\t" + "\t" + "submitter_seqid" + "\t" + "@@[entry]@@" + "\n" + \
                "\t" + "\t" + "\t" + "country" + "\t" + country)       

                for list in assembly_gaps[scaffold]:
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
                out_strand = "complement("
                out_strand_close = ")"
            else :
                out_strand=""
                out_strand_close=""

            count = 0 # set count of CDS to 0 when entering new RNA
            position="" # output
            count2 = 0 # set count of mRNA to 0 when entering new RNA
            position2="" # output

            # more than one CDS/exon within a feature calls for join in MSS format
            out_joint = ""
            out_joint_close=""
            out_joint2 = ""
            out_joint_close2=""

            for CDS_f in mRNA_f.sub_features: # CDS/mRNA within a subfeature

                if str(CDS_f.type) == "CDS": # handling the CDS records
                    count, transl_table, position, out_joint, out_joint_close = CDSmRNA(count, CDS_f, position, out_joint, out_joint_close)      

                elif str(CDS_f.type) == "exon": # handling the exon records forming mRNA
                    count2, transl_table, position2, out_joint2, out_joint_close2 = CDSmRNA(count2, CDS_f, position2, out_joint2, out_joint_close2)

                else :
                    transl_table = CDS_f.qualifiers.get("transl_table", ["1"])

            #printing CDS and mRNA 
            if count != 0 and count2 != 0: 
            
                printCDS(out_strand, out_joint, position, out_joint_close, out_strand_close, locus_tag_ID, mRNA_ID, product_name, transl_table)
                printmRNA(out_strand, out_joint2, position2, out_joint_close2, out_strand_close, locus_tag_ID, mRNA_ID, product_name, transl_table)


# iterating through the rest of the scaffolds and assembly gaps without features

for j in range (i, cap):   
        scaffold = scaffolds[j]
        print(scaffold + "\t" + "source" + "\t" + str(1) + ".." + lengths[scaffold] + "\t" + "ff_definition" + "\t" + "@@[organism]@@ DNA, @@[submitter_seqid]@@" + "\n" + \
                "\t" + "\t" + "\t" + "mol_type" + "\t" + mol_type + "\n" + \
                "\t" + "\t" + "\t" + "organism" + "\t" + 	organism + "\n" + \
                "\t" + "\t" + "\t" + "strain" + "\t" + strain + "\n" + \
                "\t" + "\t" + "\t" + "submitter_seqid" + "\t" + "@@[entry]@@" + "\n" + \
                "\t" + "\t" + "\t" + "country" + "\t" + country)        

        for list in assembly_gaps[scaffold]:
            print("\t" + "assembly_gap" + "\t" + str(list[0]+1) + ".." + str(list[1]) + "\t" + "estimated_length" + "\t" + "known" )
            print("\t" + "\t" + "\t" + "gap_type" + "\t" + "within scaffold" )
            print("\t" + "\t" + "\t" + "linkage_evidence" + "\t" + "paired-ends" )   



in_handle.close()
 



