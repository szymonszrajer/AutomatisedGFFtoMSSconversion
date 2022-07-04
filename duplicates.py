import re
import fileinput


filename = "slines"
fileW = open(filename, "w")

with fileinput.FileInput("duplicates.txt") as file:
    for line in file:
        match  = re.search(r"(JP0080:WAR:STX:ANN:Lines )(\[.+?\])( and )(\[.+?\])",line)
        if match:
            match = str(match.group(4))
            match = match.replace(u'\xa0',u'')
            match = match.replace(u'[',u'')
            match = match.replace(u']',u'')
            fileW.write(str(int(match) + 2) + "\n")

fileW.close()


with open("slines") as f:
    sl = f.read().splitlines()
sl = [int(s) for s in sl]

index = 1

filename = "GBI_IDS"
fileW = open(filename, "w")

with fileinput.FileInput("mss_fixing_start99_9_re") as file:
    for line in file:
        if index in sl:
            match  = re.search(r"GBI_[0-9]+-R[A-Z]",line)
            fileW.write(match.group() + "\n")
        index += 1

fileW.close()


with open("GBI_IDS") as f:
    ids = f.read().splitlines()


filename = "Gb_V4.onlyGenes.sorted2_edited99.gff_re"
fileW = open(filename, "w")

with fileinput.FileInput("Gb_V4.onlyGenes.sorted2_edited99.gff") as file:
    for line in file:
        match  = re.findall(r"GBI_[0-9]+-R[A-Z]",line)
        if match:    
            if match[0] not in ids:
                fileW.write(line)
        else:
            fileW.write(line)

fileW.close()






#---------------------------------------------------------
filename = "slines"
fileW = open(filename, "w")

with fileinput.FileInput("duplicates2.txt") as file:
    for line in file:
        match  = re.search(r"(JP0080:WAR:STX:ANN:Lines )(\[.+?\])( and )(\[.+?\])",line)
        if match:
            match = str(match.group(4))
            match = match.replace(u'\xa0',u'')
            match = match.replace(u'[',u'')
            match = match.replace(u']',u'')
            fileW.write(str(int(match)) + "\n")

fileW.close()


with open("slines") as f:
    sl = f.read().splitlines()
sl = [int(s) for s in sl]

index = 1

filename = "GBI_IDS"
fileW = open(filename, "w")

with fileinput.FileInput("mss_duplicates_11_re") as file:
    for line in file:
        if index in sl:
            print(index)
            match  = re.search(r"GBI_[0-9]+-R[A-Z]",line)
            fileW.write(match.group() + "\n")
            print(line)
        index += 1

fileW.close()


with open("GBI_IDS") as f:
    ids = f.read().splitlines()


filename = "Gb_V4.onlyGenes.sorted2_edited99.gff_re_re"
fileW = open(filename, "w")

with fileinput.FileInput("Gb_V4.onlyGenes.sorted2_edited99.gff_re") as file:
    for line in file:
        match  = re.findall(r"GBI_[0-9]+-R[A-Z]",line)
        if match:    
            if match[0] not in ids:
                fileW.write(line)
        else:
            fileW.write(line)

fileW.close()