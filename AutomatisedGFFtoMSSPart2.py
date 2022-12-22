#!/usr/bin/python
# coding: UTF-8

import fileinput
import re
import sys

#------------------------------------------
# Part 2:  fixing translation errors -> fixing syntax erros
# Translation errors

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



args = sys.argv
transChecker_errors = args[1] # output of transChecker
mss_file = args[2] # converted mss file 


filename = "what_to_fix"    # generated file with information how to solve things
fileW = open(filename, "w")

with open(transChecker_errors) as in_file:
    for line in in_file:
        if line[0] == '>': # new location of error
            location = line.split()[1]
            oldLocation = location
            while line[0:2] != "//": # until errors for this location dont end
                line = in_file.readline()
                if (line.split()[0] == "TC0017:ER2:"): # first codon error
                    location = first(location)
                if (line.split()[0] == "TC0018:ER2:"): # final codon error
                    location = final(location)
            uninterrupted = oldLocation + ";" + location + "\n" # write to fix file in "bad_location;location_fix" manner
            fileW.write(uninterrupted)
fileW.close()


with fileinput.FileInput(mss_file, inplace=True) as out_file:
    x = 0 # position in replacement file
    for line in out_file: # for line in mss file
            splitted = re.split(r'\t+', line) # splitted line by tabs
            flag = 0 # flag to check if replacement in the line is being made
            if splitted[1] == "CDS": # if given line is CDS of mRNA header
                with open(filename) as in_file: # with open file with things that need to be fixed
                    for linee in in_file.readlines()[x:]: # opening the file from a new position as to not overwrite the fixed errors
                        replacement = linee.split(";") # version giving an error and the new version are separated my a ';' in a file
                        r1 = replacement[0] # old version
                        r2 = replacement[1] # new version
                        if r1 == splitted[2]: # if location is approptiate (location in the file = the one to be replaced)
                            x += 1 # new file opening position
                            flag = 1 # replacement needs to be done                           
                            savedR1 = r1.rstrip('\n') # stripping the line of the newline     
                            savedR2 = r2.rstrip('\n')                                         
                            break;
            if flag == 1 : print(line.replace(savedR1,savedR2), end='') # replacement to be done ###
            else :  print(line, end = '') # no replacement


#---------------------------------------------------
# Syntax errors - some of the common issues with DDBJ file format that can be quickly fixed


with fileinput.FileInput(mss_file, inplace = True) as file:
    for line in file:
        line = re.sub(r"\(\w+\s\w+.+\)","",line) # remove species names
        line = re.sub(r"\\","",line) # remove backslashes
        
        
        """
        match  = re.search(r"(Scaffold)([0-9]+)",line)  # used to changing scaffold names to DDBJ names
        if match:
            replacement = "BOPP0" + str(int(match.group(2)) + 1000000) 
            line = re.sub(match.group(0),replacement,line) # replace "scaffold" with names appointed by the DDBJ
        """
        
        """ # can be used for quick fixes
        match  = re.search(r"GBIMM",line) 
        if match:
            replacement = "GBIM" 
            line = re.sub(match.group(0),replacement,line) 

        
        match  = re.search(r" Gryllus bimacaulatus",line) 
        if match:
            replacement = "Gryllus bimacaulatus" 
            line = re.sub(match.group(0),replacement,line) 
        """
        """
        line  = re.sub(r"Similar to ", "", line) 
        """
        """
        match  = re.search(r"NCBI",line) 
        if match:
            replacement = "[NCBI]" 
            line = re.sub(match.group(0),replacement,line) 
        """


        print(line, end = '')
        
        
#-----------------  extend the correction to mRNA

# Open the input file in read mode
with open(mss_file, "r") as f:  
    # Read all the lines in the file
    g = f.readlines()

# Initialize a flag to False
flag = False

# Open the output file in write mode
with open("changed_ann.ann", 'w') as output:
    # Iterate through all the lines in the input file
    for line in g:
        # Search for a specific pattern in the line
        a = re.search("(\d+\.\.)(>)(\d+)", line)

        # If the pattern is found, set the flag to True
        if a:
            reference = line
            flag = True

        # If the flag is True and the line contains the string "mRNA\t"
        elif flag and re.search("mRNA\t", line):
            # Find all substrings that match the pattern "(\d+\.\.)(\d+)" in the line
            change = re.findall("(\d+\.\.)(\d+)", line)
            # Remove the last element from the list and assign it to the variable 'pop'
            pop = change.pop()
            # Assign the values of the tuple to variables 'old' and 'new'
            old = pop[0] + pop[1]
            new = pop[0] + ">" + pop[1]
            # Replace the old string with the new string in the line
            line = re.sub(old, new, line)
            # Set the flag to False
            flag = False

        # Write the modified line to the output file
        output.write(line)

# Open the intermediate output file in read mode
with open("changed_ann.ann", "r") as f:  
    # Read all the lines in the file
    g = f.readlines()

# Initialize a flag to False
flag = False

# Open the final output file in write mode
with open("final_ann.ann", 'w') as output:
    # Iterate through all the lines in the intermediate output file
    for line in g:
        # Search for a specific pattern in the line
        a = re.search("(<)(\d+\.\.)(>?)(\d+)", line)

        # If the pattern is found, set the flag to True
        if a:
            flag = True

        # If the flag is True and the line contains the string "mRNA\t"
        elif flag and re.search("mRNA\t", line):
            # Find all substrings that match the pattern "(\d+\.\.)(>?\d+)" in the line
            change = re.findall("(\d+\.\.)(>?\d+)", line)
            # Get the first element from the list and assign it to the variable 'pop'
            pop = change[0]
            # Assign the values of the tuple to the variable 'old'
            old = pop[0] + pop[1]
            # Create a new string by adding the character '<' to the beginning of 'old'
            new = "<" + old
            # Replace the old string with the new string in the line
            line = re.sub(old, new, line)
            # Set the flag to False
            flag = False

        # Write the modified line to the final output file
        output.write(line)


