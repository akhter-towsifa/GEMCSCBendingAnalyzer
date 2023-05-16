'''
This file reads a .list file or .txt file that already has name of files
line by line. All this file does is concatenates the string path name to
each file line by line.

--Towsifa Akhter
'''

f = open("Run2023B-MuAlCalIsolatedMu-PromptReco-v1-exp-2023-08-05.list", "r")

#print(f.read(10)) #print first 10 characters of the file

list_of_original_lines = f.readlines()
list_of_new_lines = ["root://cms-xrd-global.cern.ch/" + line.strip() +"\n" for line in list_of_original_lines]

f = open("Run2023B-MuAlCalIsolatedMu-PromptReco-v1-exp-2023-08-05.list", "w")
f.writelines(list_of_new_lines)
f.close()
