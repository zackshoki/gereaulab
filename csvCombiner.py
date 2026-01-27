# READ: all files must be placed in a folder called "toCombine", outputs a single file called "combinedHALOdata.csv"

import csv
import os
from os import listdir
from os.path import isfile, join

folderPath = os.path.dirname(__file__) + "/toCombine" # edit as needed
outputName = "combinedHALOdata.csv"

filesToCombine = [f for f in listdir(folderPath) if isfile(join(folderPath, f))] # get a list of the names of all csv files to combine


combinedData = []; 

for index, csvName in enumerate(filesToCombine): # iterate through the list of excel sheet names

    csvPath = folderPath + "/" + csvName # build path to each csv
    with open(csvPath, newline='', encoding='utf-8') as csvFile: 
        reader = csv.reader(csvFile, delimiter=",") # open and read a excel sheet, as a list of lists (representing the rows of the sheet)
        if index == 0:
            combinedData.insert(0, next(reader)) # set the header of the data 
        else: 
            next(reader)
        for row in reader:
            combinedData.append(row) # add a new row of data


with open(outputName, 'w', newline='') as csvfile: # writes data to csv
    writer = csv.writer(csvfile)
    writer.writerows(combinedData)
    print('success!')

