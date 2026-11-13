#PCAogThiagoScript.py (plus some perosnal edits)

# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 16:25:17 2024

@author: PEIXOTT
"""
import os
import argparse
import numpy as np


def execute(command, logFile=""):
    os.system(command)
    print(command)



def convertPfileToBfile(pfile, name, folder, plink2, logFile, setVarIds=True):
    outputPrefix = f"{folder}/{name}"
    varIdFlag = " --set-all-var-ids @:#:\\$r:\\$a --new-id-max-allele-len 96" if setVarIds else ""
    commandLine = f"{plink2} --pfile {pfile} --make-bed --out {outputPrefix} --double-id --keep-allele-order --chr 1-22{varIdFlag}"

    execute(commandLine)

    seen = set()
    duplicates = []
    with open(f"{outputPrefix}.bim") as bimFile:
        for line in bimFile:
            varID = line.strip().split()[1]
            if varID in seen:
                duplicates.append(varID)
            else:
                seen.add(varID)
    if duplicates:
        print(f"WARNING: {len(duplicates)} duplicate variant IDs found in {outputPrefix}.bim")
        for dup in duplicates[:10]:
            print(f"  Duplicate ID: {dup}")
        if len(duplicates) > 10:
            print(f"  ... and {len(duplicates) - 10} more")
    else:
        print(f"OK: No duplicate variant IDs found in {outputPrefix}.bim")

    return outputPrefix


def mergeCommonFiles(file1, file2, name, folder, plink1, logFile):
    outputPrefix = f"{folder}/{name}"
    commandLine = f"{plink1} --bfile {file1} --bmerge {file2} --make-bed --keep-allele-order --out {outputPrefix}"

    execute(commandLine)
    return outputPrefix

def extractVariants(bfile, extractFile, name, folder, plink1, logFile):
    outputPrefix = f"{folder}/{name}"
    commandLine = f"{plink1} --bfile {bfile} --make-bed --keep-allele-order --out {outputPrefix} --extract {extractFile}"

    execute(commandLine, logFile)

    return outputPrefix

def removeLDAndMAFToPCA(merged, target, ref, folder, name, plink1, logFile):
    outputPrefix = f"{folder}/{name}"
    command = f"{plink1} --bfile {merged} --indep-pairwise 200 50 0.2 --out {folder}/{name} --maf 0.01"
    execute(command, logFile)

    targetLD = extractVariants(target, f"{folder}/{name}.prune.in", f"{name}_Target", folder, plink1, logFile)
    if ref == "":
        return targetLD, ""

    refLD = extractVariants(ref, f"{folder}/{name}.prune.in", f"{name}_Ref", folder, plink1, logFile)

    return targetLD, refLD

def mergeRefAndTarget(bfileTarget, bfileRef, folder, name, plink1, logFile):
    print(f"Merging reference and target")

    execute(f"mkdir {folder}/{name}_MergeRefAlt")

    command = f"{plink1} --bfile {bfileTarget} --make-bed --keep-allele-order --out {folder}/{name}_MergeRefAlt/TargetMAF --maf 0.01"
    execute(command, logFile)
    execute(f"cp {folder}/{name}_MergeRefAlt/TargetMAF.bed {folder}/{name}_MergeRefAlt/TargetMAFID.bed")
    execute(f"cp {folder}/{name}_MergeRefAlt/TargetMAF.fam {folder}/{name}_MergeRefAlt/TargetMAFID.fam")


    execute(f"cp {bfileRef}.bed {folder}/{name}_MergeRefAlt/Ref.bed")
    execute(f"cp {bfileRef}.fam {folder}/{name}_MergeRefAlt/Ref.fam")


    print(f"Changing the bim ID for target file")
    #Open bim file to change varID to chrom:pos:A1:A2
    dictVar = {}
    bimFileIn = open(f"{folder}/{name}_MergeRefAlt/TargetMAF.bim")
    bimFileOut = open(f"{folder}/{name}_MergeRefAlt/TargetMAFID.bim", "w")
    for line in bimFileIn:
        chrom, ID, poscM, posBP, A1, A2 = line.strip().split()
        newID = f"{chrom}:{posBP}:{A1}:{A2}"
        bimFileOut.write(f"{chrom}\t{newID}\t{poscM}\t{posBP}\t{A1}\t{A2}\n")

        if newID not in dictVar:
            dictVar[newID] = True
        else:
            dictVar[newID] = False
    bimFileIn.close()
    bimFileOut.close()



    print(f"Changing the bim ID for ref file and generating the list of variants in common")
    # Open bim file to change varID to chrom:pos
    bimFileIn = open(f"{bfileRef}.bim")

    bimFileOut = open(f"{folder}/{name}_MergeRefAlt/Ref.bim", "w")
    inCommon = open(f"{folder}/{name}_MergeRefAlt/inCommon.txt", "w")
    for line in bimFileIn:
        chrom, ID, poscM, posBP, A1, A2 = line.strip().split()
        newID = f"{chrom}:{posBP}:{A1}:{A2}"
        bimFileOut.write(f"{chrom}\t{newID}\t{poscM}\t{posBP}\t{A1}\t{A2}\n")

        if newID in dictVar:
            #print(f"In common {newID}")
            if dictVar[newID]:
                #print(f"not added because it was false")
                inCommon.write(f"{newID}\n")

    inCommon.close()
    bimFileIn.close()
    bimFileOut.close()

    refFile = f"{folder}/{name}_MergeRefAlt/Ref"
    targetFile = f"{folder}/{name}_MergeRefAlt/TargetMAFID"
    commonList = f"{folder}/{name}_MergeRefAlt/inCommon.txt"
    targetCommon = extractVariants(targetFile, commonList, f"Target_common", f"{folder}/{name}_MergeRefAlt/", plink1, logFile)
    refCommon = extractVariants(refFile, commonList, f"Ref_common", f"{folder}/{name}_MergeRefAlt/", plink1, logFile)

    return mergeCommonFiles(targetCommon, refCommon, "Merged", f"{folder}/{name}_MergeRefAlt/", plink1, logFile), targetCommon, refCommon


def getProjectedPCA(target, reference, gcta, X, threads, name, folder, plink1, removeList, logFile):
    command = f"{gcta} --bfile {reference} --maf 0.01 --make-grm --out {folder}/{name} --thread-num {threads}"
    if X:
        command = f"{command} --chr 23 --autosome-num 25"
    else:
        command = f"{command} --autosome"

    if removeList != "":
        command = f"{command} --remove {removeList}"
    execute(command, logFile)

    command = f"{gcta} --grm {folder}/{name} --maf 0.01 --pca 50 --out {folder}/{name}_PCs --thread-num {threads}"
    if X:
        command = f"{command} --chr 23 --autosome-num 25"
    else:
        command = f"{command} --autosome"

    execute(command, logFile)

    command = f"{gcta} --bfile {reference} --maf 0.01 --pc-loading {folder}/{name}_PCs --out {folder}/{name}_VariantLoading --thread-num {threads}"
    if X:
        command = f"{command} --chr 23 --autosome-num 25"
    else:
        command = f"{command} --autosome"
    execute(command, logFile)
    if target == "":
        return f"{folder}/{name}_PCs.eigenvec"


    command = f"{gcta} --bfile {target} --maf 0.01 --project-loading {folder}/{name}_VariantLoading 50 --out {folder}/{name}_TargetPC --thread-num {threads}"
    if X:
        command = f"{command} --chr 23 --autosome-num 25"
    else:
        command = f"{command} --autosome"

    if removeList != "":
        command = f"{command} --remove {removeList}"
    execute(command, logFile)

    return f"{folder}/{name}_TargetPC.proj.eigenvec"


#======================================================================== Covar =========================================================================
def modelSelection(covar, modelSelection, folder, name, logFile):
    commandLine = f"Rscript {modelSelection} {covar} {folder}/{name}"

    execute(commandLine)

    print(f"Model variables file {folder}/{name}_variables.tsv")

def buildCovarFile(covarDict, remove, PCList, PCASource, bfile, folder, name):
    removeList = []
    if remove != "":
        removeFile = open(remove)
        for line in removeFile:
            FID, IID = line.strip().split()
            removeList.append(IID)
        removeFile.close()

    covarList = []
    famFile = open(f"{bfile}.fam")
    covarFile = open(f"{folder}/{name}.tsv", "w")

    #Header -> IID <covar non PC> <Covar PC in Source_PCNumber format>
    covarFile.write("IID")
    for ind in covarDict:
        for covar in covarDict[ind]:
            if covar != "PCA":
                covarFile.write(f"\t{covar}")
                covarList.append(covar)
        break
    for source in PCASource:
        for PC in PCList:
            covarFile.write(f"\t{source}_{PC}")
    covarFile.write("\n")

    for line in famFile:
        FID, IID, mother, father, sex, pheno = line.strip().split()
        if IID not in removeList and IID in covarDict:
            covarFile.write(f"{IID}")
            for covar in covarList:
                covarFile.write(f"\t{covarDict[IID][covar]}")
            for source in PCASource:
                for PC in PCList:
                    covarFile.write(f"\t{covarDict[IID]['PCA'][source][PC]}")
            covarFile.write("\n")
    covarFile.close()
    famFile.close()
    return f"{folder}/{name}.tsv"

def addPCAToCovarDict(filePCA, covarDict, dataSource):
    for ind in covarDict:
        #print(f"For error: {covarDict[ind]}")
        if "PCA" not in covarDict[ind]:
            covarDict[ind]["PCA"] = {}
        if dataSource not in covarDict[ind]["PCA"]:
            covarDict[ind]["PCA"][dataSource] = {}

    print(f"Open {filePCA} to add to covar")
    file = open(f"{filePCA}")

    #Project PCA has no header, while the previous PCA had
    PCList = []

    for line in file:
        #The projection is tab separated, the single PCA is space
        data = line.strip().split()
        print(line)
        ind = data[1]
        if ind in covarDict:
            for i in range(2, len(data)):
                PC = f"PC{i-1}"
                if PC not in PCList:
                    PCList.append(PC)
                covarDict[ind]["PCA"][dataSource][PC] = data[i]


    return covarDict, PCList


def readCovarFile(covarTable):
    print('We are reading the covar table. We are assuming that the ID is the first col')
    print('We are also assuming that there is the column SEX and the Phenotype column is named DISEASE')
    print('We are also checking if there is any covar field that is empty or NA')

    file = open(covarTable)

    covarDict = {}
    missing = []
    header = True

    # doNotFix is a flag to put the phenotype in PLINK2 format (1 control, 2 case)
    doNotFix = False
    for line in file:
        if header:
            header = False
            splitHeader = line.strip().split()
        else:
            split = line.strip().split()

            nonNA = True
            for i in range(len(split)):
                if split[i] == "NA" or split[i] == "" or split[i] == " " or split[i] == "nan":
                    nonNA = False
                    ind = split[0]
                    data = split[i]
                    headerName = splitHeader[i]

                    print(f'Removing the ind {ind} because there is missing data ({data}) on the field {headerName}')
                    missing.append(ind)

            if nonNA:
                covarDict[split[0]] = {}
                for i in range(1, len(split)):
                    if splitHeader[i].upper() == "DISEASE":
                        if split[i] == "0" or split[i] == 0:
                            split[i] = "1"
                        elif split[i] == "1" or split[i] == 1:
                            split[i] = "2"
                        elif split[i] == "2" or split[i] == 2:
                            split[i] = "3"
                            doNotFix = True
                        else:
                            print(f"Unknown status value to sample {split[0]}: {split[i]} ")
                    covarDict[split[0]][splitHeader[i].upper()] = split[i]

    if doNotFix:
        for sample in covarDict:
            if covarDict[sample]["DISEASE"] == "2":
                covarDict[sample]["DISEASE"] = "1"
            elif covarDict[sample]["DISEASE"] == "3":
                covarDict[sample]["DISEASE"] = "2"
            else:
                print(f"Unknown sample DISEASE field:  {covarDict[sample]['DISEASE']}")

    return covarDict, missing

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='PCA and regression')

    data = parser.add_argument_group("Data arguments")
    data.add_argument('-A', '--autosomal', help='Genotyped file name for autosomal chromosome', required=True)
    data.add_argument('-t', '--tableCovar', help='File with covariatives to be added to the model', required=True)
    data.add_argument('-r', '--remove', help='List of samples to be removed', required=False, default = "")
    data.add_argument('--threads', help='Number of processors to be used (default = 1)',
                      required=False, default = 1)

    refData = parser.add_argument_group("Reference data arguments")
    refData.add_argument('-R', '--AutosomalRef', help='Data with parental reference data to run autosomal PCA (optional)', required=False, default="")

    output = parser.add_argument_group("Output arguments")
    output.add_argument('-n', '--name', help='Analysis name', required=True)
    output.add_argument('-f', '--folder', help='Folder to output files', required=True)

    pcaArg = parser.add_argument_group("PCA arguments")
    pcaArg.add_argument('--model', help='Regression model. If you do not provide a model the script will'
                                        'use the selectModel.R to build the model', required=False, default="", nargs="+")


    programs = parser.add_argument_group("Programs")
    programs.add_argument('--plink2', help='Path of PLINK 2 (default = plink2)', required=False, default="plink2")
    programs.add_argument('--plink1', help='Path of PLINK 1 (default = plink)', required=False, default="plink")
    programs.add_argument('--gcta', help='Path of gcta', required=True)
    programs.add_argument('--selectModel', help='Path of selectModel script', required=False)

    args = parser.parse_args()
    logFile = ""

    covarDict, missing = readCovarFile(args.tableCovar)
    os.system(f"mkdir {args.folder}")
    
    if os.path.isfile(f"{args.autosomal}.pgen"):
        args.autosomal = convertPfileToBfile(args.autosomal, "Target", args.folder, args.plink2, logFile)


    if args.AutosomalRef != "":
        # Projected PCA
        if os.path.isfile(f"{args.AutosomalRef}.pgen"):
            args.AutosomalRef = convertPfileToBfile(args.AutosomalRef, "Ref", args.folder, args.plink2, logFile, setVarIds=False)


        bfileMerged, targetCommon, refCommon = mergeRefAndTarget(args.autosomal, args.AutosomalRef, args.folder, f"{args.name}_AutosomalPCA", args.plink1, logFile)
        targetLD, refLD = removeLDAndMAFToPCA(bfileMerged, targetCommon, refCommon, args.folder, f"{args.name}_LD", args.plink1, logFile)
        autosomalPCA = getProjectedPCA(targetLD, refLD, args.gcta, False, args.threads, "AutosomalPCA",
                                       args.folder, args.plink1, args.remove, logFile)
    else:
        #PCA without referece
        #Send ref as empty to cancel the LD removal for ref
        targetLD, refLD = removeLDAndMAFToPCA(args.autosomal , args.autosomal, "", args.folder, f"{args.name}_LD",args.plink1, logFile)

        #Send target as empty to perform PCA only on target (that will be the ref on the function)
        autosomalPCA = getProjectedPCA("", targetLD, args.gcta, False, args.threads, "AutosomalPCA",
                                       args.folder, args.plink1, args.remove, logFile)

    covarDict, PCList = addPCAToCovarDict(autosomalPCA, covarDict, "GCTA")
    
    PCASource = ["GCTA"]
    outlierList = []
    covar = buildCovarFile(covarDict, args.remove, PCList, PCASource, args.autosomal, args.folder, f"{args.name}")
    modelSelection(covar, args.selectModel, args.folder, f"{args.name}", logFile)