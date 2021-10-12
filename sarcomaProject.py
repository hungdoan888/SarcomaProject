# -*- coding: utf-8 -*-
"""
Created on Sun Oct 10 20:14:41 2021

@author: hungd
"""

#%% imports

import pandas as pd
import os
import numpy as np

#%% Output Path

output = os.getcwd() + r"\SarcomaProject.xlsx"

#%% Read Files

def getFiles():
    myFiles = []
    root = os.getcwd()
    for path, subdirs, files in os.walk(root):
        for name in files:
            myFiles.append(os.path.join(path, name))
    return myFiles
        
#%% Clinical Data

def clinicalData(nameOfSarcoma, myFile):
           
    # Read Data
    df_cd = pd.read_csv(myFile)
    
    # Number of Samples
    numberOfSamples = df_cd["Patient ID"].count()
    
    # Number of Patients
    numberOfPatients = df_cd["Patient ID"].drop_duplicates().count()
    
    # Mutations <= 2
    mutationsLT2 = df_cd[df_cd["Mutation Count"] < 2]["Patient ID"].count()
    
    # Mutations 3-5
    mutations3To5 = df_cd[(df_cd["Mutation Count"] >= 3) & 
                          (df_cd["Mutation Count"] <= 5)]["Patient ID"].count()
    
    # Mutations >= 6
    mutationsGT6 = df_cd[df_cd["Mutation Count"] >= 6]["Patient ID"].count()
    
    # Rename age column
    df_cd = df_cd.rename(columns={"Age at Which Sequencing was Reported": "age"})
    
    # Get unique patients and age
    df_age = df_cd[["Patient ID", "age"]].drop_duplicates().reset_index()
    
    # Convert Age to number
    df_age.loc[df_age["age"] == "<18", "age"] = 17
    df_age.loc[df_age["age"] == ">89", "age"] = 90
    
    # Convert Age to a Number if it can, otherwise convert to -99
    for i in range(len(df_age)):
        try:
            df_age.loc[i, "age"] = int(df_age.loc[i, "age"])
        except:
            df_age.loc[i, "age"] = np.nan
    
    # Age < 18
    ageLT18 = df_age[df_age["age"] < 18]["Patient ID"].count()
    
    # Age 18-40
    age18To40 = df_age[(df_age["age"] >= 18) & 
                      (df_age["age"] <= 40)]["Patient ID"].count()
    
    # Age 41-60
    age41To60 = df_age[(df_age["age"] >= 41) & 
                      (df_age["age"] <= 60)]["Patient ID"].count()
    
    # Age > 60
    ageGT60 = df_age[df_age["age"] > 60]["Patient ID"].count()
    
    # Sex
    df_sex = df_cd[["Patient ID", "Sex"]].drop_duplicates().reset_index()
    numberOfMales = df_sex[df_sex["Sex"] == "Male"]["Patient ID"].count()
    numberOfFemales = df_sex[df_sex["Sex"] == "Female"]["Patient ID"].count()
    
    # Race (White, Black, Asian, Native American, Other)
    df_race = df_cd[["Patient ID", "Primary Race"]].drop_duplicates().reset_index()
    numberOfWhite = df_race[df_race["Primary Race"] == "White"]["Patient ID"].count()
    numberOfBlack = df_race[df_race["Primary Race"] == "Black"]["Patient ID"].count()
    numberOfAsian = df_race[df_race["Primary Race"] == "Asian"]["Patient ID"].count()
    numberOfNativeAmerican = df_race[df_race["Primary Race"] == "Native American"]["Patient ID"].count()
    numberOfRaceOther = df_race[df_race["Primary Race"] == "Other"]["Patient ID"].count()
    
    # Type of Cancer
    numberOfPrimary = df_cd[df_cd["Sample Type"] == "Primary"]["Patient ID"].count()
    numberOfMetastasis = df_cd[df_cd["Sample Type"] == "Metastasis"]["Patient ID"].count()
    
    # Create Table
    df_cd_out = pd.DataFrame({"Sarcoma Type": [nameOfSarcoma],
                              "Samples (n)": [numberOfSamples],
                              "Patients (n)": [numberOfPatients],
                              "Mutations (<=2)": [mutationsLT2],
                              "Mutations (3-5)": [mutations3To5],
                              "Mutations (>=6)": [mutationsGT6],
                              "Age (<18)": [ageLT18],
                              "Age (18-40)": [age18To40],
                              "Age (41-60)": [age41To60],
                              "Age (>=61)": [ageGT60],
                              "Male": [numberOfMales], 
                              "Female": [numberOfFemales],
                              "White": [numberOfWhite],
                              "Black": [numberOfBlack], 
                              "Asian": [numberOfAsian], 
                              "Native American": [numberOfNativeAmerican], 
                              "Other Race": [numberOfRaceOther],
                              "Primary": [numberOfPrimary],
                              "Metastasis": [numberOfMetastasis]})
    return df_cd_out
    
#%% Mutated Genes

def mutatedGenes(nameOfSarcoma, myFile):
    df_mg = pd.read_csv(myFile)
    df_mg["Sarcoma Type"] = nameOfSarcoma
    df_mg_out = df_mg[["Sarcoma Type", "Gene", "#"]].sort_values("#", ascending=False)
    return df_mg_out

#%% Structural Variant Genes

def structuralVariantGenes(nameOfSarcoma, myFile):
    df_sv = pd.read_csv(myFile)
    df_sv["Sarcoma Type"] = nameOfSarcoma
    df_sv_out = df_sv[["Sarcoma Type", "Gene", "#"]].sort_values("#", ascending=False)
    return df_sv_out

#%% CNA Genes

def CNAGenes(nameOfSarcoma, myFile):
    df_cna = pd.read_csv(myFile)
    df_cna["Sarcoma Type"] = nameOfSarcoma
    df_cna_out = df_cna[["Sarcoma Type", "Gene", "CNA", "#"]].sort_values("#", ascending=False)
    return df_cna_out

#%% Run Script

if __name__ == "__main__":
    
    # Define Tables
    df_cd_master = pd.DataFrame()
    df_mg_master = pd.DataFrame()
    df_sv_master = pd.DataFrame()
    df_cna_master = pd.DataFrame()
    
    # Get Files
    myFiles = getFiles()
    
    for i in range(len(myFiles)):
        print("Working on", myFiles[i])
        
        # Name of Sarcoma
        myFile = os.path.split(myFiles[i])[1]
        firstUnderscoreIdx = myFile.find("_")
        nameOfSarcoma = myFile[:firstUnderscoreIdx]
        
        # Clinical Data
        if myFile.endswith("_clinical_data.csv"):
            df_cd_out = clinicalData(nameOfSarcoma, myFiles[i])
            df_cd_master = pd.concat([df_cd_master, df_cd_out])
        # Mutated Ganes
        elif myFile.endswith("_Mutated_Genes.csv"):
            df_mg_out = mutatedGenes(nameOfSarcoma, myFiles[i])
            df_mg_master = pd.concat([df_mg_master, df_mg_out])
        # Structural Variant Genes.csv
        elif myFile.endswith("_Structural_Variant_Genes.csv"):
            df_sv_out = structuralVariantGenes(nameOfSarcoma, myFiles[i])
            df_sv_master = pd.concat([df_sv_master, df_sv_out])
        # CNA Genes
        elif myFile.endswith("_CNA_Genes.csv"):
            df_cna_out = CNAGenes(nameOfSarcoma, myFiles[i])
            df_cna_master = pd.concat([df_cna_master, df_cna_out])
        else:
            print("\tFile is not in proper format:", myFile)
            
    # Write to excel
    with pd.ExcelWriter(output) as writer:
        df_cd_master.to_excel(writer, sheet_name="Clinical Data", index=False)
        df_mg_master.to_excel(writer, sheet_name="Mutated Genes", index=False)
        df_sv_master.to_excel(writer, sheet_name="Structural Variant Genes", index=False)
        df_cna_master.to_excel(writer, sheet_name="CNA Genes", index=False)
    
    print("Complete")
    
            
            
            
    