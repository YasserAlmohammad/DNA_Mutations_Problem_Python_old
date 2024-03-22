"""
read Nucleotides and combine them into Codons, with the Codon table at hand we use it to detect all Mutations in the sequence
the mutations have two types Insertion and PointBase, so here we will print all related statistics about these mutations

since it is unknown how many mutations or insertions the child sequnce could have we will have the following rules in identifying
the Insertion and Point Mutations:

PointBase Mutation is defined as: a child Codon with one different base from parent Codon that results in a different Codon
Insertion Mutation is defined as:  a child Codon with more than one different base, that results in another Codon, or any extra Codons
                                   in the DNA sequence found at the end of the child sequence.
"""
import os

# codon table
codonTable = {       'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
                     'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C',
                     'tta': 'L', 'tca': 'S', 'taa': '*', 'tga': '*',
                     'ttg': 'L', 'tcg': 'S', 'tag': '*', 'tgg': 'W',
                     'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R',
                     'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
                     'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R',
                     'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R',
                     'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S',
                     'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S',
                     'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R',
                     'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',
                     'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G',
                     'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G',
                     'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
                     'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G'
             }
"""
having two sample input files of the child and parent DNA sequence, we read the data into a dictionary strucutre
"""
def readSamples():
    sampleData={} #two arrays, child and parent
    childFile=open("ChildDna.txt","r")
    parentFile=open("ParentDna.txt","r")
    childData=[] #child array
    parentData=[] #parent array of dictionary values

    #for each line get child sample data
    for childLine in childFile:
        words=childLine.split(" ") #words seperated by spaces
        #first word is Sex, second Name, Third is Sequence
        childData.append({"Sex":words[0],"Name":words[1],"DNA":words[2]})
    #for each line get parent data
    for parentLine in parentFile: #tokenize the string line
        words=parentLine.split(" ")
        parentData.append({"Sex":words[0],"Name":words[1],"DNA":words[2]})

    #add data subsets to the data dictionary and return it
    sampleData["Child"]=childData
    sampleData["Parent"]=parentData

    #close files
    childFile.close()
    parentFile.close()
    print ("-------------------------------------------------------")
    print ("{0:7} {1:7} {2:8} {3}".format("Owner","Sex","Name","DNA Sequence"))
    print ("-------------------------------------------------------")
    for line in range(0,len(childData)):  #-1?
        childLine=childData[line]
        parentLine=parentData[line]
        print ("{0:7} {1:7} {2:8} {3}".format("CHILD:",childData[line]["Sex"],childData[line]["Name"],childData[line]["DNA"]))
        print ("{0:7} {1:7} {2:8} {3}".format("PARENT:",parentData[line]["Sex"],parentData[line]["Name"],parentData[line]["DNA"]))
    print ("-------------------------------------------------------")
    return sampleData

"""
 compare the child and parent DNA Sequence to see if there is any mutation,
 all mutations are extracted and stored in an array with dictionary elements with all sample's information
 only mutations are stored with the relative information, as dictionary values to allow searching and iteration
 to get statistics.
"""
def extractMutations():
    data=readSamples()
    childData=data["Child"]
    parentData=data["Parent"]
    i=0
    list=[] # array of dict of{name,sex,mutation index,desc} key is the name
    for line in range(0,len(childData)):  #-1?
        childLine=childData[line]
        parentLine=parentData[line]
        childDNA=childLine["DNA"]
        parentDNA=parentLine["DNA"]
        index=0 #codon index based on three nucls
        codonCharIndex=0 #character index in the dna sequence
        length=len(childDNA)
        plength=len(parentDNA)
        if length>plength:
            print("There Are [{0:2}] Extra Bases in the child [{1:8}] Sequence [{2}]".format(length-plength,childLine["Name"],childDNA[plength-1:length-1]))
        while codonCharIndex < length-1: #-1 or -2 will work since we get each three chars
            #two stages, first one match both sequences when they are equal, then continue for the child sequence alone
            #stage 1
            if codonCharIndex <plength-1:
                childCodon=childDNA[codonCharIndex:codonCharIndex+3] #codon = three character
                parentCodon=parentDNA[codonCharIndex:codonCharIndex+3]# parent codon
                codonCharIndex=codonCharIndex+3  #each 3 chars is a codon

                if childCodon!=parentCodon: #we have a mutation
                    mut={}  #mutations are stored as a dictionary added later as an array item
                    mut["Sex"]=childLine["Sex"]
                    mut["Name"]=childLine["Name"]
                    mut["Index"]=index
                    if codonTable[childCodon]!=codonTable[parentCodon]:  #verfify mutation type based on codon table
                        #check to see if only one base has changed
                        change=0;
                        if childCodon[0]!=parentCodon[0]:
                            change+=1
                        if childCodon[1]!=parentCodon[1]:
                            change+=1
                        if childCodon[2]!=parentCodon[2]:
                            change+=1
                        if change==1:
                            mut["Type"]="PointBase"
                        else:
                            mut["Type"]="Insertion"
                        mut["Description"]=parentCodon+"==>"+childCodon+"["+codonTable[parentCodon]+"==>"+codonTable[childCodon]+"]"
                        list.append(mut)
                index=index+1 #codon index increase by 1 only where are character index is by 3
            else:
                #extra codons)
                childCodon=childDNA[codonCharIndex:codonCharIndex+3] #codon = three character
                codonCharIndex=codonCharIndex+3
                mut={}  #mutations are stored as a dictionary added later as an array item
                mut["Sex"]=childLine["Sex"]
                mut["Name"]=childLine["Name"]
                mut["Index"]=index
                mut["Type"]="Insertion"
                mut["Description"]="["+childCodon+"]"
                list.append(mut)
                index=index+1
    return list #construct array contraining all mutations in the sample data


"""
extract the statistics needed and store it in a dictionary
"""
def calcTotalMutationStatistics(mutations):
    dictionary={}
    totalPoint=0
    totalInsertion=0
    totalPointMale=0
    totalPointFemale=0
    totalInsertionMale=0
    totalInsertionFemale=0
    for mutation in mutations:
       if mutation["Type"]=="PointBase":
           totalPoint=totalPoint+1
       if mutation["Type"]=="Insertion":
           totalInsertion=totalInsertion+1
       if (mutation["Type"]=="PointBase") & (mutation["Sex"]=="Male"):
            totalPointMale=totalPointMale+1
       if (mutation["Type"]=="PointBase") & (mutation["Sex"]=="Female"):
            totalPointFemale=totalPointFemale+1
       if (mutation["Type"]=="Insertion") & (mutation["Sex"]=="Male"):
            totalInsertionMale=totalInsertionMale+1
       if (mutation["Type"]=="Insertion") & (mutation["Sex"]=="Female"):
            totalInsertionFemale=totalInsertionFemale+1

    dictionary["TotalPoint"]=totalPoint
    dictionary["TotalInsertion"]=totalInsertion
    dictionary["TotalPointMale"]=totalPointMale
    dictionary["TotalPointFemale"]=totalPointFemale
    dictionary["TotalInsertionMale"]=totalInsertionMale
    dictionary["TotalInsertionFemale"]=totalInsertionFemale

    return dictionary

"""
    having our mutation data and family name we extract the statistics needed and store it in a dictionary
"""
def calcTotalFamilyMutation(mutations,familyName):
    totalPointForFamily=0
    totalInsertionForFamily=0
    totalMutationsForFamily=0
    dictionary={}
    for mutation in mutations:
        if (mutation["Type"]=="PointBase") & (mutation["Name"]==familyName):
            totalPointForFamily=totalPointForFamily+1
        if (mutation["Type"]=="Insertion") & (mutation["Name"]==familyName):
            totalInsertionForFamily=totalInsertionForFamily+1
    totalMutationsForFamily=totalPointForFamily+totalInsertionForFamily
    dictionary["TotalPointForFamily"]=totalPointForFamily
    dictionary["TotalInsertionForFamily"]=totalInsertionForFamily
    dictionary["TotalMutationsForFamily"]=totalMutationsForFamily
    return dictionary

"""
view Mutation information for a specific family
"""
def viewPointFamilyMutation(mutations,familyName):
    print ("{0:9} {1:9} {2:9} {3:9} {4}".format("Sex","Name","Index","Type","Description"))
    for mutation in mutations:
        if (mutation["Type"]=="PointBase") & (mutation["Name"]==familyName):
            print ("{0:9} {1:9} {2:9} {3:9} {4}".format(mutation["Sex"],mutation["Name"],mutation["Index"],mutation["Type"],mutation["Description"]))

"""
view Mutation information for a specific family
"""
def viewInsertionFamilyMutation(mutations,familyName):
    print ("{0:9} {1:9} {2:9} {3:9} {4}".format("Sex","Name","Index","Type","Description"))
    for mutation in mutations:
        if (mutation["Type"]=="Insertion") & (mutation["Name"]==familyName):
            print ("{0:9} {1:9} {2:9} {3:9} {4}".format(mutation["Sex"],mutation["Name"],mutation["Index"],mutation["Type"],mutation["Description"]))

"""
view Mutation information for a specific family
"""
def viewAllFamilyMutations(mutations,familyName):
    print ("{0:9} {1:9} {2:9} {3:9} {4}".format("Sex","Name","Index","Type","Description"))
    for mutation in mutations:
        if mutation["Name"]==familyName:
            print ("{0:9} {1:9} {2:9} {3:9} {4}".format(mutation["Sex"],mutation["Name"],mutation["Index"],mutation["Type"],mutation["Description"]))


muts=extractMutations()
statistics=calcTotalMutationStatistics(muts)  #all mutations in our sample data

try:

    while True:
        #test program choices for user input
        print("-------------------------------------------------------------------")
        print("                             Menu Options")
        print("-------------------------------------------------------------------")
        print("\tX: EXIT the program")
        print ("\t0: view mutation statistics")
        print ("\t1: Print Total Point mutations")
        print ("\t2: Print Total Insertion mutations")
        print ("\t3: Print Total Point mutations occurred for Male Child population")
        print ("\t4: Print Total Point mutations occurred for Female Child population")
        print ("\t5: Print Total Insertion mutations occurred for Male Child population")
        print ("\t6: Print Total Insertion mutations occurred for Female Child population")
        print ("\t7: Point mutation for any user inputted Family Name")
        print ("\t8: Insertion mutation for any user inputted Family Name")
        print ("\t9: All mutation for any user inputted Family Name ")

        option=input("input your choice:") #read user input
        if str(option)=="X" or str(option)=="x":
            break
        if str(option)=="0":
            print ("* Total Point Base Mutations:"+str(statistics["TotalPoint"]))
            print ("* Total Insertion Mutations:"+str(statistics["TotalInsertion"]))
            print ("* Total Point Base Male Mutations:"+str(statistics["TotalPointMale"]))
            print ("* Total Point Base Female Mutations:"+str(statistics["TotalPointFemale"]))
            print ("* Total Insertion Male Mutations:"+str(statistics["TotalInsertionMale"]))
            print ("* Total Insertion Female Mutations:"+str(statistics["TotalInsertionFemale"]))
        elif str(option)=="1":
            print ("Total Point Base Mutations ="+str(statistics["TotalPoint"]))
        elif str(option)=="2":
            print ("Total Insertion Mutations="+str(statistics["TotalInsertion"]))
        elif str(option)=="3":
            print ("Total Point Base Male Mutations="+str(statistics["TotalPointMale"]))
        elif str(option)=="4":
            print ("Total Point Base Female Mutations="+str(statistics["TotalPointFemale"]))
        elif str(option)=="5":
            print ("Total Insertion Male Mutations="+str(statistics["TotalInsertionMale"]))
        elif str(option)=="6":
            print ("Total Insertion Female Mutations="+str(statistics["TotalInsertionFemale"]))
        elif str(option)=="7":
            familyName=input("input family or child name:")
            familyData=calcTotalFamilyMutation(muts,familyName)
            print ("Total Point Base Mutations for "+familyName+" = "+str(familyData["TotalPointForFamily"]))
            viewPointFamilyMutation(muts,familyName)
        elif str(option)=="8":
            familyName=input("input family name:")
            calcTotalFamilyMutation(muts,familyName)
            familyData=calcTotalFamilyMutation(muts,familyName)
            print ("Total Insertion Mutations for "+familyName+" = "+str(familyData["TotalInsertionForFamily"]))
            viewInsertionFamilyMutation(muts,familyName)
        elif str(option)=="9":
            familyName=input("input family name:")
            calcTotalFamilyMutation(muts,familyName)
            familyData=calcTotalFamilyMutation(muts,familyName)
            print ("Total Mutations for "+familyName+" = "+str(familyData["TotalMutationsForFamily"]))
            viewAllFamilyMutations(muts,familyName)
        else:
            print ("wrong choice, check menu options, Case sensitive")
except:
    print ("Some Error Occured in the program")
