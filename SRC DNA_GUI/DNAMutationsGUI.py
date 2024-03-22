"""
having a codon table we extract all mutations from a DNA sequence.
the mutations we are looking at are: Point and Interstion mutations.
Point mutations happen when on base changes between two codons of parent and child
Insertion mutations happen when extra bases are inserted into the sequence or more than one base changes between parent and child

the sample data is altered by adding few more bases to show the results.
"""
import os
import tkinter
import tkinter.messagebox
import codontable

codonstable=codontable.codonstable
"""
we read sample data files, parent and child, and store the information into a dicitonary for each line
"""
def parseFile():
    sampleData={} 
    cFile=open("ChildDna.txt","r")
    pFile=open("ParentDna.txt","r")
    cData=[] #child
    pData=[] #parent

    #for each line get child information
    for cLine in cFile:
        words=cLine.split(" ") #words seperated by spaces
        #Sex Name Sequence
        cData.append({"SEX":words[0],"NAME":words[1],"SEQUENCE":words[2]})
    #for each line get parent information
    for pLine in pFile:
        words=pLine.split(" ")
        pData.append({"SEX":words[0],"NAME":words[1],"SEQUENCE":words[2]})

    sampleData["CHILD"]=cData
    sampleData["PARENT"]=pData

    cFile.close()
    pFile.close()

    return sampleData

"""
 compare child and parent sequences and extract mutations
"""
def findMutations():
    data=parseFile()
    childData=data["CHILD"]
    parentData=data["PARENT"]
    list=[] # array of dict of{name,sex,mutation index,desc} key is the name
    for line in range(0,len(childData)):  #-1?
        childLine=childData[line]
        parentLine=parentData[line]
        childDNA=childLine["SEQUENCE"]
        parentDNA=parentLine["SEQUENCE"]
        index=0 #codon index based on three bases
        codonCharIndex=0 #character index in the sequence
        length=len(childDNA)
        plength=len(parentDNA)
        if length>plength:
            print("{0:2} Extra Bases in the child: {1:8} Sequence at the end:  ({2})".format(length-plength,childLine["NAME"],childDNA[plength-1:length-1]))
        while codonCharIndex < length-1:
            if codonCharIndex <plength-1:
                childCodon=childDNA[codonCharIndex:codonCharIndex+3] #codon = three character
                parentCodon=parentDNA[codonCharIndex:codonCharIndex+3]# parent codon
                codonCharIndex=codonCharIndex+3  #each 3 chars is a codon

                if childCodon!=parentCodon: #a mutation
                    mut={}
                    mut["SEX"]=childLine["SEX"]
                    mut["NAME"]=childLine["NAME"]
                    mut["Index"]=index
                    if codonstable[childCodon]!=codonstable[parentCodon]:  #verfify mutation type
                        #check to see if only one base has changed
                        difference=0;
                        if childCodon[0]!=parentCodon[0]:
                            difference+=1
                        if childCodon[1]!=parentCodon[1]:
                            difference+=1
                        if childCodon[2]!=parentCodon[2]:
                            difference+=1
                        if difference==1:
                            mut["TYPE"]="POINT"
                        else:
                            mut["TYPE"]="INSERTION"
                        mut["DESCRIPTION"]=parentCodon+"==>"+childCodon+"["+codonstable[parentCodon]+"==>"+codonstable[childCodon]+"]"
                        list.append(mut)
                index=index+1 # 3 chars basis
            else:
                #remaining codons
                childCodon=childDNA[codonCharIndex:codonCharIndex+3] #codon=3 bases
                codonCharIndex=codonCharIndex+3
                mut={}  #mutations are stored as a dictionary added later as an array item
                mut["SEX"]=childLine["SEX"]
                mut["NAME"]=childLine["NAME"]
                mut["Index"]=index
                mut["TYPE"]="INSERTION"
                mut["DESCRIPTION"]="["+childCodon+"]"
                list.append(mut)
                index=index+1
    return list #construct array contraining all mutations in the sample data


"""
extract statistics from mutations, once time calculation
"""
def computeTotalMutationStatistics(mutations):
    dictionary={}
    totalPoint=0
    totalInsertion=0
    totalPointMale=0
    totalPointFemale=0
    totalInsertionMale=0
    totalInsertionFemale=0
    for mutation in mutations:
       if mutation["TYPE"]=="POINT":
           totalPoint=totalPoint+1
       if mutation["TYPE"]=="INSERTION":
           totalInsertion=totalInsertion+1
       if (mutation["TYPE"]=="POINT") & (mutation["SEX"]=="Male"):
            totalPointMale=totalPointMale+1
       if (mutation["TYPE"]=="POINT") & (mutation["SEX"]=="Female"):
            totalPointFemale=totalPointFemale+1
       if (mutation["TYPE"]=="INSERTION") & (mutation["SEX"]=="Male"):
            totalInsertionMale=totalInsertionMale+1
       if (mutation["TYPE"]=="INSERTION") & (mutation["SEX"]=="Female"):
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
def computeTotalFamilyMutation(mutations,familyName):
    totalPointForFamily=0
    totalInsertionForFamily=0
    totalMutationsForFamily=0
    dictionary={}
    for mutation in mutations:
        if (mutation["TYPE"]=="POINT") & (mutation["NAME"]==familyName):
            totalPointForFamily=totalPointForFamily+1
        if (mutation["TYPE"]=="INSERTION") & (mutation["NAME"]==familyName):
            totalInsertionForFamily=totalInsertionForFamily+1
    totalMutationsForFamily=totalPointForFamily+totalInsertionForFamily
    dictionary["TotalPointForFamily"]=totalPointForFamily
    dictionary["TotalInsertionForFamily"]=totalInsertionForFamily
    dictionary["TotalMutationsForFamily"]=totalMutationsForFamily
    return dictionary

"""
create a class for GUI set up and calls
"""
class AppWindow(tkinter.Tk):
    def __init__(self,parent):
        tkinter.Tk.__init__(self,parent)
        self.parent = parent
        self.initialize()

    #set up the GUI widgets and register their events
    def initialize(self):
        self.grid()

        label = tkinter.Label(self,text="Family Name:",anchor="w")
        label.grid(column=0,row=0,sticky='EW')

        #entry is a text field for holding the family name
        self.entryVariable = tkinter.StringVar()
        #self.entryVariable.set("Deep")
        self.entry = tkinter.Entry(textvariable=self.entryVariable)
        self.entry.grid(column=1,row=0,sticky='EW')
        self.entry.bind("<Return>", self.OnPressEnter)

        statisticsButton = tkinter.Button(self,text=u"  View Total Mutation(Point-Insertion) Statistics  ",command=self.OnStatisticsButtonClick)
        statisticsButton.grid(column=0,row=1,sticky='EW')

        familyMutationsButton = tkinter.Button(self,text=u"View Total Mutation(Point-Insertion) For Family",command=self.OnFamilyButtonClick)
        familyMutationsButton.grid(column=1,row=1,sticky='EW')

        familyPointMutationsButton = tkinter.Button(self,text=u"View Point Mutations For Family",command=self.OnFamilyPointButtonClick)
        familyPointMutationsButton.grid(column=0,row=2,sticky='EW')

        familyInsertionMutationsButton = tkinter.Button(self,text=u"View Insertion Mutations For Family",command=self.OnFamilyInsertionButtonClick)
        familyInsertionMutationsButton.grid(column=1,row=2,sticky='EW')

        familyStatisticsMutationsButton = tkinter.Button(self,text=u"View Mutations Details For Family",command=self.OnFamilyStatisticsButtonClick)
        familyStatisticsMutationsButton.grid(column=0,row=3,sticky='EW')

        clearTextAreaButton = tkinter.Button(self,text=u"Clear Text Area",command=self.OnclearTextAreaButtonClick)
        clearTextAreaButton.grid(column=1,row=3,sticky='EW')

        #==================
        self.totalPointMutationsVariable = tkinter.StringVar()
        self.totalPointMutationsEntry = tkinter.Entry(textvariable=self.totalPointMutationsVariable)
        self.totalPointMutationsEntry.grid(column=1,row=4,sticky='EW')

        totalPointMutationsButton = tkinter.Button(self,text=u"View Total Point Mutations",command=self.OnTotalPointButtonClick)
        totalPointMutationsButton.grid(column=0,row=4,sticky='EW')
        #===================
        self.totalInsertionMutationsVariable = tkinter.StringVar()
        self.totalInsertionMutationsEntry = tkinter.Entry(textvariable=self.totalInsertionMutationsVariable)
        self.totalInsertionMutationsEntry.grid(column=1,row=5,sticky='EW')

        totalInsertionMutationsButton = tkinter.Button(self,text=u"View Total Insertion Mutations",command=self.OnTotalInsertionButtonClick)
        totalInsertionMutationsButton.grid(column=0,row=5,sticky='EW')
        #===================
        self.totalPointMutationsMaleVariable = tkinter.StringVar()
        self.totalPointMutationsMaleEntry = tkinter.Entry(textvariable=self.totalPointMutationsMaleVariable)
        self.totalPointMutationsMaleEntry.grid(column=1,row=6,sticky='EW')

        totalPointMutationsMaleButton = tkinter.Button(self,text=u"View Total Point Mutations For Males",command=self.OnTotalPointMaleButtonClick)
        totalPointMutationsMaleButton.grid(column=0,row=6,sticky='EW')
        #====================
        self.totalPointMutationsFemaleVariable = tkinter.StringVar()
        self.totalPointMutationsFemaleEntry = tkinter.Entry(textvariable=self.totalPointMutationsFemaleVariable)
        self.totalPointMutationsFemaleEntry.grid(column=1,row=7,sticky='EW')

        totalPointMutationsFemaleButton = tkinter.Button(self,text=u"View Total Point Mutations For Females",command=self.OnTotalPointFemaleButtonClick)
        totalPointMutationsFemaleButton.grid(column=0,row=7,sticky='EW')
        #====================
        self.totalInsertionMutationsMaleVariable = tkinter.StringVar()
        self.totalInsertionMutationsMaleEntry = tkinter.Entry(textvariable=self.totalInsertionMutationsMaleVariable)
        self.totalInsertionMutationsMaleEntry.grid(column=1,row=8,sticky='EW')

        totalInsertionMutationsMaleButton = tkinter.Button(self,text=u"View Total Insertion Mutations For Males",command=self.OnTotalInsertionMaleButtonClick)
        totalInsertionMutationsMaleButton.grid(column=0,row=8,sticky='EW')
        #====================
        self.totalInsertionMutationsFemaleVariable = tkinter.StringVar()
        self.totalInsertionMutationsFemaleEntry = tkinter.Entry(textvariable=self.totalInsertionMutationsFemaleVariable)
        self.totalInsertionMutationsFemaleEntry.grid(column=1,row=9,sticky='EW')

        totalInsertionMutationsFemaleButton = tkinter.Button(self,text=u"View Total Insertion Mutations For Females",command=self.OnTotalInsertionFemaleButtonClick)
        totalInsertionMutationsFemaleButton.grid(column=0,row=9,sticky='EW')
        #====================


        #setup text area at the end to view results
        self.resultsText = tkinter.StringVar()
        self.results=tkinter.Text(self,height=20,width="50")
        self.results.grid(column=0,row=10,columnspan=2,sticky='EW')

        #results.insert(1,"hello")

        #setup a scroll bar for the text

        scroll = tkinter.Scrollbar(self)
        scroll.grid(column=1,row=10,sticky='E')
        #scroll.pack(side=tkinter.RIGHT, fill=tkinter.Y)
        #self.results.pack(side=tkinter.LEFT, fill=tkinter.BOTH)
        ##################################################################
        scroll.config(command=self.results.yview)
        self.results.config(yscrollcommand=scroll.set)


        self.grid_columnconfigure(1,weight=1) #resize Second column when window is resized
        self.resizable(True,False) #vertically and horizontally
        self.update()
        self.geometry(self.geometry()) #dont change window size when input is changed
        self.entry.focus_set()
        self.entry.selection_range(0, tkinter.END)

    def OnPressEnter(self,event):
        if self.entryVariable.get().strip()=="":
            tkinter.messagebox.showinfo("Info","You have to enter Family Name first")
        #self.labelVariable.set( self.entryVariable.get()+" (You clicked the button)" )
    def OnStatisticsButtonClick(self):
        self.results.insert(tkinter.END, "==============================\n")
        self.results.insert(tkinter.END,"- Total Point Base Mutations:"+str(mutationResults["TotalPoint"])+"\n")
        self.results.insert(tkinter.END,"- Total Insertion Mutations:"+str(mutationResults["TotalInsertion"])+"\n")
        self.results.insert(tkinter.END,"- Total Point Base Male Mutations:"+str(mutationResults["TotalPointMale"])+"\n")
        self.results.insert(tkinter.END,"- Total Point Base Female Mutations:"+str(mutationResults["TotalPointFemale"])+"\n")
        self.results.insert(tkinter.END,"- Total Insertion Male Mutations:"+str(mutationResults["TotalInsertionMale"])+"\n")
        self.results.insert(tkinter.END,"- Total Insertion Female Mutations:"+str(mutationResults["TotalInsertionFemale"])+"\n")

        #self.labelVariable.set( self.entryVariable.get()+" (You pressed ENTER)" )
    def OnFamilyButtonClick(self):
        if self.entryVariable.get().strip()=="":
            tkinter.messagebox.showinfo("Info","You have to enter Family Name first")
        else:
            familyName=self.entryVariable.get()
            data=computeTotalFamilyMutation(mutations,familyName)
            self.results.insert(tkinter.END, "==============================\n")
            self.results.insert(tkinter.END, "* Total Point Base Mutations for "+familyName+" :"+str(data["TotalPointForFamily"])+"\n")
            self.results.insert(tkinter.END,"* Total Insertion Mutations for "+familyName+" :"+str(data["TotalInsertionForFamily"])+"\n")
            self.results.insert(tkinter.END,"* Total Mutations for "+familyName+" :"+str(data["TotalMutationsForFamily"])+"\n")

    def OnFamilyPointButtonClick(self):
        if self.entryVariable.get().strip()=="":
            tkinter.messagebox.showinfo("Info","You have to enter Family Name first")
        else:
            familyName=self.entryVariable.get()
            self.results.insert(tkinter.END, "==============================\n")
            self.results.insert(tkinter.END,"{0:9} {1:9} {2:9} {3:9} {4}\n".format("SEX","NAME","Index","TYPE","DESCRIPTION"))
            for mutation in mutations:
                if (mutation["TYPE"]=="POINT") & (mutation["NAME"]==familyName):
                    self.results.insert(tkinter.END,"{0:9} {1:9} {2:9} {3:9} {4}\n".format(mutation["SEX"],mutation["NAME"],mutation["Index"],mutation["TYPE"],mutation["DESCRIPTION"]))

    def OnFamilyInsertionButtonClick(self):
        if self.entryVariable.get().strip()=="":
            tkinter.messagebox.showinfo("Info","You have to enter Family Name first")
        else:
            familyName=self.entryVariable.get()
            self.results.insert(tkinter.END, "==============================\n")
            self.results.insert(tkinter.END,"{0:9} {1:9} {2:9} {3:9} {4}\n".format("SEX","NAME","Index","TYPE","DESCRIPTION"))
            for mutation in mutations:
                if (mutation["TYPE"]=="INSERTION") & (mutation["NAME"]==familyName):
                    self.results.insert(tkinter.END,"{0:9} {1:9} {2:9} {3:9} {4}\n".format(mutation["SEX"],mutation["NAME"],mutation["Index"],mutation["TYPE"],mutation["DESCRIPTION"]))

    def OnFamilyStatisticsButtonClick(self):
        if self.entryVariable.get().strip()=="":
            tkinter.messagebox.showinfo("Info","You have to enter Family Name first")
        else:
            familyName=self.entryVariable.get()
            self.results.insert(tkinter.END, "==============================\n")
            self.results.insert(tkinter.END,"{0:9} {1:9} {2:9} {3:9} {4}\n".format("SEX","NAME","Index","TYPE","DESCRIPTION"))
            for mutation in mutations:
                if mutation["NAME"]==familyName:
                    self.results.insert(tkinter.END,"{0:9} {1:9} {2:9} {3:9} {4}\n".format(mutation["SEX"],mutation["NAME"],mutation["Index"],mutation["TYPE"],mutation["DESCRIPTION"]))
            
    def OnclearTextAreaButtonClick(self):
        self.results.delete(1.0,tkinter.END)
    def setVars(self,m,s):
        self.mutations=m
        self.mutationResults=s

    def OnTotalPointButtonClick(self):
            self.totalPointMutationsVariable.set(mutationResults["TotalPoint"])
    def OnTotalInsertionButtonClick(self):
            self.totalInsertionMutationsVariable.set(mutationResults["TotalInsertion"])
    def OnTotalPointMaleButtonClick(self):
            self.totalPointMutationsMaleVariable.set(mutationResults["TotalPointMale"])
    def OnTotalPointFemaleButtonClick(self):
            self.totalPointMutationsFemaleVariable.set(mutationResults["TotalPointFemale"])
    def OnTotalInsertionMaleButtonClick(self):
            self.totalInsertionMutationsMaleVariable.set(mutationResults["TotalInsertionMale"])
    def OnTotalInsertionFemaleButtonClick(self):
            self.totalInsertionMutationsFemaleVariable.set(mutationResults["TotalInsertionFemale"])


mutations=findMutations()
mutationResults=computeTotalMutationStatistics(mutations)  #all mutations in our sample data

def setupGUI():
    app = AppWindow(None)
    app.title('DNA Mutations')
    #app.setVars(mutations,mutationresults)
    app.mainloop()

setupGUI()
