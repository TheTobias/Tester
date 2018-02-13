
from FunLover import *
import xlrd



############## Input complex in fasta ########################




#################### Isolate true peptide ##################################
# The peptide we use (isolated from complex)
#ThePeptideList = ['1g6r','4','A','-QSVTQPDARVTVSEGASLQLRCKYSYSATP-------YLFWYVQY-PRQGLQLLLKYYSGD----PVVQ-GVNGFEAEF-S-KSNSSFHLRKASVHWSDSAVYFCAVSGFA-----------SALTFGSGTKVIVLP','B','-AVTQSPRNKVAVTGGKVTLSCNQTNNHN------NMYWYRQDTGHGLRLIHYSYGAGS----TEKG--DIPD-GYKASRPSQENFSLILELATPSQTSVYFCASGGGGT----------------LYFGAGTRLSVL','M','H2-Kb','GPHSLRYFVTAVSRPGLGEPRYMEVGYVDDTEFVRFDSDAENPRYEPRARWMEQEG-PEYWERETQKAKGNEQSFRVDLRTLLGYYNQSKGGSHTIQVISGCEVGSDGRLLRGYQQYAYDGCDYIALNEDLKTWTAADMAALITKHKWEQAGEAER-LRAYLEGTCVEWLRRYLKNGNATL','P','SIYRYYGL']


def yeho():
    complexExcel = '/Users/TobiasHO/MasterPrograms/Anaconda3/MineScripts/database_tpm_final1.xlsx'
    complexworkbook = xlrd.open_workbook(complexExcel)
    complexworksheet = complexworkbook.sheet_by_index(0)

    ThePeptideList = complexworksheet.cell(2, 0).value
    ThePeptideList = ThePeptideList.split(',')



    a=PepSwitch1(ThePeptideList)
    peptide=(PepSwitch2(a[0], a[1], a[2], ThePeptideList))
    peptide1=(peptide[0][1],peptide[1][1],peptide[2][1])

    complexpath = '/Users/TobiasHO/MasterPrograms/Anaconda3/MineScripts/1g6r.fsa'
    # Open and Save complex file
    cf = open(complexpath, 'r')
    complex = cf.read()
    cf.close()

    # Switch
    PepSwitch3(complex,peptide1,complexpath)





def yeha():

    complexExcel = '/Users/TobiasHO/MasterPrograms/Anaconda3/MineScripts/database_tpm_final1.xlsx'
    complexworkbook = xlrd.open_workbook(complexExcel)
    complexworksheet = complexworkbook.sheet_by_index(0)

    #maxi = complexworksheet.nrows
    maxi = 20
    j=1
    for i in range(maxi):
        if (complexworksheet.cell(i, 0).value) != '':
            ThePeptideList = complexworksheet.cell(i, 0).value
            ThePeptideList = ThePeptideList.split(',')
            print(j)
            j+=1
            print(PepSwitch1(ThePeptideList))



def PepSwitch1(ThePeptideList):

#################### Find source peptide ##################################

# Define true peptide from input
    pep1=PeptideFun(name=ThePeptideList[0],peptide=ThePeptideList[10])

#Blasts a peptide and finds ID of the best match for source protein
    sourceP=pep1.blastIDfinder()

# Uses the found ID to find the source protein (saved in pep2)
    sourceP = sourceP.NCBIseqGet()


################### Use NetMHCpan4 to find peptides from source ################

# Define inputs (source protein, HLA type and length of cleaved peptides)
    HLAtype = ThePeptideList[7]
    pepLength = len(ThePeptideList[10])

    return(HLAtype,pepLength,sourceP)

# Get results from netMHCpan4
def PepSwitch2(HLAtype,pepLength,sourceP,ThePeptideList):

# Define true peptide from input
    pep1=PeptideFun(name=ThePeptideList[0],peptide=ThePeptideList[10])

# Run netMHCpan4 with inputs
    #print(HLAtype)
    #print(pepLength)
    #print(sourceP.p)

    df = '/Users/TobiasHO/MasterPrograms/Anaconda3/MineScripts/32124_NetMHCpan.xls'
    workbook = xlrd.open_workbook(df)
    worksheet = workbook.sheet_by_index(0)


############## Use BLOSUM62 to score every peptide against true one ###############
############## Pick the 3 most relevant ones ######################################

    bloresult = pep1.MHCpanBLOSUM(worksheet)


    bestresult = bestBLOSUM(bloresult[0],bloresult[1],3)

    return(bestresult)

########### Switch the found peptides out with the real one ####################
def PepSwitch3(complex,peptides,complexpath):


# Make amount of switches equal to amount of peptides you can switch
    for i in range (len(peptides)):
        names = '>P'+'{}'.format(i+1)
        pep1 = PeptideFun(name=names, peptide=peptides[i])

# Use the pepswitch function to do the switch
        newcomplex = pep1.pepswitch(complex)

# Define output path for the new complex
        new_path = complexpath.rstrip('.fsa')
        new_path = new_path + 'fake' + str(i+1) + '.fsa'

# Write new complex with switched peptide
        ncf = open(new_path, 'w')
        ncf.write('\n'.join(newcomplex))
        ncf.close()

