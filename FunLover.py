

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.SubsMat import MatrixInfo
blosum = MatrixInfo.blosum62


from Bio import Entrez
Entrez.email = 'ohsaibot@mail.com'

import random


class PeptideFun(object):

    def __init__(self, name='', peptide='', proID=''):

        self.n = name
        self.p = peptide
        self.id = proID

############## Blasts a peptide and finds ID of the best match for source protein ######
############ Only works for small peptide ############################
    def blastIDfinder(self):

        # Perform BLASTp on peptide (Adjusted for small peptides)
        blast_handle = NCBIWWW.qblast('blastp', 'nr', self.p, matrix_name='PAM30', gapcosts='9 1', expect=200000,
                                      word_size=2, perc_ident=100)

        pep2=PeptideFun('','')

        # Read and find peptides
        records = NCBIXML.parse(blast_handle)
        item = next(records)
        for alignment in item.alignments:
            for hsp in alignment.hsps:

                # Decide which source peptide to use
                if hsp.align_length == len(self.p):  # Alignment has no gaps
                    if hsp.identities == len(self.p):  # Every aa match (doesn't account for gaps)
                        if alignment.length > 100:  # Picks relevant source peptides (not too small)
                            IDID = alignment.hit_id
                            pep2.n = alignment.hit_def
                            IDID = IDID.split('|')
                            pep2.id = IDID[1]
                            return pep2

        print('No source protein with 100% identity and length over 100 was found')
        return pep2

        # Define output path # This part is only here to check middle output is correct
        ##name = self.n  # here we can decide what parts of the name is included (fasta > can't be there)
        ##new_path = '/Users/TobiasHO/MasterPrograms/Anaconda3/MineScripts/Blast'
        ##new_path = new_path + name + '.xml'
        # Save FULL BLASTp result (Not necessary - can be deleted)
        ##blast_handle.seek(0)
        ##blast_file = open(new_path, 'w')
        ##blast_file.write(blast_handle.read())
        ##blast_file.close()

############### Switches your peptide out with another in a given complex ###########
    def pepswitch(self, complex):


        # Split complex file and peptide file and insert the peptide
        newcomplex = complex.split()
        newcomplex[2] = self.n
        newcomplex[3] = self.p

        return newcomplex

#################### Trade ID number for protein sequence from NCBI #######################
    def NCBIseqGet(self):

        # Float where output is saved
        pep2 = PeptideFun('', '')

        if self.id == '':
            print('You need an ID in self.ID to use the NCBIseqGet function')
            return pep2



        # The use of Entrez to find and get sequence for ID
        request = Entrez.epost("protein", id=self.id)
        result = Entrez.read(request)
        webEnv = result["WebEnv"]
        queryKey = result["QueryKey"]
        handle = Entrez.efetch(db="protein", retmode="xml", webenv=webEnv, query_key=queryKey)
        for r in Entrez.parse(handle):
            # Grab the GI
            try:
                gi = int([x for x in r['GBSeq_other-seqids'] if "gi" in x][0].split("|")[1])
            except ValueError:
                gi = None
            pep2.id = gi
            pep2.p = r["GBSeq_sequence"]

            return pep2

############### Finds peptides with binding affinity and gives blosum62#############
###compared to original peptide we are looking for (uses netMHCpan4 data as input)######
    def MHCpanBLOSUM(self, worksheet):

        score = []
        vals = []
        fakepep = []
        maxi = worksheet.nrows

        for i in range(maxi):

            if worksheet.cell(i, 9).value == 1 and str(worksheet.cell(i, 1).value) != self.p:
                seq2 = str(worksheet.cell(i, 1).value)
                score.append(nogapscore(self.p, seq2, blosum))
                fakepep.append(seq2)

                if float(worksheet.cell(i, 7).value) > 1:
                    vals.append(float(worksheet.cell(i, 7).value / 10000))
                    i += 1

                else:
                    vals.append(float(worksheet.cell(i, 7).value))
                    i += 1

            else:
                i += 1


        return(score, fakepep)

##################### Extra functions ##################
# Score a pair
def score_match(pair, matrix):
    if pair not in matrix:
        return matrix[(tuple(reversed(pair)))]
    else:
        return matrix[pair]

# Score two sequences with no gaps
def nogapscore(seq1,seq2,matrix):
    score = 0

    for i in range(len(seq1)):
        pair = (seq1[i], seq2[i])
        score+=score_match(pair, matrix)

    return score

# Choose best BLOSUM scores and order them
def bestBLOSUM(blosumscore,peptides,numberofpeps):

    bests = sorted(zip(blosumscore, peptides), reverse=True)[:numberofpeps]

    return bests