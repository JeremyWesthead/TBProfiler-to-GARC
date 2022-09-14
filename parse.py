'''Parse the TBProfiler DB to GARC for use with `piezo`
'''
import re
import pickle
from datetime import date

import gumpy
import pandas as pd


def parseResistant(resistant: pd.DataFrame, reference: gumpy.Genome) -> dict:
    '''Parse the resistant CSV

    Args:
        resistant (pd.DataFrame): Dataframe containing the contents of `tbdb.csv`
        reference (gumpy.Genome): Genome object of the reference

    Returns:
        dict: Dictionary mapping (mutation in GARC, drug) --> 'R'
    '''
    garc = {}
    drugDict = {
            'amikacin': 'AMI', 'aminoglycosides': ('AMI', 'KAN', 'STM'), 'bedaquiline': 'BDQ', 'capreomycin': 'CAP',
            'ciprofloxacin': 'CIP', 'clofazimine': 'CFZ', 'cycloserine': 'CYC', 'delamanid': 'DLM', 
            'ethambutol': 'EMB', 'ethionamide': 'ETO', 'fluoroquinolones': 'FQS', 'isoniazid': 'INH', 
            'kanamycin': 'KAN', 'levofloxacin': 'LEV', 'linezolid': 'LZD', 'moxifloxacin': 'MXF', 
            'ofloxacin': 'OFX', 'para-aminosalicylic_acid': 'PAS', 'pyrazinamide': 'PZA', 
            'rifampicin': 'RIF', 'streptomycin': 'STM'
        }

    for _, row in resistant.iterrows():
        mutation = row['Mutation']
        gene = row['Gene']
        drug = drugDict[row['Drug']]

        if mutation[0] in ['n', 'c']:
            m = nucleotideVariants(gene, mutation)
            garc[(m, drug)] = 'R'

        elif mutation[0] == "p":
            m = proteinMutation(gene, mutation)
            garc[(m, drug)] = 'R'

        elif "any_missense_codon" in mutation:
            pos = mutation.split("_")[-1]
            g = reference.build_gene(gene)
            m = gene + "@" + g.amino_acid_sequence[g.amino_acid_number == int(pos)][0] + pos + "?"
            garc[(m, drug)] = 'R'
        
        elif mutation == "frameshift":
            m = gene + "@*_fs"
            garc[(m, drug)] = 'R'

        else:
            print("Skipping: ", gene, mutation, drug)

    #Checking for aminoglycosides tuple and expanding as required
    toRemove = []
    toAdd = {}
    for mutation, drug in garc.keys():
        if isinstance(drug, tuple):
            for d in drug:
                toAdd[(mutation, d)] = garc[(mutation, drug)]
            toRemove.append((mutation, drug))
    for r in toRemove:
        del garc[r]
    garc = {**garc, **toAdd}
    return garc

def parseOther(other: pd.DataFrame, reference: gumpy.Genome) -> dict:
    '''Parse the other annotations CSV

    Args:
        other (pd.DataFrame): Dataframe containing the contents of `tbdb.other_annotations.csv`
        reference (gumpy.Genome): Genome object of the reference

    Returns:
        dict: Dictionary mapping (mutation in GARC, drug) --> prediction
    '''
    garc = {}
    drugDict = {
        'amikacin': 'AMI', 'aminoglycosides': ('AMI', 'KAN', 'STM'), 'bedaquiline': 'BDQ', 'capreomycin': 'CAP',
        'ciprofloxacin': 'CIP', 'clofazimine': 'CFZ', 'cycloserine': 'CYC', 'delamanid': 'DLM', 
        'ethambutol': 'EMB', 'ethionamide': 'ETO', 'fluoroquinolones': 'FQS', 'isoniazid': 'INH', 
        'kanamycin': 'KAN', 'levofloxacin': 'LEV', 'linezolid': 'LZD', 'moxifloxacin': 'MXF', 
        'ofloxacin': 'OFX', 'para-aminosalicylic_acid': 'PAS', 'pyrazinamide': 'PZA', 
        'rifampicin': 'RIF', 'streptomycin': 'STM'
        }
    confidenceToVal = {
        'Assoc w R': 'R', 'Assoc w R - interim': 'R', 
        'Not assoc w R': 'S', 'Not assoc w R - Interim': 'S',
        'Uncertain significance': 'U'
    }
    for _, row in other.iterrows():
        info = row['Info']
        gene = row['Gene']
        mutation_ = row['Mutation']

        #Pull the values out of the info column
        values = {}
        for val in info.split(";"):
            key, v = val.split("=")
            values[key] = v


        #We can only add a prediction if we have both a drug and a prediction
        if values.get('drug') and values.get('who_confidence'):
            drug = drugDict[values.get('drug')]
            #We have both so we can continue
            for mutation in mutation_.split("|"):
                if mutation[0] in ['n', 'c']:
                    m = nucleotideVariants(gene, mutation)
                    garc[(m, drug)] = confidenceToVal[values.get('who_confidence')]

                elif mutation[0] == 'p':
                    m = proteinMutation(gene, mutation)
                    garc[(m, drug)] = confidenceToVal[values.get('who_confidence')]

                elif "any_missense_codon" in mutation:
                    pos = mutation.split("_")[-1]
                    g = reference.build_gene(gene)
                    m = gene + "@" + g.amino_acid_sequence[g.amino_acid_number == int(pos)][0] + pos + "?"
                    garc[(m, drug)] = confidenceToVal[values.get('who_confidence')]
                
                elif mutation == "frameshift":
                    m = gene + "@*_fs"
                    garc[(m, drug)] = confidenceToVal[values.get('who_confidence')]
                
                else:
                    print("Skipping: ", gene, mutation, info)

    #Checking for aminoglycosides tuple and expanding as required
    toRemove = []
    toAdd = {}
    for mutation, drug in garc.keys():
        if isinstance(drug, tuple):
            for d in drug:
                toAdd[(mutation, d)] = garc[(mutation, drug)]
            toRemove.append((mutation, drug))
    for r in toRemove:
        del garc[r]
    garc = {**garc, **toAdd}
    return garc



def proteinMutation(gene: str, mutation: str) -> str:
    '''Convert a protein mutation from iHGVS to GARC

    Args:
        gene (str): Gene name
        mutation (str): Mutation in iHGVS

    Returns:
        str: Mutation in GARC
    '''
    #We only have to worry about SNPs here
    #But we have to conver the AAs to be single letter
    snp = re.compile(r"""
                    p\. #Leading type
                    ([A-Za-z*]+) #Ref
                    ([0-9]+) #Pos
                    ([A-Za-z*]+) #Alt
                    """, re.VERBOSE)
    if snp.fullmatch(mutation):
        ref, pos, alt = snp.fullmatch(mutation).groups()
        #Convert AAs (dict adapted from tbdb's `parse_db.py`)
        aaDict = {
            "Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C",
            "Gln":"Q","Glu":"E","Gly":"G","His":"H","Ile":"I",
            "Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P",
            "Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V",
            "*":"!"
            }
        return gene + "@" + aaDict[ref] + pos + aaDict[alt]
    
    assert False, f"Nothing found for {gene} @ {mutation}"

def nucleotideVariants(gene: str, mutation: str) -> str:
    '''Convert nucelotide variants from iHGVS to GARC

    Args:
        gene (str): Gene name
        mutation (str): Mutation in iHGVS

    Returns:
        str: Mutation in GARC
    '''
    #First lets look at nucleotide SNPs
    snp = re.compile(r"""
                    [nc]\. #Leading type
                    (-?[0-9]+) #Gene position
                    ([AGCT])>([AGCT]) #Ref>Alt
                    """, re.VERBOSE)
    if snp.fullmatch(mutation):
        pos, ref, alt = snp.fullmatch(mutation).groups()
        return gene + "@" + ref.lower() + pos + alt.lower()
    
    #Check for ins
    ins = re.compile(r"""
                    [nc]\. #Leading type
                    (-?[0-9]+)_-?[0-9]+ #Pos is always 1 apart so ignore second
                    ins #Ins
                    ([ACGT]+) #Bases
                    """, re.VERBOSE)
    if ins.fullmatch(mutation):
        pos, ins = ins.fullmatch(mutation).groups()
        return gene + "@" + pos + "_ins_" + ins.lower()
    
    #Check for del
    #2 possible del forms:
    # `<pos>del<base>` or `<start>_<end>del<bases>`
    del_ = re.compile(r"""
                    [nc]\. #Leading type
                    (-?[0-9]+) #Start pos
                    _?(-?[0-9]+)? #Maybe end pos
                    del #Del
                    ([AGCT]+)? #Maybe bases deleted
                    """, re.VERBOSE)
    if del_.fullmatch(mutation):
        pos, end, bases = del_.fullmatch(mutation).groups()
        if bases is None:
            #No deleted bases given, so imply from range
            if end:
                #End value given, so find number of bases deleted
                bases = str(int(end) - int(pos) + 1)
            else:
                #No end value so just 1 base
                bases = "1"
        return gene + "@" + pos + "_del_" + bases.lower()

    #Check for duplications
    #2 possible dup forms:
    # `<pos>dup<base>` or `<start>_<end>dup<bases>`
    dup = re.compile(r"""
                    [nc]\. #Leading type
                    (-?[0-9]+) #Start pos
                    _?(-?[0-9]+)? #Maybe end pos
                    dup #Dup
                    ([AGCT]+) #Bases duplicated
                    """, re.VERBOSE)
    if dup.fullmatch(mutation):
        pos, end, bases = dup.fullmatch(mutation).groups()
        if end:
            #We want the ins after the end of the sequence, so prioritise this if existing
            pos = end
        return gene + "@" + pos + "_ins_" + bases.lower()

    assert False, f"Nothing found for: {gene} @ {mutation}"        

if __name__ == "__main__":
    resistant = pd.read_csv("tbdb/tbdb.csv")
    other = pd.read_csv("tbdb/tbdb.other_annotations.csv")
    # reference = gumpy.Genome("H37rV_v3.gbk", show_progress_bar=True)
    # pickle.dump(reference, open("reference.pkl", "wb"))
    reference = pickle.load(open("reference.pkl", "rb"))

    resistant = parseResistant(resistant, reference)

    print()
    print("******************")
    print()

    other = parseOther(other, reference)

    today = date.today()
    with open(f"tbdb-{today}.GARC.csv", "w") as f:
        f.write("GENBANK_REFERENCE,CATALOGUE_NAME,CATALOGUE_VERSION,CATALOGUE_GRAMMAR,PREDICTION_VALUES,DRUG,MUTATION,PREDICTION,SOURCE,EVIDENCE,OTHER\n")
        common = f"NC_000962.3,tbdb-{today},1.0,GARC1,RUS,"

        #All of these are R
        for mutation, drug in resistant.keys():
            f.write(common+drug+","+mutation+",R,{},{},{}\n")

        #Add the others
        for mutation, drug in other.keys():
            if (mutation, drug) not in resistant.keys():
                f.write(common+drug+","+mutation+","+other[(mutation, drug)]+",{},{},{}\n")   
        
