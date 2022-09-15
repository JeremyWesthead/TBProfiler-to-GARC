'''Parse the TBProfiler DB to GARC for use with `piezo`
'''
import copy
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

def addExtras(reference: gumpy.Genome) -> None:
    '''(Adapted from the code written to parse WHO cat to GARC)
    Once the catalogue has been parsed correctly, there will be some mutations which also lie within other genes
    This finds them and adds them to the catalogue. Specifically, this checks for promoter SNPs which could be attributed
    to other genes, especially in cases where the promoter position is beyond the arbitrary internal limits of gumpy.
    Args:
        reference (gumpy.Genome): Reference genome
    '''
    today = date.today()
    catalogue = pd.read_csv(f"tbdb-{today}.GARC.csv")
    toAdd = {column: [] for column in catalogue}
    
    #Track the new resistance genes this introduces to add default rules
    newGenes = set()
    previousGenes = set([mutation.split("@")[0] for mutation in catalogue['MUTATION']])
    for _, row in catalogue.iterrows():
        mut = row['MUTATION']
        #Check for promoter
        if "-" not in mut:
            continue
        #Check for default rules/multi for skipping
        if "*" in mut or "?" in mut or "&" in mut or "indel" in mut:
            continue
        promoter = re.compile(r"""
                            ([a-zA-Z0-9_]+)@ #Leading gene name
                            ([a-z])(-[0-9]+)([a-z])
                            """, re.VERBOSE)
        if promoter.fullmatch(mut):
            gene, ref, pos, alt = promoter.fullmatch(mut).groups()
            pos = int(pos)
            sample = copy.deepcopy(reference)
            print(gene, pos, ref, alt)
            
            #Place the mutation within the genome based on the gene coordinates
            #Then regardless of what gene it started in, we can pull out others
            if reference.genes[gene]['reverse_complement']:
                #Revcomp genes' promoters will be past the `gene end`
                geneEnd = reference.genes[gene]['end']
                ref_ = ''.join(gumpy.Gene._complement(ref))
                alt_ = ''.join(gumpy.Gene._complement(alt))
                pos_ = geneEnd-pos-1
                assert reference.nucleotide_sequence[reference.nucleotide_index == geneEnd-pos-1] == ref_, "Ref does not match the genome..."
                sample.nucleotide_sequence[reference.nucleotide_index == geneEnd-pos-1] = alt_
            else:
                geneStart = reference.genes[gene]['start']
                pos_ = geneStart+pos
                assert reference.nucleotide_sequence[reference.nucleotide_index == geneStart+pos] == ref, "Ref does not match the genome..."
                sample.nucleotide_sequence[reference.nucleotide_index == geneStart+pos] = alt
            
            #The only mutations between ref and sample are this SNP
            #So pull out all available mutations (ignoring the original gene)
            mutations = []
            #Get genes at this position
            possible = [reference.stacked_gene_name[i][pos_] for i in range(len(reference.stacked_gene_name)) if reference.stacked_gene_name[i][pos_] != '']
            for g in possible:
                if g == gene:
                    continue
                if row['PREDICTION'] == "R" and g not in previousGenes:
                    newGenes.add((g, row['DRUG']))
                diff = reference.build_gene(g) - sample.build_gene(g)
                m = diff.mutations
                if m:
                    for mut_ in m:
                        mutations.append(g+"@"+mut_)
            
            #Make them neat catalouge rows to add
            for m in mutations:
                for col in catalogue:
                    if col == "MUTATION":
                        toAdd[col].append(m)
                    else:
                        toAdd[col].append(row[col])
    for gene, drug in newGenes:
        #These are new resistance genes, so add default rules as appropriate
        defaults = [
            (gene+"@*?", 'U'), (gene+"@-*?", 'U'),
            (gene+"@*_indel", "U"), (gene+"@-*_indel", 'U')
            ]
        if reference.genes[gene]['codes_protein']:
            defaults.append((gene+"@*=","S"))
        for g, predict in defaults:
            for col in catalogue:
                if col == "MUTATION":
                    toAdd[col].append(g)
                elif col == "DRUG":
                    toAdd[col].append(drug)
                elif col == "PREDICTION":
                    toAdd[col].append(predict)
                elif col in ["SOURCE", "EVIDENCE", "OTHER"]:
                    toAdd[col].append("{}")
                else:
                    #Others should be constant
                    toAdd[col].append(toAdd[col][-1])
    #Convert toAdd to dataframe and concat with catalogue
    toAdd = pd.DataFrame(toAdd)
    catalogue = pd.concat([catalogue, toAdd])
    catalogue.to_csv(f"tbdb-{today}.GARC.csv", index=False)
    

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

    #Set of (gene, drug)
    resistanceGenes = set()
    for mutation, drug in resistant.keys():
        resistanceGenes.add((mutation.split("@")[0], drug))
    for mutation, drug in other.keys():
        if other[(mutation, drug)] == 'R':
            resistanceGenes.add((mutation.split("@")[0], drug))

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
        
        #Add default rules for resistance genes
        for gene, drug in resistanceGenes:
            f.write(common + drug + "," + gene + "@*?,U,{},{},{}\n")
            f.write(common + drug + "," + gene + "@-*?,U,{},{},{}\n")
            f.write(common + drug + "," + gene + "@*_indel,U,{},{},{}\n")
            f.write(common + drug + "," + gene + "@-*_indel,U,{},{},{}\n")
            if reference.genes[gene]['codes_protein']:
                f.write(common + drug + "," + gene + "@*=,S,{},{},{}\n")
    
    #Add the alternate forms of some mutations which would otherwise be missed by gumpy
    addExtras(reference)
        
