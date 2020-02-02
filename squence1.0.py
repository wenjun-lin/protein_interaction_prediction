from Bio.PDB.PDBParser import PDBParser
import numpy as np


parser = PDBParser(PERMISSIVE=1)

# read structure
structure_id = "1a5f"
filename = "1a5f.pdb"
structure = parser.get_structure(structure_id, filename)

model = structure[0]
chain_A = model["L"]
chain_B = model["H"]

# function for calculating residue_distance
residue_list_a = []
residue_list_b = []
count_list_a = []
count_list_b = []
list_a = []
list_b = []


def residue_distance(residue_a, residue_b):
    global residue_list_a
    global residue_list_b
    if(residue_a["CA"] - residue_b["CA"]) < 9.5:
        residue_list_a.append(residue_a["CA"].get_full_id()[3][1])
        residue_list_b.append(residue_b["CA"].get_full_id()[3][1])


for residue1 in chain_A.get_residues():
    for residue2 in chain_B.get_residues():
        try:
            residue_distance(residue1, residue2)
        except Exception:
            pass
        continue

sequence_a = list(set(residue_list_a))
sequence_b = list(set(residue_list_b))


def feature_extraction_sequence(sequence_a2):

    global chain_A

    sequence_a_list1 = []
    sequence_a_list2 = []
    sequence_a_list3 = []
    sequence_a_list4 = []
    for res_id1 in sequence_a2:
        ratio = res_id1/len(chain_A)
        if ratio >=0 and ratio < 0.25:
            sequence_a_list1.append(res_id1)
        elif ratio > 0.25 and ratio <= 0.5:
            sequence_a_list2.append(res_id1)
        elif ratio > 0.5 and ratio <=0.75:
            sequence_a_list3.append(res_id1)
        elif ratio > 0.75 and ratio <= 1:
            sequence_a_list4.append(res_id1)
        else:
            continue

    return sequence_a_list1,sequence_a_list2,sequence_a_list3,sequence_a_list4


def feature_extraction_sequence_b(sequence_b2):

    global chain_B

    sequence_b_list1 = []
    sequence_b_list2 = []
    sequence_b_list3 = []
    sequence_b_list4 = []
    for res_id1 in sequence_b2:
        ratio = res_id1/len(chain_B)
        if ratio >=0 and ratio < 0.25:
            sequence_b_list1.append(res_id1)
        elif ratio > 0.25 and ratio <= 0.5:
            sequence_b_list2.append(res_id1)
        elif ratio > 0.5 and ratio <=0.75:
            sequence_b_list3.append(res_id1)
        elif ratio > 0.75 and ratio <= 1:
            sequence_b_list4.append(res_id1)
        else:
            continue

    return sequence_b_list1,sequence_b_list2,sequence_b_list3,sequence_b_list4


def list2result(sequence_aid,sequence_bid):
    global chain_A
    global chain_B
    count_list_a1 = []
    count_list_b1 = []
    list_a1 = []
    list_b1 = []
    for id_a in sequence_aid:
        count_list_a1.append(chain_A[id_a].get_resname())
    for id_b in sequence_bid:
        count_list_b1.append(chain_B[id_b].get_resname())
    for res1 in count_list_a1:
        if (res1 == 'ALA') or (res1 == 'GLY') or (res1 == 'VAL'):
            list_a1.append('type0')
        elif (res1 == 'ILE') or (res1 == 'LEU') or (res1 == 'PHE') or (res1 == 'PRO'):
            list_a1.append('type1')
        elif (res1 == 'TYR') or (res1 == 'MET') or (res1 == 'THR') or (res1 == 'SER'):
            list_a1.append('type2')
        elif (res1 == 'HIS') or (res1 == 'ASN') or (res1 == 'GLN') or (res1 == 'TRP'):
            list_a1.append('type3')
        elif (res1 == 'ARG') or (res1 == 'LYS'):
            list_a1.append('type4')
        elif (res1 == 'ASP') or (res1 == 'GLU'):
            list_a1.append('type5')
        elif res1 == 'CYS':
            list_a1.append('type6')
    sum2 = len(list_a1)
    type0a = float('%.6f' % ((list_a1.count('type0')) / sum2))
    type1a = float('%.6f' % ((list_a1.count('type1')) / sum2))
    type2a = float('%.6f' % ((list_a1.count('type2')) / sum2))
    type3a = float('%.6f' % ((list_a1.count('type3')) / sum2))
    type4a = float('%.6f' % ((list_a1.count('type4')) / sum2))
    type5a = float('%.6f' % ((list_a1.count('type5')) / sum2))
    type6a = float('%.6f' % ((list_a1.count('type6')) / sum2))
    string_a = str(type0a) + ',' + str(type1a) + ',' + str(type2a) + ',' + str(type3a) + ',' + str(type4a) + ',' + str(
        type5a) + ',' + str(type6a)

    for res2 in count_list_b1:
        if (res2 == 'ALA') or (res2 == 'GLY') or (res2 == 'VAL'):
            list_b1.append('type0')
        elif (res2 == 'ILE') or (res2 == 'LEU') or (res2 == 'PHE') or (res2 == 'PRO'):
            list_b1.append('type1')
        elif (res2 == 'TYR') or (res2 == 'MET') or (res2 == 'THR') or (res2 == 'SER'):
            list_b1.append('type2')
        elif (res2 == 'HIS') or (res2 == 'ASN') or (res2 == 'GLN') or (res2 == 'TRP'):
            list_b1.append('type3')
        elif (res2 == 'ARG') or (res2 == 'LYS'):
            list_b1.append('type4')
        elif (res2 == 'ASP') or (res2 == 'GLU'):
            list_b1.append('type5')
        elif res2 == 'CYS':
            list_b1.append('type6')
    sum1 = len(list_b1)
    type0b = float('%.6f' % ((list_b1.count('type0')) / sum1))
    type1b = float('%.6f' % ((list_b1.count('type1')) / sum1))
    type2b = float('%.6f' % ((list_b1.count('type2')) / sum1))
    type3b = float('%.6f' % ((list_b1.count('type3')) / sum1))
    type4b = float('%.6f' % ((list_b1.count('type4')) / sum1))
    type5b = float('%.6f' % ((list_b1.count('type5')) / sum1))
    type6b = float('%.6f' % ((list_b1.count('type6')) / sum1))
    string_b = str(type0b) + ',' + str(type1b) + ',' + str(type2b) + ',' + str(type3b) + ',' + str(type4b) + ',' + str(
        type5b) + ',' + str(type6b)

    return string_a, string_b


seq_a1, seq_a2, seq_a3, seq_a4 = feature_extraction_sequence(sequence_a)
seq_b1, seq_b2, seq_b3, seq_b4 = feature_extraction_sequence_b(sequence_b)


string1a,string1b = list2result(seq_a1,seq_b1)
string2a,string2b = list2result(seq_a2,seq_b2)
string3a,string3b = list2result(seq_a3,seq_b3)
string4a,string4b = list2result(seq_a4,seq_b4)
string_aa = string1a + ',' + string2a + ',' + string3a + ',' + string4a
string_bb = string1b + ',' + string2b + ',' + string3b + ',' + string4b
print(string_aa)
print(string_bb)



