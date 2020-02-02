import freesasa
from Bio.PDB.PDBParser import PDBParser
import os
import numpy as np

path = "C:\\Users\\kitty\\Documents\\protein2.1"

PDBList = []
FileList = os.listdir(path)

for i in FileList:                           # 循环读取路径下的文件并筛选输出
    if os.path.splitext(i)[1] == ".pdb":   # 筛选pdb文件
        PDBList.append(i)


def surface_list(file1):

    maximum_area = {'ALA': 120.56, 'CYS': 143.79, 'ASP': 157.04, 'GLU': 188.42, 'PHE': 227.46, 'GLY': 89.41,
                    'HIS': 200.14, 'ILE': 96.42, 'LYS': 213.74, 'LEU': 206.32, 'MET': 216.63, 'ASN': 149.85,
                    'PRO': 155.07, 'GLN': 186.83, 'ARG': 229.51, 'SER': 128.27, 'THR': 138.58, 'VAL': 169.82,
                    'TRP': 269.35, 'TYR': 241.54}

    global chain_A
    global chain_B

    surface_list_a1 = []
    surface_list_b1 = []

    structure = freesasa.Structure(file1)
    result = freesasa.calc(structure)

    for residue1 in chain_A.get_residues():
        try:
            res_id = residue1["CA"].get_full_id()[3][1]
            select_word = str(res_id) + ", " + "chain H and resi " + str(res_id)
            selections = freesasa.selectArea((select_word,), structure, result)
            for key in selections:
                if float('%.3f' % (
                        selections[key] / maximum_area[chain_A[residue1.get_full_id()[3][1]].get_resname()])) > 0.05:
                    surface_list_a1.append(res_id)
        except Exception:
            pass
        continue

    for residue2 in chain_B.get_residues():
        try:
            res_id = residue2["CA"].get_full_id()[3][1]
            select_word = str(res_id) + ", " + "chain L and resi " + str(res_id)
            selections = freesasa.selectArea((select_word,), structure, result)
            for key in selections:
                if float('%.3f' % (
                        selections[key] / maximum_area[chain_B[residue2.get_full_id()[3][1]].get_resname()])) > 0.05:
                    surface_list_b1.append(res_id)
        except Exception:
            pass
        continue

    return surface_list_a1, surface_list_b1


def feature_extraction_sequence2(sequence_a2):

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
    a1 = len(sequence_a_list1)
    a2 = len(sequence_a_list2)
    a3 = len(sequence_a_list3)
    a4 = len(sequence_a_list4)
    num2 = a1 + a2 + a3 + a4
    ratio1a = float('%.6f' % (a1 / num2))
    ratio2a = float('%.6f' % (a2 / num2))
    ratio3a = float('%.6f' % (a3 / num2))
    ratio4a = float('%.6f' % (a4 / num2))
    string_a = str(ratio1a) + ',' + str(ratio2a) + ',' + str(ratio3a) + ',' + str(ratio4a) + ',' + '0'

    return string_a


def feature_extraction_sequence_b2(sequence_b2):

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
    b1 = len(sequence_b_list1)
    b2 = len(sequence_b_list2)
    b3 = len(sequence_b_list3)
    b4 = len(sequence_b_list4)
    num1 = b1 + b2 + b3 + b4
    ratio1b = float('%.6f' % (b1/num1))
    ratio2b = float('%.6f' % (b2/num1))
    ratio3b = float('%.6f' % (b3/num1))
    ratio4b = float('%.6f' % (b4/num1))
    string_b = str(ratio1b) + ',' + str(ratio2b) + ',' + str(ratio3b) + ',' + str(ratio4b) + ',' + '0'
    return string_b


for file in PDBList:
    try:
        parser = PDBParser(PERMISSIVE=1)
        structure_id = os.path.splitext(file)[0]
        filename = file
        structure1 = parser.get_structure(structure_id, filename)

        model = structure1[0]
        chain_A = model["H"]
        chain_B = model["L"]
        surface_list_a, surface_list_b = surface_list(file)

        string_aa = feature_extraction_sequence2(surface_list_a)
        string_bb = feature_extraction_sequence_b2(surface_list_b)

        with open("result_surface_sequence2.txt", "a") as f:
            f.write(string_aa + "\n")
            f.write(string_bb + "\n")

        string1 = structure_id + ',' + str(len(surface_list_a)) + ',' + str(
            float('%.3f' % (len(surface_list_a) / len(chain_A)))) + ',' + str(len(surface_list_b)) + ',' + \
                  str(float('%.3f' % (len(surface_list_b) / len(chain_B))))

        with open("sequence2_surface_test03.txt", "a") as f:
            f.write(string1 + "\n")
    except Exception:
        pass
    continue
