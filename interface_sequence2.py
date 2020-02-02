from Bio.PDB.PDBParser import PDBParser
import os
import numpy as np

path = "C:\\Users\\kitty\\Documents\\protein2.1"

PDBList = []
FileList = os.listdir(path)

for i in FileList:                           # 循环读取路径下的文件并筛选输出
    if os.path.splitext(i)[1] == ".pdb":   # 筛选pdb文件
        PDBList.append(i)


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
    string_a = str(ratio1a) + ',' + str(ratio2a) + ',' + str(ratio3a) + ',' + str(ratio4a) + ',' + '1'

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
    string_b = str(ratio1b) + ',' + str(ratio2b) + ',' + str(ratio3b) + ',' + str(ratio4b) + ',' + '1'
    return string_b


def residue_distance(residue_a, residue_b):
    global residue_list_a
    global residue_list_b
    if (residue_a["CA"] - residue_b["CA"]) < 9.5:
        residue_list_a.append(residue_a["CA"].get_full_id()[3][1])
        residue_list_b.append(residue_b["CA"].get_full_id()[3][1])


for file in PDBList:
    try:
        parser = PDBParser(PERMISSIVE=1)

        # read structure
        structure_id = os.path.splitext(file)[0]
        filename = file
        structure = parser.get_structure(structure_id, filename)

        model = structure[0]
        chain_A = model["A"]
        chain_B = model["B"]

        # function for calculating residue_distance
        residue_list_a = []
        residue_list_b = []

        for residue1 in chain_A.get_residues():
            for residue2 in chain_B.get_residues():
                try:
                    residue_distance(residue1, residue2)
                except Exception:
                    pass
                continue

        sequence_a = list(set(residue_list_a))
        sequence_b = list(set(residue_list_b))

        string_aa = feature_extraction_sequence2(sequence_a)
        string_bb = feature_extraction_sequence_b2(sequence_b)
        with open("result_interface_sequence2.txt", "a") as f:
            f.write(string_aa + "\n")
            f.write(string_bb + "\n")
        string1 = structure_id + ',' + str(len(sequence_a)) + ',' + str(
            float('%.3f' % (len(sequence_a) / len(chain_A)))) + ',' + str(len(sequence_b)) + ',' + str(
            float('%.3f' % (len(sequence_b) / len(chain_B))))
        with open("sequence2_test_interface03.txt", "a") as f:
            f.write(string1 + "\n")

    except Exception:
        pass
    continue

