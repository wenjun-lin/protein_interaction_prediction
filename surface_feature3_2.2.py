import freesasa
from Bio.PDB.PDBParser import PDBParser
import os
import numpy as np

path = "C:\\Users\\kitty\\Documents\\protein2.0"

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


def feature_extraction_three(sequence_a2):
    sum_x = 0.0
    sum_y = 0.0
    sum_z = 0.0
    global chain_A
    for res in sequence_a2:
        sum_x = sum_x + chain_A[res]["CA"].get_coord()[0]
        sum_y = sum_y + chain_A[res]["CA"].get_coord()[1]
        sum_z = sum_z + chain_A[res]["CA"].get_coord()[2]

    x1 = sum_x / len(sequence_a2)
    y1 = sum_y / len(sequence_a2)
    z1 = sum_z / len(sequence_a2)
    vector_center = np.array([x1, y1, z1])

    distance_list_a = []
    for res_id in sequence_a2:
        distance_a = np.linalg.norm(chain_A[res_id]["CA"].get_vector() - vector_center)
        distance_list_a.append(distance_a)
    dis = np.array(distance_list_a)
    distance_a_max = dis.max()

    distance_a_list1 = []
    distance_a_list2 = []
    distance_a_list3 = []
    for res_id1 in sequence_a2:
        distance_a1 = np.linalg.norm(chain_A[res_id1]["CA"].get_vector() - vector_center)
        if distance_a1 > 0 and distance_a1 <= distance_a_max / 3:
            distance_a_list1.append(res_id1)
        elif distance_a1 > distance_a_max / 3 and distance_a1 <= (2 * distance_a_max) / 3:
            distance_a_list2.append(res_id1)
        else:
            distance_a_list3.append(res_id1)

    return distance_a_list1,distance_a_list2,distance_a_list3


def feature_extraction_three_b(sequence_b2):
    sum_xb = 0.0
    sum_yb = 0.0
    sum_zb = 0.0
    global chain_B
    for res in sequence_b2:
        sum_xb = sum_xb + chain_B[res]["CA"].get_coord()[0]
        sum_yb = sum_yb + chain_B[res]["CA"].get_coord()[1]
        sum_zb = sum_zb + chain_B[res]["CA"].get_coord()[2]

    x1 = sum_xb / len(sequence_b2)
    y1 = sum_yb / len(sequence_b2)
    z1 = sum_zb / len(sequence_b2)
    vector_center = np.array([x1, y1, z1])

    distance_list_b = []
    for res_id in sequence_b2:
        distance_b = np.linalg.norm(chain_B[res_id]["CA"].get_vector() - vector_center)
        distance_list_b.append(distance_b)
    dis = np.array(distance_list_b)
    distance_b_max = dis.max()

    distance_b_list1 = []
    distance_b_list2 = []
    distance_b_list3 = []
    for res_id1 in sequence_b2:
        distance_b1 = np.linalg.norm(chain_B[res_id1]["CA"].get_vector() - vector_center)
        if distance_b1 > 0 and distance_b1 <= distance_b_max / 3:
            distance_b_list1.append(res_id1)
        elif distance_b1 > distance_b_max / 3 and distance_b1 <= (2 * distance_b_max) / 3:
            distance_b_list2.append(res_id1)
        else:
            distance_b_list3.append(res_id1)

    return distance_b_list1, distance_b_list2, distance_b_list3


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
        dis_a1, dis_a2, dis_a3 = feature_extraction_three(surface_list_a)
        dis_b1, dis_b2, dis_b3 = feature_extraction_three_b(surface_list_b)

        string1a, string1b = list2result(dis_a1, dis_b1)
        string2a, string2b = list2result(dis_a1, dis_b1)
        string3a, string3b = list2result(dis_a1, dis_b1)
        string_aa = string1a + ',' + string2a + ',' + string3a+','+'0'
        string_bb = string1b + ',' + string2b + ',' + string3b+','+'0'

        with open("result_surface02.txt", "a") as f:
            f.write(string_aa + "\n")
            f.write(string_bb + "\n")

        string1 = structure_id + ',' + str(len(surface_list_a)) + ',' + str(
            float('%.3f' % (len(surface_list_a) / len(chain_A)))) + ',' + str(len(surface_list_b)) + ',' + \
                  str(float('%.3f' % (len(surface_list_b) / len(chain_B))))

        with open("surface_test02.txt", "a") as f:
            f.write(string1 + "\n")
    except Exception:
        pass
    continue
