import freesasa
from Bio.PDB.PDBParser import PDBParser
import os

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


def feature_extraction(surface_list_a1, surface_list_b1):
    global chain_A
    global chain_B

    surface_list_name_a = []
    surface_list_name_b = []
    list_a1 = []
    list_b1 = []
    for id_a in surface_list_a1:
        surface_list_name_a.append(chain_A[id_a].get_resname())

    for id_b in surface_list_b1:
        surface_list_name_b.append(chain_B[id_b].get_resname())

    # conjoint triad method
    for res in surface_list_name_a:
        if (res == 'ALA') or (res == 'GLY') or (res == 'VAL'):
            list_a1.append('type0')
        elif (res == 'ILE') or (res == 'LEU') or (res == 'PHE') or (res == 'PRO'):
            list_a1.append('type1')
        elif (res == 'TYR') or (res == 'MET') or (res == 'THR') or (res == 'SER'):
            list_a1.append('type2')
        elif (res == 'HIS') or (res == 'ASN') or (res == 'GLN') or (res == 'TRP'):
            list_a1.append('type3')
        elif (res == 'ARG') or (res == 'LYS'):
            list_a1.append('type4')
        elif (res == 'ASP') or (res == 'GLU'):
            list_a1.append('type5')
        elif res == 'CYS':
            list_a1.append('type6')

    # conjoint triad method
    for res in surface_list_name_b:
        if (res == 'ALA') or (res == 'GLY') or (res == 'VAL'):
            list_b1.append('type0')
        elif (res == 'ILE') or (res == 'LEU') or (res == 'PHE') or (res == 'PRO'):
            list_b1.append('type1')
        elif (res == 'TYR') or (res == 'MET') or (res == 'THR') or (res == 'SER'):
            list_b1.append('type2')
        elif (res == 'HIS') or (res == 'ASN') or (res == 'GLN') or (res == 'TRP'):
            list_b1.append('type3')
        elif (res == 'ARG') or (res == 'LYS'):
            list_b1.append('type4')
        elif (res == 'ASP') or (res == 'GLU'):
            list_b1.append('type5')
        elif res == 'CYS':
            list_b1.append('type6')

    return list_a1, list_b1


def file_output(list_a1, list_b1):
    sum2 = len(list_a1)
    type0a = float('%.6f' % ((list_a1.count('type0')) / sum2))
    type1a = float('%.6f' % ((list_a1.count('type1')) / sum2))
    type2a = float('%.6f' % ((list_a1.count('type2')) / sum2))
    type3a = float('%.6f' % ((list_a1.count('type3')) / sum2))
    type4a = float('%.6f' % ((list_a1.count('type4')) / sum2))
    type5a = float('%.6f' % ((list_a1.count('type5')) / sum2))
    type6a = float('%.6f' % ((list_a1.count('type6')) / sum2))

    sum1 = len(list_b1)
    type0 = float('%.6f' % ((list_b1.count('type0')) / sum1))
    type1 = float('%.6f' % ((list_b1.count('type1')) / sum1))
    type2 = float('%.6f' % ((list_b1.count('type2')) / sum1))
    type3 = float('%.6f' % ((list_b1.count('type3')) / sum1))
    type4 = float('%.6f' % ((list_b1.count('type4')) / sum1))
    type5 = float('%.6f' % ((list_b1.count('type5')) / sum1))
    type6 = float('%.6f' % ((list_b1.count('type6')) / sum1))

    string = str(type0) + ',' + str(type1) + ',' + str(type2) + ',' + str(type3) + ',' + str(type4) + ',' + str(
        type5) + ',' + str(type6) + ',' + '0'
    stringA = str(type0a) + ',' + str(type1a) + ',' + str(type2a) + ',' + str(type3a) + ',' + str(type4a) + ',' + str(
        type5a) + ',' + str(type6a) + ',' + '0'
    #print(string)
    with open("result_surface01.txt", "a") as f:
        f.write(string + "\n")
        f.write(stringA + "\n")


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
        list_a, list_b = feature_extraction(surface_list_a, surface_list_b)
        file_output(list_a, list_b)
        string1 = structure_id + ',' + str(len(surface_list_a)) + ',' + str(
            float('%.3f' % (len(surface_list_a) / len(chain_A)))) + ',' + str(len(surface_list_b)) + ',' + \
                  str(float('%.3f' % (len(surface_list_b) / len(chain_B))))

        with open("surface_test01.txt", "a") as f:
            f.write(string1 + "\n")
    except Exception:
        pass
    continue
