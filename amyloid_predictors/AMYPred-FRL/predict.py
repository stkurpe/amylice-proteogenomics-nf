import os
import sys
import re
import math
import pickle
import numpy as np
import pandas as pd
from collections import Counter
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from xgboost import XGBClassifier
from tabulate import tabulate

# =============================================================================
# 1. ФУНКЦИИ ИЗВЛЕЧЕНИЯ ПРИЗНАКОВ
# =============================================================================

def read_protein_sequences(file):
    if not os.path.exists(file):
        print('Error: file %s does not exist.' % file)
        sys.exit(1)
    with open(file) as f:
        records = f.read()
    if re.search('>', records) == None:
        print('Error: the input file %s seems not in FASTA format!' % file)
        sys.exit(1)
    records = records.split('>')[1:]
    fasta_sequences = []
    
    # Минимальная длина последовательности
    MIN_LENGTH = 30 
    
    for fasta in records:
        array = fasta.split('\n')
        header_line = array[0].strip()
        if not header_line: continue 
        seq_id = header_line.split()[0] 
        
        sequence = re.sub('[^ACDEFGHIKLMNPQRSTVWY-]', '-', ''.join(array[1:]).upper())
        
        # Фильтр длины
        if len(sequence) < MIN_LENGTH:
            print(f"Warning: Sequence '{seq_id}' skipped (Length {len(sequence)} < {MIN_LENGTH}). Too short for feature generation.")
            continue
            
        fasta_sequences.append([seq_id, sequence])
        
    if len(fasta_sequences) == 0:
        print("Error: No valid sequences left after filtering! All sequences were too short or file is empty.")
        sys.exit(1)
        
    return fasta_sequences

def AAC(fastas, **kw):
    AA = 'ACDEFGHIKLMNPQRSTVWY'
    encodings = []
    header = [i for i in AA]
    for i in fastas:
        sequence = re.sub('-', '', i[1])
        count = Counter(sequence)
        for key in count:
            count[key] = count[key]/len(sequence)
        code = [count[aa] for aa in AA]
        encodings.append(code)
    return np.array(encodings, dtype=float), header

def DPC(fastas, gap, **kw):
    AA = 'ACDEFGHIKLMNPQRSTVWY'
    encodings = []
    diPeptides = [aa1 + aa2 for aa1 in AA for aa2 in AA]
    header = [] + diPeptides
    AADict = {AA[i]: i for i in range(len(AA))}
    for i in fastas:
        sequence = re.sub('-', '', i[1])
        code = []
        tmpCode = [0] * 400
        for j in range(len(sequence) - 2 + 1 - gap):
            tmpCode[AADict[sequence[j]] * 20 + AADict[sequence[j+gap+1]]] += 1
        if sum(tmpCode) != 0:
            tmpCode = [i/sum(tmpCode) for i in tmpCode]
        code = code + tmpCode
        encodings.append(code)
    return np.array(encodings, dtype=float), header

def APAAC(fastas, lambdaValue=10, w=0.05, **kw):
    dataFile = 'data/PAAC.txt'
    if not os.path.exists(dataFile):
        print(f"Error: {dataFile} not found. Ensure 'data' folder is present.")
        sys.exit(1)
    with open(dataFile) as f:
        records = f.readlines()
    AA = ''.join(records[0].rstrip().split()[1:])
    AADict = {AA[i]: i for i in range(len(AA))}
    AAProperty = []
    AAPropertyNames = []
    for i in range(1, len(records) - 1):
        array = records[i].rstrip().split() if records[i].rstrip() != '' else None
        AAProperty.append([float(j) for j in array[1:]])
        AAPropertyNames.append(array[0])
    AAProperty1 = []
    for i in AAProperty:
        meanI = sum(i) / 20
        fenmu = math.sqrt(sum([(j - meanI) ** 2 for j in i]) / 20)
        AAProperty1.append([(j - meanI) / fenmu for j in i])
    encodings = []
    header = ['Pc1.' + i for i in AA]
    for j in range(1, lambdaValue + 1):
        for i in AAPropertyNames:
            header.append('Pc2.' + i + '.' + str(j))
    for i in fastas:
        sequence = re.sub('-', '', i[1])
        theta = []
        for n in range(1, lambdaValue + 1):
            for j in range(len(AAProperty1)):
                theta.append(sum([AAProperty1[j][AADict[sequence[k]]] * AAProperty1[j][AADict[sequence[k + n]]] for k in
                                  range(len(sequence) - n)]) / (len(sequence) - n))
        myDict = {aa: sequence.count(aa) for aa in AA}
        code = [myDict[aa] / (1 + w * sum(theta)) for aa in AA]
        code = code + [w * value / (1 + w * sum(theta)) for value in theta]
        encodings.append(code)
    return np.array(encodings, dtype=float), header

def Count(seq1, seq2):
    return sum([seq2.count(aa) for aa in seq1])

def CTDC(fastas, **kw):
    group1 = {'hydrophobicity_PRAM900101': 'RKEDQN', 'hydrophobicity_ARGP820101': 'QSTNGDE', 'hydrophobicity_ZIMJ680101': 'QNGSWTDERA', 'hydrophobicity_PONP930101': 'KPDESNQT', 'hydrophobicity_CASG920101': 'KDEQPSRNTG', 'hydrophobicity_ENGD860101': 'RDKENQHYP', 'hydrophobicity_FASG890101': 'KERSQD', 'normwaalsvolume': 'GASTPDC', 'polarity': 'LIFWCMVY', 'polarizability': 'GASDT', 'charge': 'KR', 'secondarystruct': 'EALMQKRH', 'solventaccess': 'ALFCGIVW'}
    group2 = {'hydrophobicity_PRAM900101': 'GASTPHY', 'hydrophobicity_ARGP820101': 'RAHCKMV', 'hydrophobicity_ZIMJ680101': 'HMCKV', 'hydrophobicity_PONP930101': 'GRHA', 'hydrophobicity_CASG920101': 'AHYMLV', 'hydrophobicity_ENGD860101': 'SGTAW', 'hydrophobicity_FASG890101': 'NTPG', 'normwaalsvolume': 'NVEQIL', 'polarity': 'PATGS', 'polarizability': 'CPNVEQIL', 'charge': 'ANCQGHILMFPSTWYV', 'secondarystruct': 'VIYCWFT', 'solventaccess': 'RKQEND'}
    group3 = {'hydrophobicity_PRAM900101': 'CLVIMFW', 'hydrophobicity_ARGP820101': 'LYPFIW', 'hydrophobicity_ZIMJ680101': 'LPFYI', 'hydrophobicity_PONP930101': 'YMFWLCVI', 'hydrophobicity_CASG920101': 'FIWC', 'hydrophobicity_ENGD860101': 'CVLIMF', 'hydrophobicity_FASG890101': 'AYHWVMFLIC', 'normwaalsvolume': 'MHKFRYW', 'polarity': 'HQRKNED', 'polarizability': 'KMHFRYW', 'charge': 'DE', 'secondarystruct': 'GNPSD', 'solventaccess': 'MSPTHY'}
    groups = [group1, group2, group3]
    property = tuple(group1.keys())
    encodings = []
    header = []
    for p in property:
        for g in range(1, len(groups) + 1):
            header.append(p + '.G' + str(g))
    for i in fastas:
        sequence = re.sub('-', '', i[1])
        code = []
        for p in property:
            c1 = Count(group1[p], sequence) / len(sequence)
            c2 = Count(group2[p], sequence) / len(sequence)
            c3 = 1 - c1 - c2
            code = code + [c1, c2, c3]
        encodings.append(code)
    return np.array(encodings, dtype=float), header

def Count2(aaSet, sequence):
    number = 0
    for aa in sequence:
        if aa in aaSet:
            number = number + 1
    cutoffNums = [1, math.floor(0.25 * number), math.floor(0.50 * number), math.floor(0.75 * number), number]
    cutoffNums = [i if i >=1 else 1 for i in cutoffNums]
    code = []
    for cutoff in cutoffNums:
        myCount = 0
        for i in range(len(sequence)):
            if sequence[i] in aaSet:
                myCount += 1
                if myCount == cutoff:
                    code.append((i + 1) / len(sequence) * 100)
                    break
        if myCount == 0:
            code.append(0)
    return code

def CTDD(fastas, **kw):
    group1 = {'hydrophobicity_PRAM900101': 'RKEDQN', 'hydrophobicity_ARGP820101': 'QSTNGDE', 'hydrophobicity_ZIMJ680101': 'QNGSWTDERA', 'hydrophobicity_PONP930101': 'KPDESNQT', 'hydrophobicity_CASG920101': 'KDEQPSRNTG', 'hydrophobicity_ENGD860101': 'RDKENQHYP', 'hydrophobicity_FASG890101': 'KERSQD', 'normwaalsvolume': 'GASTPDC', 'polarity': 'LIFWCMVY', 'polarizability': 'GASDT', 'charge': 'KR', 'secondarystruct': 'EALMQKRH', 'solventaccess': 'ALFCGIVW'}
    group2 = {'hydrophobicity_PRAM900101': 'GASTPHY', 'hydrophobicity_ARGP820101': 'RAHCKMV', 'hydrophobicity_ZIMJ680101': 'HMCKV', 'hydrophobicity_PONP930101': 'GRHA', 'hydrophobicity_CASG920101': 'AHYMLV', 'hydrophobicity_ENGD860101': 'SGTAW', 'hydrophobicity_FASG890101': 'NTPG', 'normwaalsvolume': 'NVEQIL', 'polarity': 'PATGS', 'polarizability': 'CPNVEQIL', 'charge': 'ANCQGHILMFPSTWYV', 'secondarystruct': 'VIYCWFT', 'solventaccess': 'RKQEND'}
    group3 = {'hydrophobicity_PRAM900101': 'CLVIMFW', 'hydrophobicity_ARGP820101': 'LYPFIW', 'hydrophobicity_ZIMJ680101': 'LPFYI', 'hydrophobicity_PONP930101': 'YMFWLCVI', 'hydrophobicity_CASG920101': 'FIWC', 'hydrophobicity_ENGD860101': 'CVLIMF', 'hydrophobicity_FASG890101': 'AYHWVMFLIC', 'normwaalsvolume': 'MHKFRYW', 'polarity': 'HQRKNED', 'polarizability': 'KMHFRYW', 'charge': 'DE', 'secondarystruct': 'GNPSD', 'solventaccess': 'MSPTHY'}
    property = tuple(group1.keys())
    encodings = []
    header = []
    for p in property:
        for g in ('1', '2', '3'):
            for d in ['0', '25', '50', '75', '100']:
                header.append(p + '.' + g + '.residue' + d)
    for i in fastas:
        sequence = re.sub('-', '', i[1])
        code = []
        for p in property:
            code = code + Count2(group1[p], sequence) + Count2(group2[p], sequence) + Count2(group3[p], sequence)
        encodings.append(code)
    return np.array(encodings, dtype=float), header

def CTDT(fastas, **kw):
    group1 = {'hydrophobicity_PRAM900101': 'RKEDQN', 'hydrophobicity_ARGP820101': 'QSTNGDE', 'hydrophobicity_ZIMJ680101': 'QNGSWTDERA', 'hydrophobicity_PONP930101': 'KPDESNQT', 'hydrophobicity_CASG920101': 'KDEQPSRNTG', 'hydrophobicity_ENGD860101': 'RDKENQHYP', 'hydrophobicity_FASG890101': 'KERSQD', 'normwaalsvolume': 'GASTPDC', 'polarity': 'LIFWCMVY', 'polarizability': 'GASDT', 'charge': 'KR', 'secondarystruct': 'EALMQKRH', 'solventaccess': 'ALFCGIVW'}
    group2 = {'hydrophobicity_PRAM900101': 'GASTPHY', 'hydrophobicity_ARGP820101': 'RAHCKMV', 'hydrophobicity_ZIMJ680101': 'HMCKV', 'hydrophobicity_PONP930101': 'GRHA', 'hydrophobicity_CASG920101': 'AHYMLV', 'hydrophobicity_ENGD860101': 'SGTAW', 'hydrophobicity_FASG890101': 'NTPG', 'normwaalsvolume': 'NVEQIL', 'polarity': 'PATGS', 'polarizability': 'CPNVEQIL', 'charge': 'ANCQGHILMFPSTWYV', 'secondarystruct': 'VIYCWFT', 'solventaccess': 'RKQEND'}
    group3 = {'hydrophobicity_PRAM900101': 'CLVIMFW', 'hydrophobicity_ARGP820101': 'LYPFIW', 'hydrophobicity_ZIMJ680101': 'LPFYI', 'hydrophobicity_PONP930101': 'YMFWLCVI', 'hydrophobicity_CASG920101': 'FIWC', 'hydrophobicity_ENGD860101': 'CVLIMF', 'hydrophobicity_FASG890101': 'AYHWVMFLIC', 'normwaalsvolume': 'MHKFRYW', 'polarity': 'HQRKNED', 'polarizability': 'KMHFRYW', 'charge': 'DE', 'secondarystruct': 'GNPSD', 'solventaccess': 'MSPTHY'}
    property = tuple(group1.keys())
    encodings = []
    header = []
    for p in property:
        for tr in ('Tr1221', 'Tr1331', 'Tr2332'):
            header.append(p + '.' + tr)
    for i in fastas:
        sequence = re.sub('-', '', i[1])
        code = []
        aaPair = [sequence[j:j + 2] for j in range(len(sequence) - 1)]
        for p in property:
            c1221, c1331, c2332 = 0, 0, 0
            for pair in aaPair:
                if (pair[0] in group1[p] and pair[1] in group2[p]) or (pair[0] in group2[p] and pair[1] in group1[p]):
                    c1221 += 1
                elif (pair[0] in group1[p] and pair[1] in group3[p]) or (pair[0] in group3[p] and pair[1] in group1[p]):
                    c1331 += 1
                elif (pair[0] in group2[p] and pair[1] in group3[p]) or (pair[0] in group3[p] and pair[1] in group2[p]):
                    c2332 += 1
            code = code + [c1221/len(aaPair), c1331/len(aaPair), c2332/len(aaPair)]
        encodings.append(code)
    return np.array(encodings, dtype=float), header

def CalculateKSCTriad(sequence, gap, features, AADict):
    res = []
    for g in range(gap + 1):
        myDict = {f: 0 for f in features}
        for i in range(len(sequence)):
            if i + g + 1 < len(sequence) and i + 2 * g + 2 < len(sequence):
                fea = AADict[sequence[i]] + '.' + AADict[sequence[i + g + 1]] + '.' + AADict[sequence[i + 2 * g + 2]]
                myDict[fea] += 1
        maxValue, minValue = max(myDict.values()), min(myDict.values())
        for f in features:
            res.append((myDict[f] - minValue) / maxValue)
    return res

def KSCTriad(fastas, gap=0, **kw):
    AAGroup = {'g1': 'AGV', 'g2': 'ILFP', 'g3': 'YMTS', 'g4': 'HNQW', 'g5': 'RK', 'g6': 'DE', 'g7': 'C'}
    myGroups = sorted(AAGroup.keys())
    AADict = {aa: g for g in myGroups for aa in AAGroup[g]}
    features = [f1 + '.' + f2 + '.' + f3 for f1 in myGroups for f2 in myGroups for f3 in myGroups]
    encodings = []
    header = [f + '.gap' + str(g) for g in range(gap + 1) for f in features]
    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [] # ИСПРАВЛЕНО: убрали добавление name
        
        # На случай, если основной фильтр не покрыл
        if len(sequence) < 2 * gap + 3:
            return np.zeros((1, len(header))), header
            
        code = code + CalculateKSCTriad(sequence, gap, features, AADict)
        encodings.append(code)
    return np.array(encodings, dtype=float), header

def CTriad(fastas, gap = 0, **kw):
    AAGroup = {'g1': 'AGV', 'g2': 'ILFP', 'g3': 'YMTS', 'g4': 'HNQW', 'g5': 'RK', 'g6': 'DE', 'g7': 'C'}
    myGroups = sorted(AAGroup.keys())
    AADict = {aa: g for g in myGroups for aa in AAGroup[g]}
    features = [f1 + '.'+ f2 + '.' + f3 for f1 in myGroups for f2 in myGroups for f3 in myGroups]
    encodings = []
    header = features
    for i in fastas:
        sequence = re.sub('-', '', i[1])
        code = []
        if len(sequence) < 3:
            continue
        code = code + CalculateKSCTriad(sequence, 0, features, AADict)
        encodings.append(code)
    return np.array(encodings, dtype=float), header

def DDE(fastas, **kw):
    AA = 'ACDEFGHIKLMNPQRSTVWY'
    myCodons = {'A': 4, 'C': 2, 'D': 2, 'E': 2, 'F': 2, 'G': 4, 'H': 2, 'I': 3, 'K': 2, 'L': 6, 'M': 1, 'N': 2, 'P': 4, 'Q': 2, 'R': 6, 'S': 6, 'T': 4, 'V': 4, 'W': 1, 'Y': 2}
    encodings = []
    diPeptides = ['DDE_'+aa1 + aa2 for aa1 in AA for aa2 in AA]
    header = diPeptides
    myTM = [(myCodons[p[0]] / 61) * (myCodons[p[1]] / 61) for p in [aa1+aa2 for aa1 in AA for aa2 in AA]]
    AADict = {AA[i]: i for i in range(len(AA))}
    for i in fastas:
        sequence = re.sub('-', '', i[1])
        tmpCode = [0] * 400
        for j in range(len(sequence) - 1):
            tmpCode[AADict[sequence[j]] * 20 + AADict[sequence[j+1]]] += 1
        if sum(tmpCode) != 0:
            tmpCode = [x/sum(tmpCode) for x in tmpCode]
        myTV = [x * (1-x) / (len(sequence) - 1) for x in myTM]
        
        final_code = []
        for j in range(len(tmpCode)):
            val = 0
            if myTV[j] > 1e-9:
                val = (tmpCode[j] - myTM[j]) / math.sqrt(myTV[j])
            final_code.append(val)
            
        encodings.append(final_code)
    return np.array(encodings, dtype=float), header

def GAAC(fastas, **kw):
    group = {'alphatic': 'GAVLMI', 'aromatic': 'FYW', 'postivecharge': 'KRH', 'negativecharge': 'DE', 'uncharge': 'STCPNQ'}
    groupKey = group.keys()
    encodings = []
    header = list(groupKey)
    for i in fastas:
        name, sequence= i[0], re.sub('-', '', i[1])
        code = [] # ИСПРАВЛЕНО: убрали добавление name
        count = Counter(sequence)
        myDict = {key: 0 for key in groupKey}
        for key in groupKey:
            for aa in group[key]:
                myDict[key] += count[aa]
        for key in groupKey:
            code.append(myDict[key]/len(sequence))
        encodings.append(code)
    return np.array(encodings, dtype=float), header

def Rvalue(aa1, aa2, AADict, Matrix):
    return sum([(Matrix[i][AADict[aa1]] - Matrix[i][AADict[aa2]]) ** 2 for i in range(len(Matrix))]) / len(Matrix)

def PAAC(fastas, lambdaValue=5, w=0.05, **kw):
    dataFile = 'data/PAAC.txt'
    with open(dataFile) as f:
        records = f.readlines()
    AA = ''.join(records[0].rstrip().split()[1:])
    AADict = {AA[i]: i for i in range(len(AA))}
    AAProperty = []
    for i in range(1, len(records)):
        array = records[i].rstrip().split() if records[i].rstrip() != '' else None
        AAProperty.append([float(j) for j in array[1:]])
    AAProperty1 = []
    for i in AAProperty:
        meanI = sum(i) / 20
        fenmu = math.sqrt(sum([(j - meanI) ** 2 for j in i]) / 20)
        AAProperty1.append([(j - meanI) / fenmu for j in i])
    encodings = []
    header = ['Xc1.' + aa for aa in AA] + ['Xc2.lambda' + str(n) for n in range(1, lambdaValue + 1)]
    for i in fastas:
        sequence = re.sub('-', '', i[1])
        theta = []
        for n in range(1, lambdaValue + 1):
            theta.append(sum([Rvalue(sequence[j], sequence[j + n], AADict, AAProperty1) for j in range(len(sequence) - n)]) / (len(sequence) - n))
        myDict = {aa: sequence.count(aa) for aa in AA}
        code = [myDict[aa] / (1 + w * sum(theta)) for aa in AA]
        code = code + [(w * j) / (1 + w * sum(theta)) for j in theta]
        encodings.append(code)
    return np.array(encodings, dtype=float), header

# =============================================================================
# 2. ПОДГОТОВКА И ЗАПУСК
# =============================================================================

def generate_features_pipeline(fasta_file):
    print(f"Generating features for {fasta_file}...")
    fasta = read_protein_sequences(fasta_file)
    
    # 1. AAC
    feat0, _ = AAC(fasta)
    allfeat = feat0
    
    # 2. DPC
    feat1, _ = DPC(fasta, 0)
    allfeat = np.concatenate((allfeat, feat1), axis=1)
    
    # 3. APAAC
    feat2, _ = APAAC(fasta)
    allfeat = np.concatenate((allfeat, feat2), axis=1)
    
    # 4. CTDC
    feat3, _ = CTDC(fasta)
    allfeat = np.concatenate((allfeat, feat3), axis=1)
    
    # 5. CTDD
    feat4, _ = CTDD(fasta)
    allfeat = np.concatenate((allfeat, feat4), axis=1)
    
    # 6. CTDT
    feat5, _ = CTDT(fasta)
    allfeat = np.concatenate((allfeat, feat5), axis=1)
    
    # 7. GAAC
    feat6, _ = GAAC(fasta)
    # GAAC теперь не возвращает имя, поэтому берем целиком, без среза [:, 1:]
    allfeat = np.concatenate((allfeat, feat6), axis=1)
    
    # 8. KSCTriad
    feat7, _ = KSCTriad(fasta)
    # KSCTriad теперь не возвращает имя, берем целиком
    allfeat = np.concatenate((allfeat, feat7), axis=1)
    
    # 9. CTriad
    feat8, _ = CTriad(fasta)
    allfeat = np.concatenate((allfeat, feat8), axis=1)
    
    # 10. DDE
    feat9, _ = DDE(fasta)
    allfeat = np.concatenate((allfeat, feat9), axis=1)
    
    # 11. PAAC
    feat10, _ = PAAC(fasta)
    allfeat = np.concatenate((allfeat, feat10), axis=1)
    
    # Вычисляем индексы столбцов для Stacking
    header_lens = [
        len(feat0[0]), len(feat1[0]), len(feat2[0]), len(feat3[0]),
        len(feat4[0]), len(feat5[0]), len(feat6[0]), len(feat7[0]),
        len(feat8[0]), len(feat9[0]), len(feat10[0])
    ]
    
    f = []
    before = 0
    # В оригинале цикл идет до 10 (последняя группа PAAC игнорируется в ансамбле)
    for length in header_lens[:-1]: 
        after = before + length
        f.append(list(range(before, after)))
        before = after
        
    return allfeat, f, fasta

if __name__ == "__main__":
    # --- 1. АРГУМЕНТЫ КОМАНДНОЙ СТРОКИ ---
    if len(sys.argv) > 1:
        TARGET_FILE = sys.argv[1]
    else:
        TARGET_FILE = 'blind_test_dataset.fasta'
    
    # ВТОРОЙ АРГУМЕНТ: Путь для сохранения (опционально)
    if len(sys.argv) > 2:
        OUTPUT_FILE = sys.argv[2]
    else:
        # Автоматическое имя, если второй аргумент не дан
        input_filename = os.path.basename(TARGET_FILE)
        OUTPUT_FILE = f"Results_{input_filename}.csv"

    print(f"Target file: {TARGET_FILE}")
    print(f"Output will be saved to: {OUTPUT_FILE}")
    
    if not os.path.exists(TARGET_FILE):
        print(f"Error: Target file {TARGET_FILE} not found!")
        sys.exit(1)

    # --- 2. ЗАГРУЗКА ТРЕНИРОВОЧНЫХ ДАННЫХ ---
    print("Loading Training Data (Required for Stacking)...")
    try:
        allfeat_pos, f, _ = generate_features_pipeline('data/TR_P_132.fasta')
        allfeat_neg, _, _ = generate_features_pipeline('data/TR_N_305.fasta')
    except Exception as e:
        print(f"Error: {e}")
        print("Make sure you are in the directory containing 'data/' folder.")
        sys.exit(1)

    X_train = np.concatenate((allfeat_pos, allfeat_neg), axis=0)
    y_train = np.concatenate((np.ones(len(allfeat_pos)), np.zeros(len(allfeat_neg))))

    # --- 3. ПОДГОТОВКА ЦЕЛЕВОГО ФАЙЛА ---
    try:
        X_target, _, target_sequences = generate_features_pipeline(TARGET_FILE)
    except Exception as e:
        print(f"Error processing target file: {e}")
        sys.exit(1)

    # --- 4. STACKING (ГЕНЕРАЦИЯ МЕТА-ПРИЗНАКОВ) ---
    print("Running Stacking Ensemble...")
    featx = []
    for i in range(10): 
        Xs = X_train[:, f[i]]
        Xts = X_target[:, f[i]]
        clfs = [
            RandomForestClassifier(n_estimators=500, random_state=0),
            ExtraTreesClassifier(n_estimators=500, random_state=0),
            SVC(probability=True, random_state=0),
            LogisticRegression(random_state=0, max_iter=5000),
            XGBClassifier(use_label_encoder=False, eval_metric='logloss'),
            KNeighborsClassifier(weights="distance", algorithm="auto")
        ]
        for clf in clfs:
            clf.fit(Xs, y_train)
            pr = clf.predict_proba(Xts)[:, 0]
            feat = np.reshape(pr, (len(Xts), 1))
            if len(featx) == 0: featx = feat
            else: featx = np.concatenate((featx, feat), axis=1)

    # --- 5. ФИНАЛЬНОЕ ПРЕДСКАЗАНИЕ ---
    print("Final Prediction...")
    mask = [2,3,4,5,9,10,12,14,15,18,24,25,28,34,46,47,48,53,56,57]
    try:
        Selected_feat = featx[:, mask]
    except IndexError:
        print("Error: Generated feature matrix has wrong dimensions.")
        sys.exit(1)

    model_path = "model/pima.pickle_model_svm_PF.dat"
    if not os.path.exists(model_path):
        print(f"Error: Model file {model_path} not found.")
        sys.exit(1)
    ldmodel = pickle.load(open(model_path, "rb"))
    final_probs = ldmodel.predict_proba(Selected_feat)
    threshold = 0.5

    results = []
    for i in range(len(target_sequences)):
        seq_id = target_sequences[i][0]
        amy_prob = final_probs[i][1] 
        prediction = "AMYLOID" if amy_prob > threshold else "Non-Amyloid"
        results.append([seq_id, amy_prob, prediction])

    # --- 6. ВЫВОД И СОХРАНЕНИЕ ---
    headers = ["Sequence ID", "Amyloid Prob", "Prediction"]
    print(tabulate(results, headers=headers, tablefmt="pretty"))

    # Создаем папку для вывода, если её нет
    output_dir = os.path.dirname(OUTPUT_FILE)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    df = pd.DataFrame(results, columns=headers)
    df.to_csv(OUTPUT_FILE, index=False)
    print(f"\nResults saved to {OUTPUT_FILE}")
