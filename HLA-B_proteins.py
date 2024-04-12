########################################################################
### Import protein sequences
B_prot = open("./B_prot.txt")
proteins = B_prot.read()

### Subset portion of sequences containing position 116
target = proteins.split("Prot")
target = target[2]

### Separate each line of sequences
target1 = target.split("B")

### Create list of amino acids at position 116
AAs = []
for i in range(len(target1)):
    if len(target1[i]) < 88:
        continue
    AAs.append(target1[i][88])

### Clean up list of amino acids
AAs = [x.replace("X", " ") for x in AAs]
AAs = [y.replace("-", "Y") for y in AAs]
AAs.remove("/n")
AAs = ' '.join(AAs).split()

### Identify unique amino acids in list
AAs_set = set(AAs)

### Count amino acid frequency at postion 116
counts = []
for z in list(AAs_set):
    counts.append(AAs.count(z))

### Create a df with amino acid frequency and physical property
import pandas as pd
data = pd.DataFrame({
    'Amino Acid': list(AAs_set), 'Counts': counts,
                     })
data = data.sort_values('Counts', ascending=False)
data['Charge'] = ['neutral', 'neutral', 'negative', 'neutral', 'neutral', 'neutral', 'neutral', 'positive', 'neutral',
           'positive', 'neutral', 'neutral']

### Plot amino acid frequency
import matplotlib.pyplot as plt
plt.hist(x=data["Amino Acid"], weights=data["Counts"])
plt.ylim(0, 3500)
plt.title("Amino Acid Frequency at Position 116", loc='left')
plt.xlabel("Amino Acid")
plt.ylabel("Count")
plt.show()

### Plot amino acid frequency with physical property
plt.hist(data[data['Charge'] == 'neutral']['Amino Acid'], weights=data[data['Charge'] == 'neutral']["Counts"],
         color = 'green', label='Neutral Charge')
plt.hist(data[data['Charge'] == 'negative']['Amino Acid'], weights=data[data['Charge'] == 'negative']["Counts"],
         color = 'red', label='Negative Charge')
plt.hist(data[data['Charge'] == 'positive']['Amino Acid'], weights=data[data['Charge'] == 'positive']["Counts"],
         color = 'blue', label='Positive Charge')
plt.ylim(0, 3500)
plt.title("Amino Acid Frequency at Position 116 w/ Charge", loc='left')
plt.xlabel("Amino Acid")
plt.ylabel("Count")
plt.legend()
plt.show()

### Create list of proteins for each amino acid at position 116
target_Y = []
for a in range(len(target1)):
    if len(target1[a]) >= 88 and target1[a][88] == "Y":
        target_Y.append(target1[a][0:3])
    if len(target1[a]) >= 88 and target1[a][88] == "-":
        target_Y.append(target1[a][0:3])
target_Y = sorted(set(target_Y))
target_Y = [str("HLA-B" + target_Y[a]) for a in range(len(target_Y))]
target_Y = [b.replace("*", "") for b in target_Y]

target_S = []
for a in range(len(target1)):
    if len(target1[a]) >= 88 and target1[a][88] == "S":
        target_S.append(target1[a][0:3])
target_S = sorted(set(target_S))
target_S = [str("HLA-B" + target_S[a]) for a in range(len(target_S))]
target_S = [b.replace("*", "") for b in target_S]

target_D = []
for a in range(len(target1)):
    if len(target1[a]) >= 88 and target1[a][88] == "D":
        target_D.append(target1[a][0:3])
target_D = sorted(set(target_D))
target_D = [str("HLA-B" + target_D[a]) for a in range(len(target_D))]
target_D = [b.replace("*", "") for b in target_D]

target_F = []
for a in range(len(target1)):
    if len(target1[a]) >= 88 and target1[a][88] == "F":
        target_F.append(target1[a][0:3])
target_F = sorted(set(target_F))
target_F = [str("HLA-B" + target_F[a]) for a in range(len(target_F))]
target_F = [b.replace("*", "") for b in target_F]

target_L = []
for a in range(len(target1)):
    if len(target1[a]) >= 88 and target1[a][88] == "L":
        target_L.append(target1[a][0:3])
target_L = sorted(set(target_L))
target_L = [str("HLA-B" + target_L[a]) for a in range(len(target_L))]
target_L = [b.replace("*", "") for b in target_L]

target_T = []
for a in range(len(target1)):
    if len(target1[a]) >= 88 and target1[a][88] == "T":
        target_T.append(target1[a][0:3])
target_T = sorted(set(target_T))
target_T = [str("HLA-B" + target_T[a]) for a in range(len(target_T))]
target_T = [b.replace("*", "") for b in target_T]

target_P = []
for a in range(len(target1)):
    if len(target1[a]) >= 88 and target1[a][88] == "P":
        target_P.append(target1[a][0:3])
target_P = sorted(set(target_P))
target_P = [str("HLA-B" + target_P[a]) for a in range(len(target_P))]
target_P = [b.replace("*", "") for b in target_P]

target_H = []
for a in range(len(target1)):
    if len(target1[a]) >= 88 and target1[a][88] == "H":
        target_H.append(target1[a][0:3])
target_H = sorted(set(target_H))
target_H = [str("HLA-B" + target_H[a]) for a in range(len(target_H))]
target_H = [b.replace("*", "") for b in target_H]

target_N = []
for a in range(len(target1)):
    if len(target1[a]) >= 88 and target1[a][88] == "N":
        target_N.append(target1[a][0:3])
target_N = sorted(set(target_N))
target_N = [str("HLA-B" + target_N[a]) for a in range(len(target_N))]
target_N = [b.replace("*", "") for b in target_N]

target_R = []
for a in range(len(target1)):
    if len(target1[a]) >= 88 and target1[a][88] == "R":
        target_R.append(target1[a][0:3])
target_R = sorted(set(target_R))
target_R = [str("HLA-B" + target_R[a]) for a in range(len(target_R))]
target_R = [b.replace("*", "") for b in target_R]

target_A = []
for a in range(len(target1)):
    if len(target1[a]) >= 88 and target1[a][88] == "A":
        target_A.append(target1[a][0:3])
target_A = sorted(set(target_A))
target_A = [str("HLA-B" + target_A[a]) for a in range(len(target_A))]
target_A = [b.replace("*", "") for b in target_A]

target_I = []
for a in range(len(target1)):
    if len(target1[a]) >= 88 and target1[a][88] == "I":
        target_I.append(target1[a][0:3])
target_I = sorted(set(target_I))
target_I = [str("HLA-B" + target_I[a]) for a in range(len(target_I))]
target_I = [b.replace("*", "") for b in target_I]

### Combine protein lists into df and export as csv
export = pd.DataFrame({'Y': pd.Series(target_Y), 'S': pd.Series(target_S), 'D': pd.Series(target_D), 'F': pd.Series(target_F), 'L': pd.Series(target_L), 'T': pd.Series(target_T),
                       'P': pd.Series(target_P), 'H': pd.Series(target_H), 'N': pd.Series(target_N), 'R': pd.Series(target_R), 'A': pd.Series(target_A), 'I': pd.Series(target_I)})
export.to_csv("./proteins_by_AA.csv", index=False)


