import pandas as pd
from pandas import read_csv
import os

os.chdir("D:/project1/Informativeness")
# Convert Excel to the CSV
read_file = pd.read_excel(
    "J_K_US_AABB_CC_imputed3_chr_inputfile.xlsx",
    "Blad3")
read_file.to_csv(
    "J_K_US_AABB_CC_imputed3_chr_inputfile.csv",
    index=None,
    header=True)



import itertools
import pandas as pd

df = pd.read_csv(
    "J_K_US_AABB_CC_imputed3_chr_inputfile.csv"
)
# Creating an empty list
res = []

# Iterating through the columns of dataframe
for column in df.columns:
    # Storing the rows of a column
    # into a temporary list
    li = df[column].tolist()

    # appending the temporary list
    res.append(li)

# Printing the final list
# print(res)

for each_list in res:
    letters = each_list
    prev = letters[0]
    result = []
    grouped = [(len(list(g)) - 1, k) for k, g in (itertools.groupby(letters))]
    weird_transitions = [grouped[0]] + [(n + 1, b) for (n, (
        a, b)) in enumerate(zip(letters, letters[1:])) if a != b and b != 'X']
    # print(weird_transitions)
    print(len(weird_transitions))
# save it as csv file
# python x.py > log5.csv

