import pandas as pd

# convert excel file to csv file
# transpose data before converting!
read_file = pd.read_excel("QC_AATT_127g_6671m.xlsx")
# Write the dataframe object into csv file
read_file.to_csv("QC_AATT_127g_6671m.csv",
                 index=False)

# transform objects to corresponding string values
df = pd.read_csv("QC_AATT_127g_6671m.csv")
df = df.astype("|S")
print(df.head())
# save
df.to_csv("QC_AATT_127g_6671m-1.csv", index=False)
# delete "b" in the saved csv file

# split words
lines = []
with open('QC_AATT_127g_6671m-1.csv', 'r') as f:
    for line in f.readlines():
        l = line.strip().split(',')
        lines.append((l))

# remove-square-brackets-from-list
for x in lines:
    # remove square brackets from list
    s = (", ".join(str(e) for e in x))
    print(list(s))
# python s.py > log55.csv