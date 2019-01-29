import pandas as pd

def convert_to_sets(dataframe):

    return dataframe.apply(
        lambda row: list(row.index.values[row.values.astype(bool)]), axis=1)

dfT = pd.read_csv("datasets/fileA.csv")

dfT= convert_to_sets(dfT)

dfT.to_csv("blabla.csv", index=False)

