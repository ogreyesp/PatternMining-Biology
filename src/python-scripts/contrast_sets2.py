import sys
import pandas as pd
import numpy as np
from itertools import combinations

def convert_to_sets(dataframe):
    return dataframe.apply(
        lambda row: frozenset(row.index.values + "=" + row.values), axis=1)

def compute_supports(pattern, serie):
    """This funtion computes the support of a pattern over a Series pandas which has a multiple frozensets"""
    
    return serie.apply(lambda set: set.issuperset(pattern)).sum() / serie.size

# Dict conatining name and path of the datasets
names_files = {
    "KICH": ("datasets/discretizaen5kich.csv", "datasets/rinonprueba.csv"),
    "KIRC": ("datasets/discretizaen5kirc.csv", "datasets/prueba0.csv"),
    "KIRP": ("datasets/discretizaen5kirp.csv", "datasets/prueba1.csv"),
}

#List of dataframes with frequent item sets
list_frequentsets = {}

#List of dataframes converted to a series of frozensets
list_all_itemsets = {}

all_patterns = set()

number_patterns = 0

# For each dataset the apriori algorithm is executed and the datasets are also converted to a more convenient way for searching
for name, fileN in names_files.items():

    dataframe = pd.read_csv(fileN[0])

    #Drop the class column if it exists
    if "Class" in dataframe.columns:
        dataframe.drop(['Class'], axis=1, inplace=True)

    list_all_itemsets[name] = convert_to_sets(dataframe)
    
    frequent_itemsets = pd.read_csv(fileN[1])

    if name == 'KICH':
        frequent_itemsets.drop("Unnamed: 0", axis=1, inplace=True)
    
    number_patterns += frequent_itemsets.shape[0]

    #Create a new column with the length of the patterns
    frequent_itemsets['items'] = frequent_itemsets['items'].apply(
        lambda x: frozenset(x[1:-1].split(',')))

    #Create a new column with the length of the patterns
    frequent_itemsets['size'] = frequent_itemsets['items'].apply(
        lambda x: len(x))

    list_frequentsets[name] = frequent_itemsets

    all_patterns.update(frequent_itemsets['items'].values)

print("The number of unique global patterns is", len(all_patterns))
print("The number of total patterns is:", number_patterns, " and there are",
      (number_patterns - len(all_patterns)), "repeated patterns")

# Computing the constrast sets

# For each global pattern

dict_support = {}
dict_contrast_set = {}
for pattern in all_patterns:
    # Compute the support of the pattern in each dataset
    for name in names_files.keys():
        support = 0
        # To speed-up the computation, first we check if the support of the pattern already exists and, therefore, the computation is avoided
        frequent_itemsets = list_frequentsets[name]
        frequent_itemsets = frequent_itemsets[frequent_itemsets['items'] ==
                                              pattern]
        # If the number of selected rows is greater than cero then the pattern exists for this dataset
        if not frequent_itemsets.empty:
            support = frequent_itemsets.iat[0, 1]  # getting the support
        else:  # compute the support on the dataset
            support = compute_supports(pattern, list_all_itemsets[name])
        dict_support[(pattern, name)] = support
    # Computing the value of the contrast set
    dict_contrast_set[pattern] = max([
        abs(dict_support[(pattern, comb[0])] - dict_support[(pattern, comb[1])])
        for comb in combinations(names_files.keys(), 2)
    ])

contrast_sets = pd.DataFrame(
    list(dict_contrast_set.items()), columns=['Pattern', 'Value'])

contrast_sets.set_index("Pattern", inplace=True)

# Insert one column for each dataset to represent the supports on these datasets
serie = pd.Series(dict_support).reset_index()
serie.columns = ["Pattern", "Dataset", "Support"]

for name_dataset in names_files.keys():
    filter_dataset = serie[serie["Dataset"] == name_dataset][[
        "Pattern", "Support"
    ]]
    filter_dataset.columns = ["Pattern", name_dataset + "-Support"]
    filter_dataset.set_index("Pattern", inplace=True)
    contrast_sets = contrast_sets.join(filter_dataset)

#Sort the patterns by their values
contrast_sets.sort_values(
    ['Value'], ascending=[False], inplace=True)

contrast_sets.reset_index(inplace=True)

contrast_sets.to_csv("contrast_setskidney.csv", index=False)

print("Done!!")
