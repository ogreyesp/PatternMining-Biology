import pandas as pd
import numpy as np
from mlxtend.frequent_patterns import apriori
from src.utils import convert_to_sets
from src.utils import compute_supports
from itertools import combinations

# Dict conatining name and path of the datasets
names_files = {
    "Basal": "datasets/discretizaen5basal.csv",
    "Her2": "datasets/discretizaen5her2.csv",
    "Luma": "datasets/discretizaen5luma.csv",
    "Lumb": "datasets/discretizaen5lumb.csv"
}

#List of dataframes with frequent item sets
list_frequentsets = {}

#List of dataframes converted to a series of frozensets
list_all_itemsets = {}

all_patterns = set()

min_support = 0.9

number_patterns = 0

# For each dataset the apriori algorithm is executed and the datasets are also converted to a more convenient way for searching
for name, file in names_files.items():

    dataframe = pd.read_csv(file)

    #Drop the class column if it exists
    if "Class" in dataframe.columns:
        dataframe.drop(['Class'], axis=1, inplace=True)

    list_all_itemsets[name] = convert_to_sets(dataframe)

    # Before running apriori the dataset must be transformed in a onehot representation with columns containing boolean values
    dataframe = pd.get_dummies(dataframe)
    dataframe = dataframe[dataframe.columns].astype(bool)

    #Running a priori algorithm in dataset i-th
    print("Running apriori algorithm on dataset ", name)
    frequent_itemsets = apriori(
        dataframe, min_support=min_support, use_colnames=True)

    print("Patterns in dataset", name, ": ", frequent_itemsets.shape[0])

    number_patterns += frequent_itemsets.shape[0]

    #Create a new column with the length of the patterns
    frequent_itemsets['size'] = frequent_itemsets['itemsets'].apply(
        lambda x: len(x))

    list_frequentsets[name] = frequent_itemsets

    all_patterns.update(frequent_itemsets['itemsets'].values)

print("The number of unique global patterns is", len(all_patterns))
print("The number of total patterns is:", number_patterns, " and there is ",
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

        frequent_itemsets = frequent_itemsets[frequent_itemsets['itemsets'] ==
                                              pattern]

        # If the number of selected rows is greater than cero then the pattern exists for this dataset
        if not frequent_itemsets.empty:

            support = frequent_itemsets.iat[0, 0]  # getting the support

        else:  # compute the support on the dataset
            support = compute_supports(pattern, list_all_itemsets[name])

        dict_support[(pattern, name)] = support

    # Computing the value of the contrast set
    dict_contrast_set[pattern] = max([
        abs(dict_support[(pattern, comb[0])] - dict_support[(pattern,
                                                             comb[1])])
        for comb in combinations(names_files.keys(), 2)
    ])

contrast_sets = pd.DataFrame(
    list(dict_contrast_set.items()), columns=['Pattern', 'Value'])

#Create a new column with the length of pattern
contrast_sets['Size'] = contrast_sets['Pattern'].apply(lambda x: len(x))

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
    ['Value', 'Size'], ascending=[False, False], inplace=True)

contrast_sets.reset_index(inplace=True)

contrast_sets.to_csv("results/contrast_sets.csv", index=False)

print("Done!!")