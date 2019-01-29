import pandas as pd
import numpy as np
from mlxtend.frequent_patterns import apriori
from utils import convert_to_sets
from utils import compute_supports

dfT = pd.read_csv("../datasets/discretizaen5T.csv")
dfN = pd.read_csv("../datasets/discretizaen5N.csv")

min_support = 0.32

#Drop the class column
dfT.drop(['Class'], axis=1, inplace=True)
dfN.drop(['Class'], axis=1, inplace=True)

#Running a priori algorithm in class T
# Before running apriori the dataset must be transformed in a onehot representation with columns containing boolean values
dfT_OneHot = pd.get_dummies(dfT)

dfT_OneHot = dfT_OneHot[dfT_OneHot.columns].astype(bool)

frequent_itemsets_T = apriori(
    dfT_OneHot, min_support=min_support, use_colnames=True)

#Create a new column with the length of pattern
frequent_itemsets_T['length'] = frequent_itemsets_T['itemsets'].apply(
    lambda x: len(x))

# Transform the other dataset
dfN = convert_to_sets(dfN)

#Computing the support of each pattern encountered in the first class in the second class
nrowN = dfN.size

# Create the dataframe of emerging patterns
emerging_patterns = frequent_itemsets_T[['itemsets', 'length', 'support']]
#Changing the name of the columns
emerging_patterns.columns = ['Pattern', 'Size', 'supportT']

#One way to insert new columns
#emerging_patterns = emerging_patterns.reindex(columns=['Pattern', 'Size', 'supportT', 'supportN', 'GrowthRatio'])
emerging_patterns = emerging_patterns.assign(supportN=0.0, GrowthRatio=np.inf)

for pattern in emerging_patterns.itertuples():

    supportN = compute_supports(pattern.Pattern, dfN)

    if supportN > 0:
        emerging_patterns.at[pattern.Index, 'supportN'] = supportN
        emerging_patterns.at[
            pattern.Index,
            'GrowthRatio'] = emerging_patterns.at[pattern.Index,
                                                  'supportT'] / supportN

#Sort the patterns by their support values and then by their sizes
emerging_patterns.sort_values(
    ['GrowthRatio', 'supportT', 'Size'],
    ascending=[False, False, False],
    inplace=True)

emerging_patterns.to_csv("../results/emerging_patterns.csv", index=False)

print("Done!")