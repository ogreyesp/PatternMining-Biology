import pandas as pd

def convert_to_sets(dataframe):
    return dataframe.apply(
        lambda row: frozenset(row.index.values + "_" + row.values), axis=1)


def compute_supports(pattern, serie):
    """This funtion computes the support of a pattern over a Series pandas which has a multiple frozensets"""
    return serie.apply(lambda set: set.issuperset(pattern)).sum() / serie.size
