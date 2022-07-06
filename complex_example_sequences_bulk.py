import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#using the dataset compiled by Si and Yan 2022 (https://doi.org/10.1101/2021.12.21.473437 )

#TODO: get those PDB chain / Uniprot query scripts I wrote from AWS & Use those to scrape RCSB PDB for the relevant complex info + uniprot IDs to see what ones can be queried from OMA DB

test_set_complexes = pd.read_excel('test_complexes.xlsx', 'PPIs_after_AF2')
print (test_set_complexes)