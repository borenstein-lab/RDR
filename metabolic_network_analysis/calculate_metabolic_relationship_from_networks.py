# !/usr/bin/env python

# from __future__ import print_function
import sys
import SeedSet
import NetCooperate
import csv
import pandas as pd
import numpy as np
import argparse
import os
import csv

parser = argparse.ArgumentParser(description="Calculate pairwise metabolic relationships")

parser.add_argument("networks_folder",
                    help="Path to the folder of metabolic network files, should end with network.tab")
args = parser.parse_args()

# Our graph reading script
#  Reads a tab-separated file of edges

def readGraph(graphFile):
    graph = {}
    for line in open(graphFile):

        vals = line.strip().split("\t");

        if len(vals[0]) == 0 or vals[0][0] == '#':
            continue
        elif len(vals) != 2:
            print("Bad line: " + line)
            continue

        graph.setdefault(vals[0], []).append(vals[1])

    return graph




# Import all files in the folder
directory = os.fsencode(args.networks_folder)

seed_list_all_taxa = []

# Iterate over network files in specified folder (does not look recursively into subdirectories)
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    full_path = os.path.join(directory, file)
    if filename.endswith(".tab"):
        graph = readGraph(full_path)
        seed, seed_groups, non_seeds, pruned, nodes = SeedSet.calculate_seeds(graph)
        seed_taxa = {"taxa_name": str.replace(filename, ".tab", ""),
                     "seed": seed,
                     "seed_groups": seed_groups,
                     "non_seeds": non_seeds,
                     "pruned": pruned,
                     "nodes": nodes}
        seed_list_all_taxa.append(seed_taxa)
else:
    (print("Seed set calculated for " + str(len(seed_list_all_taxa)) + " taxa\n"))

with open('results.csv', mode='w', newline='') as f:
    keys = ['taxa_one', 'taxa_two', 'complementarity', 'competition', "complementarity_support", "competition_support"]
    writer = csv.DictWriter(f, fieldnames=keys)
    writer.writeheader()    # add column names in the CSV file
    for i in range(len(seed_list_all_taxa)):
        #print(seed_list_all_taxa[i].get('taxa_name'))
        for j in range(len(seed_list_all_taxa)):
            complementarity,  complementarity_support = NetCooperate.compute_single_interaction_score(seed_list_all_taxa[i].get('seed_groups'), seed_list_all_taxa[j].get('non_seeds'))
            competition, competition_support = NetCooperate.compute_single_interaction_score(seed_list_all_taxa[i].get('seed_groups'), seed_list_all_taxa[j].get('seed'))
            taxa_one = seed_list_all_taxa[i].get('taxa_name')
            taxa_two = seed_list_all_taxa[j].get('taxa_name')
            writer.writerow({'taxa_one': taxa_one,
                             'taxa_two': taxa_two,
                             'complementarity' : complementarity,
                             'competition' : competition,
                             'complementarity_support' : complementarity_support,
                             'competition_support': competition_support})
f.close()
