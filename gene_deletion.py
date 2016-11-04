from psamm.lpsolver import generic

solver = generic.Solver()

import sys
sys.path.insert(0, '.')
import sys

import moma_test

gene_dict = {}


# The medium that we are limiting in each of the experiments
medium_ids = ['R_EX_glc_e_', 'R_EX_glc_e_', 'R_EX_gal_e_', 'R_EX_glyc_e_', 'R_EX_etoh_e_']

aerobe = True # Our first glucose exp will be aerobe

for exp in medium_ids:
    # We are working on a glucose exp
    if exp == 'R_EX_glc_e_':
        if aerobe == True:
            mm.limits['R_EX_o2_e_'].lower = -2.0
            mm.limits['R_EX_co2_e_'].lower = 0.0
            aerobe = False # set up for the next experiment
        else:
            mm.limits['R_EX_o2_e_'].lower = -2.0
            mm.limits['R_EX_co2_e_'].lower = -10.0
    # Every other exp we give both
    else:
        mm.limits['R_EX_o2_e_'].lower = -2.0
        mm.limits['R_EX_co2_e_'].lower = -10.0
    # We are giving some the exp nutrient
    mm.limits[exp].lower = -10.0
    print("Working on exp " + exp)
    # Run the test
    simulation = dict(moma_test.gene_deletion(model, mm, solver))
    # Put the experiment in a dictionary for output later
    for gene, flux_percent in simulation.iteritems():
        # We have already seen this gene in an exp
        if gene in gene_dict:
            gene_dict[gene] += [round(flux_percent)] # Append to the end of the list
        # This is only for the first experiment
        else:
            gene_dict[gene] = [round(flux_percent)]
    mm.limits[exp].lower =  0.0 # Reset the exp medium to 0

import csv

def output_results(gene_dict):
    writer = csv.writer(open('results.csv', 'w+'))
    # Header for the table
    writer.writerow([" ", "No.", "Deletion", "Glucose (aerobe)", "Glucose (anaerobe)", "Galactose", "Glycerol", "Ethanol"])
    i = 1 # Gives each gene deletion a number
    for gene, flux_percent in gene_dict.iteritems():
        writer.writerow([" ", i, gene, flux_percent[0], flux_percent[1], flux_percent[2], flux_percent[3], flux_percent[4]])
        i += 1 # iterate reaction counter
