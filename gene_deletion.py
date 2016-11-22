# Notes for me

# Removing the second glucose exp for anaerobe. It yields a biomass of 0 and I am
# not sure why.

from psamm.lpsolver import generic

solver = generic.Solver()

import sys
sys.path.insert(0, '.')
import sys

import moma_test

import csv

# Run the FBA tests with a certain cutoff for viable/not viable
def run_tests(cutoff):
    gene_dict = {}
    # The medium that we are limiting in each of the experiments
    medium_ids = ['R_EX_glc_e_', 'R_EX_gal_e_', 'R_EX_glyc_e_', 'R_EX_etoh_e_']

    aerobe = True # Our first glucose exp will be aerobe

    for exp in medium_ids:
        # We are working on a glucose exp
        if exp == 'R_EX_glc_e_':
            if aerobe == True:
                mm.limits['R_EX_o2_e_'].lower = -2.0
                mm.limits['R_EX_co2_e_'].lower = 0.0
                aerobe = False # set up for the next experiment
            # Exp for anaerobe
            else:
                mm.limits['R_EX_o2_e_'].lower = 0.0
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

            # If the flux is greater then the cutoff, then we have a viable
            if flux_percent >= cutoff:
                flux = "1"
            # The knockout was not viable
            else:
                flux = "0"
            # We have already seen this gene in an exp
            if gene in gene_dict:
                gene_dict[gene] += [flux] # Append to the end of the lise. Use round() to round up and down
            # This is only for the first experiment
            else:
                gene_dict[gene] = [flux]
        mm.limits[exp].lower =  0.0 # Reset the exp medium to 0
    return gene_dict

def output_results(gene_dict):
    writer = csv.writer(open('raw_results.csv', 'w+'))
    # Header for the table
    writer.writerow([" ", "No.", "Deletion", "Glucose (aerobe)", "Glucose (anaerobe)",  "Galactose", "Glycerol", "Ethanol"])
    i = 1 # Gives each gene deletion a number
    for gene, flux_percent in gene_dict.iteritems():
        writer.writerow([" ", i, gene, flux_percent[0], " ", flux_percent[1], flux_percent[2], flux_percent[3]])
        i += 1 # iterate reaction counter
import csv

# Currently need to fix the float and round issue. Going to round with the exp python script
def read_full_gene(filename):
    gene_dictionary = {}
    f = open(filename)
    csv_f = csv.reader(f)
    next(csv_f, None) # skip the header
    for row in csv_f:
        # key: gene
        # value: [glucose aerobe, glucose anaerobe, galactose, glycerol, ethanol]
        # row[3] and row[9] are for anaerobe
        #print(row)
        gene_dictionary[row[1]] = [row[2], row[4], row[5], row[6], row[8], row[10],row[11],row[12]]
    return gene_dictionary

def clean_dictionary(gene_dictionary):
    clean_dictionary = {}
    for gene, flux in gene_dictionary.iteritems():
        clean_dictionary[gene] = []
        for exp_value in flux:
            if exp_value == '0' or exp_value == '0.0':
                clean_dictionary[gene] += [0]
            elif exp_value == "1" or exp_value == "1.0":
                clean_dictionary[gene] += [1]
            elif exp_value == "m" or exp_value == "-":
                clean_dictionary[gene] += ["-"]
            else:
                clean_dictionary[gene] += [1]
    #print(clean_dictionary)
    return clean_dictionary

def calc_stat(gene_dictionary):
    tests = [0, 1, 2, 3]
    medium_dictionary = {}
    for medium in tests:
        true_live = 0
        false_live = 0
        true_death = 0
        false_death = 0
        for gene, flux in gene_dictionary.iteritems():
            if flux[medium] == flux[medium + 4]:
                if flux[medium] == 0:
                    true_death += 1
                else:
                    true_live += 1
            elif flux[medium + 4] == "-":
                i = 0 # do nothing
            else:

                if flux[medium] == 0:
                    false_death += 1
                else:
                    false_live += 1
        medium_dictionary[medium] = [true_live, false_live, true_death, false_death]
    return medium_dictionary

def read_gene_file(filename):
    gene_dictionary = {}
    f = open(filename)
    csv_f = csv.reader(f)
    next(csv_f, None) # skip the header
    for row in csv_f:
        # key: gene
        # value: [glucose aerobe, glucose anaerobe, galactose, glycerol, ethanol]
        gene_dictionary[row[2]] = [row[3], row[5], row[6], row[7]]
    return gene_dictionary

def combine_gene(my_exp, compare_exp):
    for gene, exp in my_exp.iteritems():
        if gene in compare_exp:
            my_exp[gene] += compare_exp[gene]
        else:
            my_exp[gene] += ["-", "-", "-", "-", "-"] # Put in a dash to show this was never tested
    return my_exp

#gene_dictionary = read_gene_file("compare_exp.csv")
compare_exp = read_gene_file("paper_exp.csv")

for i in range(100):
    # Test every possible cutoff between 0 and 1
    cutoff = float(i / 100)
    gene_dict = run_tests(cutoff)
    gene_dictionary = combine_gene(gene_dict, compare_exp)
    gene_dictionary = clean_dictionary(gene_dictionary)
    medium_dictionary = calc_stat(gene_dictionary)
    print("Cutoff " + str(cutoff) + ": ")
    print(medium_dictionary)


#output_results(gene_dict)
