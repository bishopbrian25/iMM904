import csv

# Currently need to fix the float and round issue. Going to round with the exp python script
def read_gene_file(filename):
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
    print(clean_dictionary)
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
                    print("HI")
                    false_death += 1
                else:
                    false_live += 1
        medium_dictionary[medium] = [true_live, false_live, true_death, false_death]
    return medium_dictionary


def export(medium_dictionary):
    writer = csv.writer(open('stats.csv', 'w+'))
    # Header for the table
    writer.writerow(["Medium", "True Live", "False Live", "True Death", "False Death"])
    i = 1 # Gives each gene deletion a number
    for medium, stats in medium_dictionary.iteritems():
        writer.writerow([medium, stats[0], stats[1], stats[2], stats[3]])
        i += 1 # iterate reaction counter

gene_dictionary = read_gene_file("compare_exp.csv")
gene_dictionary = clean_dictionary(gene_dictionary)
medium_dictionary = calc_stat(gene_dictionary)



print(medium_dictionary)
export(medium_dictionary)
