import csv


def read_gene_file(filename):
    gene_dictionary = {}
    f = open(filename)
    csv_f = csv.reader(f)
    next(csv_f, None) # skip the header
    for row in csv_f:
        # key: gene
        # value: [glucose aerobe, glucose anaerobe, galactose, glycerol, ethanol]
        gene_dictionary[row[2]] = [row[3], row[4], row[5], row[6], row[7]]
    return gene_dictionary

def compare(compare_exp, my_exp):
    for gene, exp in my_exp.iteritems():
        my_exp[gene] += [" "] # This is just to create a divide between our exp and theirs
        if gene in compare_exp:
            my_exp[gene] += compare_exp[gene]
        else:
            my_exp[gene] += ["-", "-", "-", "-", "-"] # Put in a dash to show this was never tested

    return my_exp

def output_results(gene_dict):
    writer = csv.writer(open('compare_exp.csv', 'w+'))
    # Header for the table
    writer.writerow(["No.", "Deletion", "Glucose (aerobe)", "Glucose (anaerobe)", "Galactose", "Glycerol", "Ethanol", " ", "Glucose (aerobe)", "Glucose (anaerobe)", "Galactose", "Glycerol", "Ethanol"])
    i = 1 # Gives each gene deletion a number
    for gene, flux_percent in gene_dict.iteritems():
        writer.writerow([i, gene, flux_percent[0], flux_percent[1], flux_percent[2], flux_percent[3], flux_percent[4],flux_percent[5], flux_percent[6], flux_percent[7], flux_percent[8], flux_percent[9], flux_percent[10]])
        i += 1 # iterate reaction counter


my_exp = read_gene_file("results.csv")
compare_exp = read_gene_file("compare.csv")

compared_exp = compare(compare_exp, my_exp)

output_results(compared_exp)
