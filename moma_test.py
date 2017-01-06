import csv

from six import iteritems

from psamm.expression import boolean
from psamm import fluxanalysis
import moma

# Read in an expirement
# I have never used this before
def read_experiments(f):
    header = f.readline().strip().split('\t')
    experiments = {}
    for row in csv.reader(f, delimiter='\t'):
        gene = row[0]
        for e, v in zip(header[1:], row[1:]):
            experiments.setdefault(e, {})[gene] = v
    return experiments


# Gene deletion script
def gene_deletion(model, mm, solver, algorithm):
    genes = {} # Genes -> reaction1, reaction2
    gene_assoc = {} # Reaction -> gene1, gene2

    # Associate genes with all the reactions
    for reaction in model.parse_reactions():

        # Ignore anything that doesn't have any genes
        if reaction.genes is None:
            continue
        # Express the genes controling a reaction as the Expression datatype
        expr = boolean.Expression(reaction.genes)
        # Associates each reaction (reaction.id) with all the genes that control it (expr)
        gene_assoc[reaction.id] = expr

        # Create a dictionary of each of the genes and all the associated reactions
        for var in expr.variables:
            genes.setdefault(var.symbol, set()).add(reaction.id)

    # Set up the MOMA solver
    p = moma.MOMAProblem(mm, solver)

    # Get the equation that the solver will use for biomass
    biomass = model.get_biomass_reaction()

    # Get the wildtype biomass for comparison later
    wt_biomass = p.get_fba_biomass(biomass)

    # Get the wildtype fluxes that we are going to constrain are model with
    wt_fluxes = p.get_fba_flux(biomass)

    #wt_biomass = p.get_flux(biomass)
    print('Wildtype biomass: {}'.format(wt_biomass))
    for gene, reactions in iteritems(genes):
        #print('Deleting {}...'.format(gene))
        deleted_reactions = set()

        # Search through all the reactions and mark each reaction associated with the gene
        for reaction in reactions:
            e = gene_assoc[reaction].substitute(
                lambda v: False if v.symbol == gene else v)
            if e.has_value() and not e.value:
                deleted_reactions.add(reaction) # Add the reaction to the deleted list

        #print(deleted_reactions)
        constr = []
        # Set all the reactions assiciated to the gene to a flux state of 0
        for reaction in deleted_reactions:
            constr.extend(p.prob.add_linear_constraints(
                p.get_flux_var(reaction) == 0))
        #print(constr)

        # Test the biomass
        try:
            if algorithm == "lp2":
                p.minimize_lp2(biomass, wt_fluxes)
            elif algorithm == "qlp2":
                p.minimize_qlp2(biomass, wt_fluxes)
            elif algorithm == "lp3":
                p.minimize_l1(biomass)
            elif algorithm == "qlp3":
                p.minimize_l2(biomass)
            yield gene, p.get_flux(biomass) / wt_biomass
        except moma.MOMAError:
            yield gene, -5.0
        # Do upon exit
        finally:
            for c in constr:
                c.delete()
