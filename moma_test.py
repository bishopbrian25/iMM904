import csv

from six import iteritems

from psamm.expression import boolean
from psamm import fluxanalysis


def read_experiments(f):
    header = f.readline().strip().split('\t')
    experiments = {}

    for row in csv.reader(f, delimiter='\t'):
        gene = row[0]
        for e, v in zip(header[1:], row[1:]):
            experiments.setdefault(e, {})[gene] = v

    return experiments


def gene_deletion(model, mm, solver):
    genes = {}
    gene_assoc = {}
    for reaction in model.parse_reactions():
        if reaction.genes is None:
            continue
        expr = boolean.Expression(reaction.genes)
        gene_assoc[reaction.id] = expr
        for var in expr.variables:
            genes.setdefault(var.symbol, set()).add(reaction.id)

    p = fluxanalysis.FluxBalanceProblem(mm, solver)
    biomass = model.get_biomass_reaction()

    p.maximize(biomass)
    wt_biomass = p.get_flux(biomass)
    print('Wildtype biomass: {}'.format(wt_biomass))

    for gene, reactions in iteritems(genes):
        #print('Deleting {}...'.format(gene))

        deleted_reactions = set()
        for reaction in reactions:
            e = gene_assoc[reaction].substitute(
                lambda v: False if v.symbol == gene else v)
            if e.has_value() and not e.value:
                deleted_reactions.add(reaction)

        constr = []
        for reaction in deleted_reactions:
            constr.extend(p.prob.add_linear_constraints(
                p.get_flux_var(reaction) == 0))

        try:
            p.maximize(biomass)
            yield gene, p.get_flux(biomass) / wt_biomass
        except fluxanalysis.FluxBalanceError:
            yield gene, 0.0
        finally:
            for c in constr:
                c.delete()
