#!/usr/bin/env python3
import sys
import json
from subprocess import call
from cameo import models
from cameo.strain_design import pathway_prediction

RESULTS = []
n_pathways = int(sys.argv[3])
n_progress = 0


def reaction_to_dict(reaction):
    return dict(
        id=reaction.id,
        name=reaction.name,
        reaction_string=reaction.build_reaction_string(use_metabolite_names=True),
    )


def save_results(pathway):
    global n_progress
    RESULTS.append([reaction_to_dict(reaction) for reaction in pathway.reactions])
    n_progress += 1
    print('{"proc.progress":%.2f}' % (n_progress / n_pathways,))
    call(['re-save', 'pathways', json.dumps(RESULTS)])


model = getattr(models.bigg, sys.argv[1])
predictor = pathway_prediction.PathwayPredictor(model)
pathways = predictor.run(product=sys.argv[2], max_predictions=n_pathways, callback=save_results)
