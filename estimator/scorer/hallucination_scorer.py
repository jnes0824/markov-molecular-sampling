import os
import pickle
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator

from ...common.chem import fingerprints_from_mol

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)


ROOT_DIR = 'MARS/estimator/scorer'
TASKS = ['hallucination']

models = {}

def mol_to_fingerprint(mol):
    if mol is None:
        return None
    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fingerprint = mfpgen.GetFingerprint(mol)
    return np.array(fingerprint)

def load_model(task):
    with open(os.path.join(ROOT_DIR, 'custom_qsar/%s.pkl' % task), 'rb') as f:
        models[task] = pickle.load(f, encoding='iso-8859-1')

def get_scores(task, mols):
    model = models.get(task)
    if model is None:
        load_model(task)
        model = models[task]
        
    fps = [mol_to_fingerprint(mol) for mol in mols]
    fps = np.stack(fps, axis=0)
    scores = models[task].predict_proba(fps)
    hallucination_probabilities = scores[:, 1]
    non_hallucination_probabilities = 1.0 - hallucination_probabilities
    return non_hallucination_probabilities.tolist()
