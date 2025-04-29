import pickle
from pdb import set_trace as bp

def read_pkl(pkl_file):
    with open(pkl_file, 'rb') as f:
        return pickle.load(f, encoding='latin1')

for i in range(4):
    data = read_pkl(f'/eos/user/j/jiehan/finalfit_102X/CMSSW_10_2_13/src/flashggFinalFit/Datacard/yields_fiducial/VBF{i}.pkl')
    print(data['nominal_yield'][:-2].sum())