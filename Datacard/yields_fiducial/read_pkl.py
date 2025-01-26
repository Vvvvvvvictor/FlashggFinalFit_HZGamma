import pickle
from pdb import set_trace as bp

def read_pkl(pkl_file):
    with open(pkl_file, 'rb') as f:
        return pickle.load(f, encoding='latin1')

data = read_pkl('/eos/user/j/jiehan/finalfit_102X/CMSSW_10_2_13/src/flashggFinalFit/Datacard/yields_fiducial/VBF3.pkl')
bp()