""" Functions used to save/load data contained in Python Dicts.
"""
import pickle

def saveData(data, filename):
    " Save a dict containing fields of data in a pickle file. "
    with open(filename + '.pickle', 'wb') as f:
        pickle.dump(data, f, 2)

def loadData(filename):
    " Load a dict containing fields of data in a pickle file. "
    with open(filename + '.pickle', 'rb') as f:
        return pickle.load(f, encoding='latin1')