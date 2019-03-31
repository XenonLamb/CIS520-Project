from rdkit import Chem
import random

def gen_split(datapath):
    suppl = Chem.SDMolSupplier(datapath)
    nummols = len(suppl)
    numtrain = int(floor(0.8 * nummols))
    numdev = int(floor(0.9*nummols))
    idxs = list(range(0,nummols))
    random.shuffle(idxs)
    trainadj = []
    trainnames = []
    trainlabels = []
    trainsmiles = []
    devadj = []
    devnames = []
    devlabels = []
    devsmiles = []
    testadj = []
    testnames = []
    testlabels = []
    testsmiles = []
    for i in range(0,numtrain):
        idx = idxs[i]
        trainadj.append(Chem.rdmolops.GetAdjacencyMatrix(suppl[idxs[idx]], True))
        trainnames.append([ato.GetAtomicNum() for ato in suppl[idxs[idx]].GetAtoms()])
        trainlabels.append(suppl[idxs[idx]].GetProp("value"))
        trainsmiles.append(Chem.MolToSmiles(suppl[idxs[idx]]))
    for i in range(numtrain,numdev):
        idx = idxs[i]
        devadj.append(Chem.rdmolops.GetAdjacencyMatrix(suppl[idxs[idx]], True))
        devnames.append([ato.GetAtomicNum() for ato in suppl[idxs[idx]].GetAtoms()])
        devlabels.append(suppl[idxs[idx]].GetProp("value"))
        devsmiles.append(Chem.MolToSmiles(suppl[idxs[idx]]))
    for i in range(numdev,nummols):
        idx = idxs[i]
        testadj.append(Chem.rdmolops.GetAdjacencyMatrix(suppl[idxs[idx]], True))
        testnames.append([ato.GetAtomicNum() for ato in suppl[idxs[idx]].GetAtoms()])
        testlabels.append(suppl[idxs[idx]].GetProp("value"))
        testsmiles.append(Chem.MolToSmiles(suppl[idxs[idx]]))

    return trainadj, trainnames, trainsmiles, trainlabels, devadj, devnames, devsmiles, devlabels, testadj,\
           testnames, testsmiles, testlabels


