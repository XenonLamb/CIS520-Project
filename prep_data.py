from rdkit import Chem
import random
import math


def gen_split(datapath, trainrate=0.8, devrate=0.9):
    suppls = Chem.SDMolSupplier(datapath, removeHs = False, sanitize = False)
    suppl = [x for x in suppls]
    nummols = len(suppl)
    numtrain = int(math.floor(trainrate * nummols))
    numdev = int(math.floor(devrate*nummols))
    idxs = list(range(0,nummols))
    random.shuffle(idxs)
    trainadj = []
    trainnames = []
    trainlabels = []
    #trainsmiles = []
    devadj = []
    devnames = []
    devlabels = []
    #devsmiles = []
    testadj = []
    testnames = []
    testlabels = []
    #testsmiles = []
    for i in range(0,numtrain):
        idx = idxs[i]
        if suppl[idx]==None:
          print("No molecule at", idx)
        else:
          trainadj.append(Chem.rdmolops.GetAdjacencyMatrix(suppl[idx], True))
          trainnames.append([ato.GetAtomicNum() for ato in suppl[idx].GetAtoms()])
          trainlabels.append(int(float(suppl[idx].GetProp("value"))))
          #trainsmiles.append(Chem.MolToSmiles(suppl[idx]))
    for i in range(numtrain,numdev):
        idx = idxs[i]
        devadj.append(Chem.rdmolops.GetAdjacencyMatrix(suppl[idx], True))
        devnames.append([ato.GetAtomicNum() for ato in suppl[idx].GetAtoms()])
        devlabels.append(int(float(suppl[idx].GetProp("value"))))
        #devsmiles.append(Chem.MolToSmiles(suppl[idx]))
    for i in range(numdev,nummols):
        idx = idxs[i]
        testadj.append(Chem.rdmolops.GetAdjacencyMatrix(suppl[idx], True))
        testnames.append([ato.GetAtomicNum() for ato in suppl[idx].GetAtoms()])
        testlabels.append(int(float(suppl[idx].GetProp("value"))))
        #testsmiles.append(Chem.MolToSmiles(suppl[idx]))

    return trainadj, trainnames, trainlabels, devadj, devnames, devlabels, testadj,\
           testnames,  testlabels


def nameArray2Dic(names):
    """
    input
        names: 1d array containing atom names of 1 molecule
    return
        namedic: dictionary in the form of "node index":"atom name"
    """
    namedic = {}
    for idx in range(0,len(names)):
        namedic[idx] = names[idx]

    return namedic


def padded_spectral(moladj, embedding_dimension=30, normalized=True):
    # Padding with zeros
    embedding = np.zeros(embedding_dimension)

    # Usage of networkx graph objects
    graph = nx.from_numpy_matrix(np.array(moladj))
    adj_matrix = nx.adj_matrix(graph)
    n_nodes, m_nodes = adj_matrix.shape
    k = min(embedding_dimension + 1, n_nodes - 1)

    if normalized:
        laplacian = nx.normalized_laplacian_matrix(graph)
    else:
        laplacian = nx.laplacian_matrix(graph)

    # Minus the eigen decomposition of minus the Laplacian is more stable than directly
    # computing the eigen decomposition of the Laplacian

    v0 = np.random.uniform(-1, 1, laplacian.shape[0])
    eigenvalues = sparse.linalg.eigsh(-laplacian, k=k, sigma=1.0, which='LM', tol=1e-6, v0=v0,
                                      return_eigenvectors=False)
    embedding[:len(eigenvalues) - 1] = sorted(-eigenvalues)[1:]

    return embedding


def gen_emb(listadj):
    listemb = []
    for i in range(len(listadj)):
        listemb.append(padded_spectral(listadj[i]))
    return listemb
