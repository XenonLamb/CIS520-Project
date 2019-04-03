"""
SVM implementation
"""

import numpy as np
#
from grakel import Graph, RandomWalk, ShortestPath, WeisfeilerLehman, MultiscaleLaplacian, Propagation
from sklearn.metrics import accuracy_score
from sklearn.svm import SVC

from prep_data import gen_split


DATA_ROOT = "../NCI_FULL/109total-connect.sdf"
NUM_TRAIN_INSTANCES = 100
NUM_TEST_INSTANCES = 50
KERNELS = [RandomWalk, ShortestPath, WeisfeilerLehman,
           MultiscaleLaplacian, Propagation]


def main():
    train_adj, train_names, train_labels, devadj, devnames, devlabels, test_adj, test_names, test_labels = gen_split(
        DATA_ROOT)

    num_train_samples = len(train_names)
    num_test_samples = len(test_names)

    print("Number of training instances: %d" % num_train_samples)
    print("Number of test instances: %d" % num_test_samples)

    training_data = [train_adj, train_names, train_labels]
    [train_adj, train_names, train_labels] = [ls[:NUM_TRAIN_INSTANCES]
                                              for ls in training_data]

    test_data = [test_adj, test_names, test_labels]
    [test_adj, test_names, test_labels] = [ls[:NUM_TEST_INSTANCES]
                                           for ls in test_data]

    train_graphs = create_graphs(
        train_adj, train_names)
    test_graphs = create_graphs(
        test_adj, test_names)

    results = {}
    for graph_kernel in KERNELS:

        gk = graph_kernel()
        kernel_name = type(gk).__name__
        print("Kernel is %s" % kernel_name)

        acc = train_SVM(gk, train_graphs, test_graphs,
                        train_labels, test_labels, kernel=graph_kernel)
        results[kernel_name] = acc

        print("-" * 50 + '\n')

    print(results)


def create_graphs(adj_mat, names):
    """Create list of Graph instances
    params:
        adj_mat: list of matrices       - adjacency matrices of graphs
        names: list of lists of ints    - names of each node (molecule) in graph

    returns:
        list of Graph instances
    """

    num_samples = len(names)
    return [Graph(adj_mat[i], names[i], graph_format="adjacency") for i in range(num_samples)]


def train_SVM(gk, train_graphs, test_graphs, train_labels, test_labels, kernel=RandomWalk):
    """Train SVM model for given kernel on NCI data
    params:
        gk: grakel.Kernel instance              - SVM graph kernel
        train_graphs: list of Graph instances   - training data
        test_graphs: list of Graph instances    - test data
        train_labels: list of ints              - training labels
        kernel: grakel.Kernel function          - kernel for SVM

    returns:
        accuracy of kernel SVM
    """

    print("Computing training kernel")
    K_train = gk.fit_transform(train_graphs)

    print("Computing test kernel")
    K_test = gk.transform(test_graphs)

    print("Classifying molecules")
    clf = SVC(kernel='precomputed')

    # Fit on the train Kernel
    clf.fit(K_train, train_labels)

    # Predict and test
    y_pred = clf.predict(K_test)

    return accuracy_score(test_labels, y_pred)


if __name__ == '__main__':
    main()
