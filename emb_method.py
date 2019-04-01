from sklearn.ensemble import AdaBoostClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC


def get_precision(y_pred, y_true):
    ## YOUR CODE HERE...
    counts = Counter(zip(y_pred, y_true))
    TP = counts[1, 1]
    FP = counts[1, 0]
    if (TP + FP > 0):
        precision = TP / float(TP + FP)
    else:
        precision = 0
    return precision


## Calculates the recall of the predicted labels
def get_recall(y_pred, y_true):
    ## YOUR CODE HERE...
    counts = Counter(zip(y_pred, y_true))
    TP = counts[1, 1]
    FN = counts[0, 1]
    if (TP + FN > 0):
        recall = TP / float(TP + FN)
    else:
        recall = 0
    return recall


## Calculates the f-score of the predicted labels
def get_fscore(y_pred, y_true):
    ## YOUR CODE HERE...
    precision = get_precision(y_pred, y_true)
    recall = get_recall(y_pred, y_true)
    if (precision + recall > 0):
        fscore = 2 * precision * recall / (precision + recall)
    else:
        fscore = 0
    return fscore


def test_predictions(y_pred, y_true):
    precision = get_precision(y_pred, y_true)
    recall = get_recall(y_pred, y_true)
    fscore = get_fscore(y_pred, y_true)

    print("Precision: {}, recall: {}, F score: {}".format(precision, recall, fscore))

    return

def eval(model, trainemb, trainlabel, valemb, vallabel):
    model.fit(trainemb, trainlabel)
    pred = model.predict(valemb)
    test_predictions(pred, vallabel)

    return
