# CIS520-Project
Graph classification on NCI

classes: 2 
maximum node tag: 37
data: 4110

https://ls11-www.cs.tu-dortmund.de/staff/morris/graphkerneldatasets


GNN 10 folds: 
Parameter containing:
tensor([0., 0., 0., 0.], device='cuda:0', requires_grad=True)
epoch: 350: 100%|██████████| 50/50 [00:00<00:00, 76.89batch/s]
loss training: 0.341898
accuracy train: 0.880746 test: 0.808252
precision test: 0.803828
recall test: 0.815534
fscore test: 0.809639

SVM:
Accuracy: 0.6155988857938719,Precision: 0.6641221374045801, recall: 0.48066298342541436, F score: 0.5576923076923076
Log regression:
Accuracy: 0.6072423398328691,Precision: 0.6298701298701299, recall: 0.5359116022099447, F score: 0.5791044776119402
Adaboost:
Accuracy: 0.6573816155988857,Precision: 0.6959459459459459, recall: 0.569060773480663, F score: 0.6261398176291793
Random forest:
Accuracy: 0.7465181058495822,Precision: 0.7960526315789473, recall: 0.6685082872928176, F score: 0.7267267267267267

accuracy: {'RandomWalk': 0.5, 'ShortestPath': 0.62}
F1:{'RandomWalk': 0.46, 'ShortestPath': 0.6}
