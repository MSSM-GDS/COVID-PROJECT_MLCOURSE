#!/usr/bin/env python

import pandas as pd
import numpy as np
import sys
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler, label_binarize
from sklearn import svm, metrics, model_selection, preprocessing
from sklearn.model_selection import GridSearchCV, train_test_split
from itertools import cycle
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, roc_auc_score
from sklearn.multiclass import OneVsRestClassifier
from scipy import interp



fulldata = pd.read_csv("VOOM_COMBATSEQ_NORM_FILTERED.csv", index_col=0)
#print(fulldata.head(10))

##transpose
fulldata_t = fulldata.T

#print(fulldata_t.head(10))
#print(fulldata_t.columns.values)
#print(fulldata_t.index.values)

##make index into another column
fulldata_t['samples'] = fulldata_t.index

##remove index col
fulldata_t.reset_index(drop=True, inplace=True)
print(fulldata_t.columns.values)

phenodata = pd.read_csv("MASTER_PHENOTYPES_sorted.csv", index_col=0)
print(phenodata.head(10))

###merge phenotype data on sampleid:
inputdata = fulldata_t.merge(phenodata, left_on='samples', right_on='GSEID_PatID')
#print(inputdata.head(10))
#print(inputdata.columns.values)

####create training and test set 80-20
#msk = np.random.rand(len(inputdata)) < 0.8
#train = inputdata[msk]
#print(train.Disease.value_counts())

#test = inputdata[~msk]
#print(test.Disease.value_counts())

####read test and train from csv file
#train, test = train_test_split(inputdata, test_size=0.2)
train = pd.read_csv("FINAL_300/train_data_ML.csv", index_col=0)
test = pd.read_csv("FINAL_300/test_data_ML.csv", index_col=0)

#print(train.head(10))
#print(train.columns.values)


print(train.columns.values)
print(train.phenotype.value_counts())
train_labels = train.phenotype
#train = train.iloc[: , :-7]
train = train.iloc[: , :-1]
print(train.columns.values)
#train_data = train.to_numpy()

print(test.columns.values)
print(test.phenotype.value_counts())
test_labels = test.phenotype
#test = test.iloc[: , :-7]
test = test.iloc[: , :-1]
print(test.columns.values)
#test_data = test.to_numpy()


###Scaling train data
scaler = preprocessing.StandardScaler().fit(train)
X_train_scaled_features = scaler.transform(train)
X_train_scaled = pd.DataFrame(X_train_scaled_features, index=train.index, columns=train.columns)
print(X_train_scaled.shape)
print(X_train_scaled)
train_data = X_train_scaled.to_numpy()

####Scaling test data
X_test_scaled_features = scaler.transform (test)
X_test_scaled = pd.DataFrame(X_test_scaled_features, index=test.index, columns=test.columns)
print(X_test_scaled.shape)
print(X_test_scaled)
test_data = X_test_scaled.to_numpy()


#parameters = [{'kernel':['linear', 'poly', 'rbf', 'sigmoid'],'class_weight':['balanced'], 'max_iter':[1000], 'gamma':[0.0000000001, 0.00000001, 0.00001,0.0001,0.001,0.01,0.1], 'C':[0.1, 1, 10, 100, 1000, 10000], 'probability':[True]}]

parameters = [{'kernel':['rbf'],'class_weight':['balanced'], 'max_iter':[1000], 'gamma':[0.001], 'C':[1000], 'probability':[True]}]

#get the labels in 1D array or list. I don't know what files you'll use so you need to figure this out. Let me know if you have issues
scores = ['precision', 'recall']

for score in scores:
	print("# Tuning hyper-parameters for %s" % score)
	print()
	svc = svm.SVC()
	clf = GridSearchCV(svc, parameters, refit=True, verbose=4, cv=5, scoring='%s_macro' % score)

	clf.fit(train_data, train_labels)
	print("BEST PARAMS:")
	print(clf.best_params_)
	print(clf.best_estimator_)

	print("Best parameters set found on development set:")
	print()
	print(clf.best_params_)
	print()
	print("Grid scores on development set:")
	print()
	means = clf.cv_results_['mean_test_score']
	stds = clf.cv_results_['std_test_score']
	for mean, std, params in zip(means, stds, clf.cv_results_['params']):
		print("%0.3f (+/-%0.03f) for %r" % (mean, std * 2, params))
	print()

	print("Detailed classification report:")
	print()
	print("The model is trained.")
	print("The scores are computed.")
	print()
	prediction2 = clf.predict(test_data)
	print(metrics.classification_report(test_labels, prediction2))
	print()

####Run model after best params
clf2 = make_pipeline(StandardScaler(), svm.SVC(C=1000, class_weight='balanced', gamma=0.01, kernel='rbf', max_iter=1000, probability=True))

clf2.fit(train_data, train_labels)

print("Prediction stats:")
predictions = clf2.predict(test_data)
print(metrics.accuracy_score(test_labels, predictions))
print(metrics.balanced_accuracy_score(test_labels, predictions))
print(metrics.confusion_matrix(test_labels, predictions))
print(metrics.classification_report(test_labels, predictions))
print(metrics.jaccard_score(test_labels, predictions, average = 'weighted'))


#### make labels numeric!
test_labels_numeric = label_binarize(test_labels, classes=['covid', 'viral', 'bacterial', 'healthy', 'others'])
train_labels_numeric = label_binarize(train_labels, classes=['covid', 'viral', 'bacterial', 'healthy', 'others'])
test_labels_uq = ['covid', 'viral', 'bacterial', 'healthy', 'others']
n_classes = len(test_labels_uq)

print("n_classes:")
print(n_classes)
print(train_labels.value_counts())
print(test_labels.value_counts())
print(len(train_labels))

#f = plt.figure()
#metrics.plot_confusion_matrix(clf2, train_labels_numeric, test_labels_numeric)
#plt.show()
#f.savefig("SVM_testCM.pdf", bbox_inches='tight')

#print("area under curve (auc): ", metrics.roc_auc_score(test_labels_numeric, predictions))

# print ROC for each label OVR
clf3 = OneVsRestClassifier(svm.SVC(C=1000, class_weight='balanced', gamma=0.01, kernel='rbf', max_iter=1000, probability=True))
y_score = clf3.fit(train_data, train_labels_numeric).decision_function(test_data)

# Compute ROC curve and ROC area for each class
fpr = dict()
tpr = dict()
roc_auc = dict()
for i in range(n_classes):
    fpr[i], tpr[i], _ = roc_curve(test_labels_numeric[:, i], y_score[:, i])
    roc_auc[test_labels_uq[i]] = auc(fpr[i], tpr[i])
    print("ROC AUC")
    print(test_labels_uq[i])
    print(roc_auc[test_labels_uq[i]])

# Compute micro-average ROC curve and ROC area
fpr["micro"], tpr["micro"], _ = roc_curve(test_labels_numeric.ravel(), y_score.ravel())
roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])
print("ROC MICRO")
print(roc_auc["micro"])


#'''you can add more stats that we need'''
#metrics.plot_precision_recall_curve(clf, test_data, test_labels)
#metrics.plot_roc_curve(clf, test_data, test_labels)

all_fpr = np.unique(np.concatenate([fpr[i] for i in range(n_classes)]))

# Then interpolate all ROC curves at this points
mean_tpr = np.zeros_like(all_fpr)
for i in range(n_classes):
    mean_tpr += interp(all_fpr, fpr[i], tpr[i])

# Finally average it and compute AUC
mean_tpr /= n_classes

fpr["macro"] = all_fpr
tpr["macro"] = mean_tpr
roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])

# Plot all ROC curves
f = plt.figure()
lw = 2
plt.plot(fpr["micro"], tpr["micro"],
         label='micro-average ROC curve (area = {0:0.2f})'
               ''.format(roc_auc["micro"]),
         color='deeppink', linestyle=':', linewidth=4)

plt.plot(fpr["macro"], tpr["macro"],
         label='macro-average ROC curve (area = {0:0.2f})'
               ''.format(roc_auc["macro"]),
         color='navy', linestyle=':', linewidth=4)

colors = cycle(['aqua', 'darkorange', 'cornflowerblue', 'lightgreen', 'purple'])
for i, color in zip(range(n_classes), colors):
    plt.plot(fpr[i], tpr[i], color=color, lw=lw,
             label='ROC curve of class {0} (area = {1:0.2f})'
             ''.format(test_labels_uq[i], roc_auc[test_labels_uq[i]]))

plt.plot([0, 1], [0, 1], 'k--', lw=lw)
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC curves for multinomial(ovr) kernel-based SVM classification')
plt.legend(loc="lower right")
plt.show()
f.savefig("ROCAUC_SVM.pdf", bbox_inches='tight')
