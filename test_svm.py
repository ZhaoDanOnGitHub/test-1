from sklearn import svm
import numpy as np
import pickle
from sklearn.cross_validation import train_test_split
def loadDataSet(filename):
    dataMat = []
    with open(filename) as fr:
        for line in fr.readlines():
            lineArr = line.strip().split('\t')
            dataMat.append([float(lineArr[0]), float(lineArr[1]), float(lineArr[2]),float(lineArr[3]),float(lineArr[4]), float(lineArr[5]), float(lineArr[6]), float(lineArr[7]), float(lineArr[8]), float(lineArr[9]), float(lineArr[10]), float(lineArr[11]), float(lineArr[12]), float(lineArr[13]), float(lineArr[14]), float(lineArr[15])])
    print "_______________________________________________"
    return dataMat


if __name__ == '__main__':
    X_test= loadDataSet("./data")
    clf=svm.SVC(gamma=0.001,C=1000)

    #load model
    with open('model.pickle', 'rb') as f:
        model = pickle.load(f)
        pre = model.predict(X_test)
    print pre

    print "++++++++++++++++++++++++++++++++++++"

