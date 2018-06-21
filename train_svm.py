from sklearn import svm
import numpy as np
import pickle
from sklearn.cross_validation import train_test_split
def loadDataSet(filename):
    dataMat = []
    labelMat = []
    with open(filename) as fr:
        for line in fr.readlines():
            lineArr = line.strip().split('\t')
            dataMat.append([float(lineArr[0]), float(lineArr[1]), float(lineArr[2]),float(lineArr[3]),float(lineArr[4]), float(lineArr[5]), float(lineArr[6]), float(lineArr[7]), float(lineArr[8]), float(lineArr[9]), float(lineArr[10]), float(lineArr[11]), float(lineArr[12]), float(lineArr[13]), float(lineArr[14]), float(lineArr[15])])
            if float(lineArr[-1]) == 1.0:
                lineArr[-1] =0
            labelMat.append(float(lineArr[-1]))
    print labelMat
    print "_______________________________________________"
    return dataMat, labelMat


if __name__ == '__main__':
    x, y = loadDataSet("./data")
    X_train, X_test, y_train, y_test = train_test_split(x, y, test_size=0.25)
    #clf=svm.SVC()
    clf=svm.SVC(gamma=0.001,C=1000)
    clf.fit(X_train, y_train)
    print "++++++++++++++++++++++++++++++++++++"
    print clf
    res=clf.predict(X_test)
    ans=clf.predict(X_train)
    print res
    len_arr = len(res)
    l = 0
    #for i in range(0, len_arr):
    #    print res[i], y_test
    for i in range(0,len_arr):
        if y_test[i] == res[i] :
            l = l+ 1
    ans_len = len(ans)
    m = 0
    for i in range(0,ans_len):
        if y_train[i] == ans[i]:
            m = m+1;
    print l
    print len_arr
    print float(l)/len_arr
    print float(m)/ans_len

    #save model
    with open('model.pickle', 'wb') as f:
        pickle.dump(clf, f)

