from numpy import *

def run(size=4,filename="input.txt"):
    A = random.rand(size,size)
    (Q,R) = linalg.qr(A)
    for i in xrange(size):
        for j in xrange(i+1,size):
            R[i][j] = 0
        R[i][i] = abs(R[i][i])
    QT = transpose(Q)
    A = dot(dot(QT,R),Q)
    X = random.rand(size)
    B = dot(A,X)
    outFile = open(filename,"w")
    for i in xrange(size):
        for j in xrange(size):
            outFile.write("%lf " %(A[i][j]))
        outFile.write("\n")
    for i in xrange(size):
        outFile.write("%lf " %(B[i]))
    outFile.write("\n")
    for i in xrange(size):
        outFile.write("%lf " %(X[i]))
    outFile.write("\n")
    outFile.close()
