from sys import argv
import numpy as np
import pandas as pd
import MatrixOperations as MO

A = np.zeros((4,6))
K = np.zeros((4,6))
print A
n = np.zeros((6,))
k = 0
for j in xrange(1,7):
    name = "results_N="+str(j)
    i = 0
    with open(name) as f:
        n[j-1] = float(f.readline())
        for line in f:
            if line[0:len('assemble_system')] == 'assemble_system':
                A[i,k] = (float(line[-22:-12]))
            if line[0:len('ksp.solve')] == 'ksp.solve':
                K[i,k] = (float(line[-22:-12]))
                i = i+1
    k = k+1
        # print A

print n
print A
print K

print "\n\n"

AssembleTable = n#[1, 2, 3, 4, 5, 6]
# Assemble = np.concatenate((n.T,A), axis=1)
LatexTable = pd.DataFrame(A, columns = AssembleTable)
pd.set_option('precision',3)
LatexTable = MO.PandasFormat(LatexTable,n[1-1],"%1.2e")
LatexTable = MO.PandasFormat(LatexTable,n[2-1],"%1.2e")
LatexTable = MO.PandasFormat(LatexTable,n[3-1],"%1.2e")
LatexTable = MO.PandasFormat(LatexTable,n[4-1],"%1.2e")
LatexTable = MO.PandasFormat(LatexTable,n[5-1],"%1.2e")
LatexTable = MO.PandasFormat(LatexTable,n[6-1],"%1.2e")
print LatexTable.to_latex()

AssembleTable = n#[ "1", "2", "3", "4", "5", "6"]
# Assemble = np.concatenate((n,A), axis=2)
LatexTable = pd.DataFrame(K, columns = AssembleTable)
pd.set_option('precision',3)
LatexTable = MO.PandasFormat(LatexTable,n[1-1],"%1.2e")
LatexTable = MO.PandasFormat(LatexTable,n[2-1],"%1.2e")
LatexTable = MO.PandasFormat(LatexTable,n[3-1],"%1.2e")
LatexTable = MO.PandasFormat(LatexTable,n[4-1],"%1.2e")
LatexTable = MO.PandasFormat(LatexTable,n[5-1],"%1.2e")
LatexTable = MO.PandasFormat(LatexTable,n[6-1],"%1.2e")
print LatexTable.to_latex()