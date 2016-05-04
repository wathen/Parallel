from sys import argv
import numpy as np
import pandas as pd
import MatrixOperations as MO
m = 7
A = np.zeros((6,m))
K = np.zeros((6,m))
print A
n = np.zeros((m,))
k = 0
for j in xrange(1,m+1):
    name = "results_N="+str(j)
    print name
    i = 0
    with open(name) as f:
        n[j-1] = float(f.readline())
        for line in f:
            if line[0:len('assemble_system')] == 'assemble_system':
                A[i,k] = (float(line[-22:-12]))
            if line[0:len('ksp.solve')] == 'ksp.solve':
                K[i,k] = (float(line[-22:-12]))
                i = i+1
            if line[0:len('=   BAD TERMINATION')] == '=   BAD TERMINATION':
                i = i+1
    k = k+1
        # print A

print n
print A
print A*8
print K

print "\n\n"

AssembleTable = n#[1, 2, 3, 4, 5, 6]
# Assemble = np.concatenate((n.T,A), axis=1)
LatexTable = pd.DataFrame(A, columns = AssembleTable)
pd.set_option('precision',3)
LatexTable = MO.PandasFormat(LatexTable,n[1-1],"%4.2f")
LatexTable = MO.PandasFormat(LatexTable,n[2-1],"%4.2f")
LatexTable = MO.PandasFormat(LatexTable,n[3-1],"%4.2f")
LatexTable = MO.PandasFormat(LatexTable,n[4-1],"%4.2f")
LatexTable = MO.PandasFormat(LatexTable,n[5-1],"%4.2f")
LatexTable = MO.PandasFormat(LatexTable,n[6-1],"%4.2f")
LatexTable = MO.PandasFormat(LatexTable,n[7-1],"%4.2f")
print LatexTable.to_latex()

AssembleTable = n#[1, 2, 3, 4, 5, 6]
# Assemble = np.concatenate((n.T,A), axis=1)
LatexTable = pd.DataFrame(8*A, columns = AssembleTable)
pd.set_option('precision',3)
LatexTable = MO.PandasFormat(LatexTable,n[1-1],"%4.2f")
LatexTable = MO.PandasFormat(LatexTable,n[2-1],"%4.2f")
LatexTable = MO.PandasFormat(LatexTable,n[3-1],"%4.2f")
LatexTable = MO.PandasFormat(LatexTable,n[4-1],"%4.2f")
LatexTable = MO.PandasFormat(LatexTable,n[5-1],"%4.2f")
LatexTable = MO.PandasFormat(LatexTable,n[6-1],"%4.2f")
LatexTable = MO.PandasFormat(LatexTable,n[7-1],"%4.2f")
print LatexTable.to_latex()

AssembleTable = n#[ "1", "2", "3", "4", "5", "6"]
# Assemble = np.concatenate((n,A), axis=2)
LatexTable = pd.DataFrame(K, columns = AssembleTable)
pd.set_option('precision',3)
LatexTable = MO.PandasFormat(LatexTable,n[1-1],"%4.2f")
LatexTable = MO.PandasFormat(LatexTable,n[2-1],"%4.2f")
LatexTable = MO.PandasFormat(LatexTable,n[3-1],"%4.2f")
LatexTable = MO.PandasFormat(LatexTable,n[4-1],"%4.2f")
LatexTable = MO.PandasFormat(LatexTable,n[5-1],"%4.2f")
LatexTable = MO.PandasFormat(LatexTable,n[6-1],"%4.2f")
LatexTable = MO.PandasFormat(LatexTable,n[7-1],"%4.2f")

print LatexTable.to_latex()


AssembleTable = n#[ "1", "2", "3", "4", "5", "6"]
# Assemble = np.concatenate((n,A), axis=2)
LatexTable = pd.DataFrame(8*K, columns = AssembleTable)
pd.set_option('precision',3)
LatexTable = MO.PandasFormat(LatexTable,n[1-1],"%4.2f")
LatexTable = MO.PandasFormat(LatexTable,n[2-1],"%4.2f")
LatexTable = MO.PandasFormat(LatexTable,n[3-1],"%4.2f")
LatexTable = MO.PandasFormat(LatexTable,n[4-1],"%4.2f")
LatexTable = MO.PandasFormat(LatexTable,n[5-1],"%4.2f")
LatexTable = MO.PandasFormat(LatexTable,n[6-1],"%4.2f")
LatexTable = MO.PandasFormat(LatexTable,n[7-1],"%4.2f")

print LatexTable.to_latex()