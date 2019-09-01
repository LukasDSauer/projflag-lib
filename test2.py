#  Developed by Lukas Sauer at the Heidelberg Institute of Theoretical Studies on 17.05.19, 15:30.
#  Contact: lukas.sauer@h-its.org.
#  Last modified on 17.05.19, 15:30.
#  (C) 2019. All rights reserved.

from LinAlgebraUtility import get_adj
import numpy as np

A = np.array([[1,9,3],[4,5,6],[7,8,9]])
print(np.linalg.det(A))
B = get_adj(A)
print(B)
C = np.matmul(A, B)
print(C)
