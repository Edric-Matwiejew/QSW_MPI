import numpy as np
import fortran_sub

print(fortran_sub.sub_2(np.array([1,2]),np.array([1,2])))
print(fortran_sub.sub_2(np.array([]),np.array([])))

print(fortran_sub.sub_3(np.array([1,2]),np.array([1,2]),np.array([1,2],dtype=float)))
print(fortran_sub.sub_3(np.array([]),np.array([]),np.array([],dtype=float)))
