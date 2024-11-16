import numpy as np
from scipy.linalg import hessenberg
from scipy.linalg import qr
from pprint import pprint
# Define the complex matrix A
'''
A = np.array([[5+1j, 9+6j, 3+2j, 5+1j, 1+1j],
              [3+2j, 3+7j, 0+4j, 2+3j, 2+4j],
              [6+3j, 2+8j, 8+6j, 9+5j, 3+7j],
              [8+4j, 7+9j, 1+8j, 8+7j, 4+2j],
              [5+5j, 2+0j, 6+0j, 7+9j, 5+5j]])
A = np.array([[5+1j, 9+6j, 3+9j, 5+8j],
              [3+2j, 3+7j, 0+7j, 2+6j],
              [6+3j, 2+8j, 8+5j, 9+4j],
              [8+4j, 7+9j, 1+3j, 8+2j]])

A = np.array([[5+0j, 9+0j, 3+0j, 5+0j],
              [3+0j, 3+0j, 0+0j, 2+0j],
              [6+0j, 2+0j, 8+0j, 9+0j],
              [8+0j, 7+0j, 1+0j, 8+0j]])
'''
A = np.array([
    [complex(2, 0), complex(3, 0), complex(5, 0), complex(7, 0), complex(11, 0)],
    [complex(13, 0), complex(17, 0), complex(19, 0), complex(23, 0), complex(29, 0)],
    [complex(31, 0), complex(37, 0), complex(41, 0), complex(43, 0), complex(47, 0)],
    [complex(53, 0), complex(59, 0), complex(61, 0), complex(67, 0), complex(71, 0)],
    [complex(73, 0), complex(79, 0), complex(83, 0), complex(89, 0), complex(97, 0)]
])
# Compute the Hessenberg form using SciPy's 'hessenberg' function
H = hessenberg(A)
#print("Upper hessenberg:")
#pprint(H)
eig1, _= np.linalg.eig(H)
#print("Eigen values of upper hessenberg matrix")
pprint(eig1)
Q, R = qr(H)
#print("Upper triangular:")
pprint(R)
#eig2, _ = np.linalg.eig(A)
#print("Eigen values of original matrix")
pprint(eig2)
print("\n")



