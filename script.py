import numpy as np
from scipy.linalg import hessenberg

# Define the complex matrix A
A = np.array([[4+0j, 1+0j, 3+0j, 4+0j],
              [1+0j, 2+0j, 4+0j, 0+0j],
              [3+0j, 4+0j, 5+0j, 1+0j],
              [4+0j, 0+0j, 1+0j, 4+0j]])

# Compute the Hessenberg form using SciPy's 'hessenberg' function
H = hessenberg(A)
eigenvalues_np = np.linalg.eigvals(H)
print("Eigenvalues using NumPy:", eigenvalues_np)
eigenvalues_np_1 = np.linalg.eigvals(A)
print("Eigenvalues using NumPy:", eigenvalues_np_1)

print(H)
