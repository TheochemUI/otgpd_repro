import numpy as np
from scipy.optimize import linear_sum_assignment

cost_matrix = np.array([
    [10, 19, 8,  15],
    [10, 18, 7,  17],
    [13, 16, 9,  14],
    [12, 19, 8,  18]
])

row_ind, col_ind = linear_sum_assignment(cost_matrix)
min_cost = cost_matrix[row_ind, col_ind].sum()
print("--- SciPy Assignment Problem Verification ---")
print("Optimal assignment (row -> col):")
for r, c in zip(row_ind, col_ind):
    print(f"  Row {r} -> Col {c} (Cost: {cost_matrix[r, c]})")

print(f"\nTotal minimum cost: {min_cost}")
