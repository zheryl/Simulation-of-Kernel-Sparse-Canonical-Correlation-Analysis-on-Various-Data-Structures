**Simulation of Kernel Sparse Canonical Correlation Analysis on Various Data Structures**
This project evaluates the performance of Kernel Sparse Canonical Correlation Analysis (Kernel Sparse CCA) compared to Classical CCA and Sparse CCA using simulated datasets with linear, non-linear, and varying sparsity structures.

**Objectives**
1. Assess canonical correlation accuracy across methods
2. Examine stability on test data
3. Evaluate variable selection capability in high-dimensional settings

**Methods**
Three approaches are compared:
1. Classical CCA (linear, no sparsity)
2. Sparse CCA (variable selection with sparsity penalties)
3. Kernel Sparse CCA (captures non-linearity + sparsity)

Synthetic datasets are generated to represent:
1. Linear relationships
2. Non-linear (kernel-based) relationships
3. Different sparsity levels

**Key Findings**
1. Kernel Sparse CCA consistently yields the highest canonical correlation (>0.99).
2. Maintains strong stability on test data (>0.99).
3. Effectively selects relevant variables without weakening canonical relationships.

**Conclusion**
Kernel Sparse CCA outperforms linear methods for complex, high-dimensional, or non-linear data, making it suitable for genomic analysis, medical imaging, and multivariate socioeconomic datasets.
