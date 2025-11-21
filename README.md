Simulation of Kernel Sparse Canonical Correlation Analysis on Various Data Structures

This project evaluates the performance of Kernel Sparse Canonical Correlation Analysis (Kernel Sparse CCA) compared to Classical CCA and Sparse CCA using simulated datasets with linear, non-linear, and varying sparsity structures.

Objectives

Assess canonical correlation accuracy across methods

Examine stability on test data

Evaluate variable selection capability in high-dimensional settings

Methods

Three approaches are compared:

Classical CCA (linear, no sparsity)

Sparse CCA (variable selection with sparsity penalties)

Kernel Sparse CCA (captures non-linearity + sparsity)

Synthetic datasets are generated to represent:

Linear relationships

Non-linear (kernel-based) relationships

Different sparsity levels

Key Findings

Kernel Sparse CCA consistently yields the highest canonical correlation (>0.99).

Maintains strong stability on test data (>0.99).

Effectively selects relevant variables without weakening canonical relationships.

Conclusion

Kernel Sparse CCA outperforms linear methods for complex, high-dimensional, or non-linear data, making it suitable for genomic analysis, medical imaging, and multivariate socioeconomic datasets.
