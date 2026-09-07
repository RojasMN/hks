# Hierarchical Kolmogorov-Smirnov (HKS) Test

Python and R code designed for performing Hierarchical (Nested) Kolmogorov-Smirnov (KS) and Anderson-Darling (AD) tests. It addresses the critical statistical challenge of pseudoreplication in grouped or nested datasets (e.g., multiple neurons recorded from the same animal).

---

## The Problem: Pseudoreplication in Nested Data

Standard statistical tests assume complete independence across all data points. In biological and experimental settings, data is frequently hierarchical—such as recording multiple cells across a handful of subjects. 

Pooling all measurements treats dependent observations as independent samples, artificially inflating the sample size. This results in severe type I error inflation, leading to false positives where significant differences are detected solely due to correlated within-subject measurements.

---

## The Solution: Hierarchical Permutation Testing

The package implements the `HierarchicalPermutationTest` class to compute valid p-values that respect nested data structures through three core mechanisms:

* **Robust Observed Statistics:** To handle highly unbalanced designs (varying numbers of subunits and observations per subject), the algorithm draws balanced subsamples and calculates the median KS or Anderson-Darling statistic, preventing outlier dominance.
* **Hierarchical Shuffling:** Rather than shuffling individual observations, it permutes top-level group labels (e.g., Subject IDs), preserving the natural correlation structure within each subject.
* **Accurate P-Values:** The null distribution is built through thousands of hierarchical permutations, ensuring the resulting p-value reflects the true number of independent subjects rather than total observations.

---

## Data Structure Requirements

Input data must be supplied as a long-format table (CSV or Pandas DataFrame) containing four required columns:

| Column | Description | Example |
| :--- | :--- | :--- |
| **`subject`** | Top-level independent unit | Animal ID, Patient, Batch |
| **`subunit`** | Nested unit within the subject | Cell, Region of Interest |
| **`group`** | Experimental condition or treatment category | Control, Treatment |
| **`variable`** | Continuous numerical metric under analysis | Measurement value |

---

## Implementation & Language Support

* **Python (Core):** Highly optimized using NumPy vectorization to bypass DataFrame overhead, making it ideal for large datasets and intensive bootstrapping.
* **R:** An equivalent implementation of the statistical procedures is located in the `/R` directory for R-based workflows.
