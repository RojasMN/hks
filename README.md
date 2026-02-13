# Hierarchical Kolmogorov-Smirnov (HKS) Test
A Python package for performing a Hierarchical (Nested) Kolmogorov-Smirnov Test. This tool is designed for datasets where observations are not independent but are grouped (e.g., multiple neurons recorded from the same animal), addressing the problem of pseudoreplication.

* **Note on Language Support:** The core of this package is implemented in Python. For users who prefer working in R, the `/R` directory contains an equivalent implementation of the statistical procedures. However, Python implementation is significantly faster. It is optimized using NumPy vectorization to bypass the overhead of DataFrame processing, making it the preferred choice for large datasets or high-iteration bootstraps.
  
## 1. Conceptual Background

### The Kolmogorov-Smirnov (KS) Test
The **Kolmogorov-Smirnov test** is a non-parametric method used to determine if two datasets differ significantly in their underlying distributions. Conceptually, it works by plotting the **Empirical Cumulative Distribution Function (ECDF)** for both datasets. The test statistic ($D$) is simply the **maximum vertical distance** between these two curves. If $D$ is large, the distributions are likely different.

<p align="center">
  <img src="images/d_statistic.png" width="500" title="KS Test - D Statistic">
</p>

### The Anderson-Darling (AD) Test
The **Anderson-Darling test** is a modification of the KS test. While the KS test looks only at the *maximum* distance (which often happens near the median/center of the distribution), the Anderson-Darling test uses a weighted distance that places more emphasis on the **tails** (the extremes) of the distribution. It is often more sensitive to differences in the start or end of the distributions.

## 2. Data Structure

This package is designed to support hierarchical (nested) data structures. It is specifically optimized for **unbalanced designs**, where:
* Each **Subject** can have a different number of **Subunits**.
* Each **Subunit** can have a different number of **Observations** (measurements).

### Input Format
Your input data should be a long-format table (CSV or Pandas DataFrame) with the following structure:

| subject | subunit | group | variable |
|:---|:---|:---|:---|
| S_B0 | CELL_S_B0_0 | B | 2.807 |
| S_B0 | CELL_S_B0_0 | B | 1.350 |
| S_B0 | CELL_S_B0_1 | B | 1.612 |
| S_T1 | CELL_S_T1_0 | T | 0.896 |

### Hierarchy Overview
- **Subject (`subject`):** The top-level experimental unit (e.g., Animal ID, Patient, or Batch).
- **Subunit (`subunit`):** The nested unit within the subject (e.g., individual Cells, Regions of Interest, or Sensors).
- **Group (`group`):** The experimental condition, genotype, or treatment group.
- **Variable (`variable`):** The continuous numerical value being analyzed.

> [!NOTE]
> The statistical procedures in this package account for the nested nature of the data. In the provided `data_example.csv`, subjects contain between 10 and 56 subunits, and each subunit contains between 25 and 200 individual observations, demonstrating the package's ability to handle highly unbalanced datasets.

## 3.. The Problem: Pseudoreplication in Nested Data

Standard statistical tests (like the standard KS test) assume that every data point is **independent**. However, in biological and scientific data, this is often false.

**Example:**
* You have 2 groups of mice (Control vs. Treatment).
* You record 100 neurons from *each* mouse.
* You pool all neurons and compare 500 control neurons vs. 500 treatment neurons.
* The test treats this as $N = 1000$ independent samples. In reality, you only have $N=10$ mice. The neurons within a single mouse are highly correlated.

This inflation of the effective sample size is called **pseudoreplication**. It causes standard tests to return artificially low p-values, leading you to detect significant differences that don't actually exist (False Positives).

## 4. The Solution: `HierarchicalPermutationTest` Class

The `HierarchicalPermutationTest` class implements a **hierarchical resampling (permutation) test** to calculate a valid p-value that respects the grouped structure of the data.

### How it works (`hks.py`)

1.  **Robust Observed Statistic & Handling Imbalance:**
    Instead of calculating a single KS statistic on the raw data, the class computes a **robust estimate**.
    * The resampling process draws a fixed number of observations from each subject (and a fixed number of subjects per group). This ensures that the test result is not dominated by subjects with higher number of measures.
    * It calculates the KS (or Anderson-Darling) statistic for each of these balanced resamples.
    * The final **Observed Statistic** is reported as the **median** of this distribution. This ensures the metric is stable and not driven by outliers.

2.  **Hierarchical Resampling (Permutation):**
    To test if this difference is significant, we need a null distribution (what the difference would look like by random chance).
    * **Standard Shuffling:** A standard test shuffles all *observations* randomly. This destroys the group structure.
    * **Hierarchical Shuffling:** The `HierarchicalPermutationTest` class shuffles the **Class Labels** (e.g., the Subject IDs) rather than the individual observations.
        * It keeps all observations from "Subject A" together.
        * It randomly assigns the *entire* "Subject A" group to either the "Control" or "Treatment" bin.

3.  **Computing the P-Value:**
    The class repeats this hierarchical shuffling thousands of times.
    * It counts how many times the shuffled groups produced a difference (KS statistic) as large or larger than the robust observed statistic calculated in Step 1.
This ensures that the p-value reflects the probability of seeing such a difference *given the number of independent groups (subjects)*, not the number of total observations.

