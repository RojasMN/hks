import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
from joblib import Parallel, delayed
import warnings
from typing import Literal, Optional, List, Dict

class HierarchicalPermutationTest:
    
    """
    Performs a hierarchical permutation test to compare the distributions of a nested 
    variable between two groups, using the Kolmogorov-Smirnov or Anderson-Darling statistic.
    This procedure handles nested data structures where subjects from two groups 
    (e.g., WT and Treatment) contain multiple measurement units (e.g., cells), 
    each yielding multiple distinct observations.
    
    Strategies for handling data structure:
    
    1. Hierarchical Resampling: Balances the weight of every cell and subject to
    mitigate the effects of unbalanced sample sizes (pseudoreplication).
    2. Subject-Level Permutation: Shuffles group labels at the subject level rather 
    than the observation level, strictly preserving intra-subject correlation.
    
    """
    
    
    def __init__(self,
                 df: pd.DataFrame,
                 value_col: str,
                 subject_col: str,
                 cell_col: str,
                 group_col: str,
                 metric: Literal['ks', 'ad'] = 'ks',
                 n_cells_per_subject: Optional[int] = None,
                 n_obs_per_cell: Optional[int] = None,
                 replace: bool = False):
        
        self.metric = metric
        self.replace = replace 
        
        self.groups = df[group_col].unique()
        if len(self.groups) != 2:
            raise ValueError(f"The column '{group_col}' must contain exactly 2 groups.")
        
        # ------------------------ 1. Setting the number of observations per cell ------------------------
        
        actual_min_obs = df.groupby([subject_col, cell_col])[value_col].count().min()
      
        if n_obs_per_cell is None:
            self.n_obs = actual_min_obs
            print(f"INFO: Auto-detected n_obs_per_cell = {self.n_obs} (minimum number across the entire dataset).")
        else:
            if not self.replace and n_obs_per_cell > actual_min_obs:
                raise ValueError(f"ERROR: Requested {n_obs_per_cell} observations, but the minimum available is {actual_min_obs}. "
                                 f"Reduce the number or set 'replace = True'.")
            elif self.replace and n_obs_per_cell > actual_min_obs:
                print(f"WARNING: Requested {n_obs_per_cell} observations (available: {actual_min_obs}). Oversampling (with replacement) will be used.")
            
            self.n_obs = n_obs_per_cell
            print(f"INFO: User defined n_obs_per_cell = {self.n_obs}.")
        
        # ------------------------ 2. Setting the number of cells per subject ------------------------
        
        actual_min_cells = df.groupby(subject_col)[cell_col].nunique().min()
      
        if n_cells_per_subject is None:
            self.n_cells = actual_min_cells
            print(f"INFO: Auto-detected n_cells_per_subject = {self.n_cells} (minimum number across the entire dataset).")
        else:
            if not self.replace and n_cells_per_subject > actual_min_cells:
                raise ValueError(f"ERROR: Requested {n_cells_per_subject} cells, but the minimum available is {actual_min_cells}. "
                                 f"Reduce the number or set 'replace = True'.")
            elif self.replace and n_cells_per_subject > actual_min_cells:
                print(f"WARNING: Requested {n_cells_per_subject} cells (available: {actual_min_cells}). Oversampling (with replacement) will be used.")

            self.n_cells = n_cells_per_subject
            print(f"INFO: User defined n_cells_per_subject = {self.n_cells}.")
            
        # ------------------------ 3. Data preprocessing ------------------------
        
        print("INFO: Optimizing data structure for fast resampling...")
        self.optimized_data = []
        self.subject_ids = df[subject_col].unique()
        self.group_map = {self.groups[0]: 0, self.groups[1]: 1}
        
        for sub in self.subject_ids:
            sub_df = df[df[subject_col] == sub]
            group_label = self.group_map[sub_df[group_col].iloc[0]]
            
            cells_data = []
            for cell in sub_df[cell_col].unique():
                vals = sub_df[sub_df[cell_col] == cell][value_col].to_numpy()
                cells_data.append(vals)
                
            self.optimized_data.append((group_label, cells_data))
            
        self.n_subjects = len(self.optimized_data)
        self.observed_distribution = []
        self.null_distribution = []
        self.p_value = None 
    
    # --------------------------- 4. Resampling ---------------------------
    
    def _fast_resample_and_compute(self, data_structure, current_labels, seed):
        
        # Initialize unique RNG for this worker
        rng = np.random.default_rng(seed)
        
        g0_samples = []
        g1_samples = []
        
        for i in range(self.n_subjects):
            _, cells_list = data_structure[i]
            assigned_group = current_labels[i]
            n_available_cells = len(cells_list)
            
            # -------------- 4.1 Select cells from each subject --------------
       
            if self.replace:
                chosen_cell_indices = rng.choice(n_available_cells, self.n_cells, replace = True)
            else:
                chosen_cell_indices = rng.choice(n_available_cells, self.n_cells, replace = False)
            
            # -------------- 4.2 Sample observations from each selected cell --------------
            
            subject_samples = []
            for c_idx in chosen_cell_indices:
                cell_vals = cells_list[c_idx]
                if self.replace:
                    obs = rng.choice(cell_vals, self.n_obs, replace = True)
                else:
                    obs = rng.choice(cell_vals, self.n_obs, replace = False)
                subject_samples.append(obs)

            subject_samples = np.concatenate(subject_samples)
            
            if assigned_group == 0:
                g0_samples.append(subject_samples)
            else:
                g1_samples.append(subject_samples)
        
        if not g0_samples or not g1_samples:
            return 0
        
        g0_arr = np.concatenate(g0_samples)
        g1_arr = np.concatenate(g1_samples)
        
        # -------------- 4.3 Compute statistic --------------
        
        if self.metric == 'ks':
            return stats.ks_2samp(g0_arr, g1_arr).statistic
        elif self.metric == 'ad':
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                return stats.anderson_ksamp([g0_arr, g1_arr]).statistic
            
    # --------------------------- 5. Bootstrap ---------------------------
    
    def run(self, n_resamples = 1000, n_permutations = 1000, n_jobs = -1, random_state = None):
        
        """
        The bootstrap resampling procedure is performed on the observed data to generate a
        distribution of KS or AD statistics and compute its median. Subsequently, a null 
        distribution is obtained by performing the same procedure on data with shuffled 
        subject labels.
        
        The p-value is calculated as the proportion of statistics in the null distribution 
        that are greater than or equal to the median statistic derived from the observed data.
        
        """
        
        if random_state is None:
            master_rng = np.random.default_rng()
        else:
            master_rng = np.random.default_rng(random_state)
            
        original_labels = np.array([item[0] for item in self.optimized_data])
        
        # -------------- 5.1 Observed distribution --------------
        
        print(f"Calculating Observed Distribution ({n_resamples} resamples)...")
        obs_seeds = master_rng.integers(0, 10**9, size = n_resamples)
        
        self.observed_distribution = Parallel(n_jobs = n_jobs)(
            delayed(self._fast_resample_and_compute)(self.optimized_data, original_labels, seed)
            for seed in tqdm(obs_seeds, total = n_resamples, desc = "Observed")
        )
        
        self.observed_distribution = np.array(self.observed_distribution)
        self.observed_stat_median = np.median(self.observed_distribution)
        
        # -------------- 5.2 Null distribution --------------
        
        print(f"Calculating Null Distribution ({n_permutations} permutations)...")
        
        null_seeds = master_rng.integers(0, 10**9, size = n_permutations)
        permuted_label_sets = [master_rng.permutation(original_labels) for _ in range(n_permutations)]
            
        self.null_distribution = Parallel(n_jobs = n_jobs)(
            delayed(self._fast_resample_and_compute)(self.optimized_data, p_labels, seed)
            for p_labels, seed in tqdm(zip(permuted_label_sets, null_seeds), total = n_permutations, desc = "Null")
        )
        
        self.null_distribution = np.array(self.null_distribution)
        
        # -------------- 5.3 P-Value --------------
        
        self.p_value = np.mean(self.null_distribution >= self.observed_stat_median)
        return self.p_value
    
    # --------------------------- 6. Plotting ---------------------------
    
    def plot_results(self, title = None):
        
        if self.p_value is None:
            raise RuntimeError("Run .run() first")
        
        plt.figure(figsize = (10, 5))
        sns.kdeplot(self.null_distribution, fill = True, color = "grey", label = "Null Distribution (Permuted)")
        sns.kdeplot(self.observed_distribution, fill = True, color = "blue", label = "Observed Distribution (Resampled)")
        plt.axvline(self.observed_stat_median, color = 'red', linestyle = '--', label = f'Observed Median: {self.observed_stat_median:.3f}')
        
        metric_name = "KS Statistic" if self.metric == 'ks' else "Anderson-Darling Statistic"
        plt.xlabel(metric_name)
        plt.title(title or f"Nested Permutation Test ({self.metric.upper()})\np-value: {self.p_value:.5f}")
        plt.legend()
        plt.show()
        