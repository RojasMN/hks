import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
from typing import Literal, Optional
import warnings


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
      
        self.df = df.copy()
        self.value_col = value_col
        self.subject_col = subject_col
        self.cell_col = cell_col
        self.group_col = group_col
        self.metric = metric 
        self.replace = replace
      
        self.groups = self.df[group_col].unique()
        if len(self.groups) != 2:
            raise ValueError(f"The column '{group_col}' must contain exactly 2 groups.")
      
        # ------------ 1. Setting the number of observations per cell ------------
        
        actual_min_obs = self.df.groupby([subject_col, cell_col])[value_col].count().min()  
      
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
        
        # ------------ 2. Setting the number of cells per subject ------------
        
        actual_min_cells = self.df.groupby(subject_col)[cell_col].nunique().min()
      
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
        
        self.observed_distribution = []
        self.null_distribution = []
        self.observed_stat_median = None
        self.p_value = None
      
    # ------------ 3. Hierarchical resample ------------
    
    def _hierarchical_resample(self, df_input):
        
        # Balancing number of observations
        balanced_cells_df = df_input.groupby([self.subject_col, self.cell_col], group_keys=True).apply(
            lambda x: x.sample(n=self.n_obs, replace=self.replace)
        )
        
        balanced_cells_df = balanced_cells_df.reset_index()
        balanced_cells_df = balanced_cells_df.loc[:, ~balanced_cells_df.columns.duplicated()]
        
        # Balancing the number of cells 
        unique_cells = balanced_cells_df[[self.subject_col, self.cell_col]].drop_duplicates()
        selected_cells = unique_cells.groupby(self.subject_col).sample(n=self.n_cells, replace=self.replace)
        
        # MERGE
        final_df = balanced_cells_df.merge(selected_cells, on=[self.subject_col, self.cell_col], how='inner')
        
        return final_df
            
    # ------------ 4. Computing the KS/AD statistic ------------
    
    def _compute_statistic(self, current_df):
        
        g1_data = current_df[current_df[self.group_col] == self.groups[0]][self.value_col].values
        g2_data = current_df[current_df[self.group_col] == self.groups[1]][self.value_col].values
        
        if len(g1_data) == 0 or len(g2_data) == 0:
            return 0
        
        if self.metric == 'ks':
            stat, _ = stats.ks_2samp(g1_data, g2_data)
        
        elif self.metric == 'ad':
            
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                res = stats.anderson_ksamp([g1_data, g2_data])
                
                try:
                    stat = res.statistic
                except AttributeError:
                    stat = res[0]
                
        return stat
      
    # ------------ 5. Bootstrap ------------
    
    def run(self, n_resamples = 100, n_permutations = 1000, random_state = None):
        
        """
        The bootstrap resampling procedure is performed on the observed data to generate a
        distribution of KS or AD statistics and compute its median. Subsequently, a null 
        distribution is obtained by performing the same procedure on data with shuffled 
        subject labels.
        
        The p-value is calculated as the proportion of statistics in the null distribution 
        that are greater than or equal to the median statistic derived from the observed data.
        
        """
        
        if random_state:
            np.random.seed(random_state)
        
        # 5.1 Observed distribution
        obs_stats = []
        for _ in tqdm(range(n_resamples), desc = "Observed Distribution"):
            resampled_df = self._hierarchical_resample(self.df)
            stat = self._compute_statistic(resampled_df)
            obs_stats.append(stat)
        
        self.observed_distribution = np.array(obs_stats)
        self.observed_stat_median = np.median(self.observed_distribution)
        
        # 5.2 Null distribution
        null_stats = []
        subjects = self.df[[self.subject_col]].drop_duplicates()
        original_labels = self.df.set_index(self.subject_col)[self.group_col].to_dict()
        
        for _ in tqdm(range(n_permutations), desc = "Null Distribution"):
            shuffled_labels = np.random.permutation(list(original_labels.values()))
            subject_mapping = dict(zip(subjects[self.subject_col], shuffled_labels))
            
            permuted_df = self.df.copy()
            permuted_df[self.group_col] = permuted_df[self.subject_col].map(subject_mapping)
            
            resampled_perm_df = self._hierarchical_resample(permuted_df)
            stat = self._compute_statistic(resampled_perm_df)
            null_stats.append(stat)
            
        self.null_distribution = np.array(null_stats)
        
        # P-value
        self.p_value = np.mean(self.null_distribution >= self.observed_stat_median)
        
        return self.p_value
      
    # ------------ 6. Plotting ------------
    
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
        