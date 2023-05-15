import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tsinfer
import tskit
import tszip

"""
QC functions for tsinfer trees
"""
class TreeInfo:
    """
    Class for storing tree information
    """
    def __init__(self, ts, chr):
        self.ts = ts
        self.trees = ts.trees
        self.num_samples = ts.num_samples
        self.sequence_length = ts.sequence_length
        self.chr = chr

    def calc_polytomy_fractions(self):
        """
        Calculates the fraction of polytomies for each tree in the tree sequence
        """
        assert self.num_samples > 2
        polytomy_fractions = []
        for tree in self.trees():
            if (tree.num_edges == 0):
                polytomy_fractions.append(None)
            else:
                polytomy_fractions.append(float((tree.num_edges - self.num_samples)/(self.num_samples - 2)))
        return polytomy_fractions

    def map_stats_to_genome(self, to_map):
        """
        Converts a list of tree-based stats to genomic coordinates
        """
        mapped = np.zeros(int(self.sequence_length))
        for i, tree in enumerate(self.trees()):
            left, right = map(int, tree.interval)     
            mapped[left: right] = to_map[i] 
        return(mapped)
    
    def make_sliding_windows(self, iterable, size, overlap=0):
        start = 0
        assert overlap < size, "overlap must be smaller then window size"
        end = size
        step = size - overlap
        
        length = len(iterable)
        while end < length:
            yield iterable[start:end]
            start += step
            end += step
        yield iterable[start:]

    def plot_polytomy_fractions(self, zoom_start=None, zoom_end=None, window_size=100_000, overlap=0):
        """
        Plots the fraction of polytomies in windows actoss the genomic sequence
        """
        if zoom_start is None:
            zoom_start = 0
        if zoom_end is None:
            zoom_end = self.sequence_length
        fig, ax = plt.subplots(figsize=(20, 5))
        polytomy_fractions = self.calc_polytomy_fractions()
        poly_fracs_by_pos = self.map_stats_to_genome(polytomy_fractions)
        poly_fracs_means = []
        poly_fracs_sd = []
        genomic_positions = []
        for poly_win in self.make_sliding_windows(poly_fracs_by_pos, window_size, overlap):
            poly_fracs_means.append(np.mean(poly_win))
            poly_fracs_sd.append(np.std(poly_win))
        for gen_win in self.make_sliding_windows(np.arange(1, self.sequence_length), window_size, overlap):
            genomic_positions.append(gen_win[0]/1_000_000)
        ax.plot(genomic_positions, poly_fracs_means, color="grey", label="mean", linewidth=0.5)
        ax.fill_between(genomic_positions, np.array(poly_fracs_means) - np.array(poly_fracs_sd), np.array(poly_fracs_means) + np.array(poly_fracs_sd), color="blue", alpha=0.3, label="mean +/- std")
        missing_vals = np.take(genomic_positions,np.where(np.isnan(poly_fracs_means)))
        ax.plot(missing_vals, np.zeros(len(missing_vals)), color="red", marker = "o", label = "missing data")
        ax.set_xlabel(f"Position on chr {self.chr}(Mb)", fontsize=10)
        ax.set_ylabel("Window mean", fontsize=10)
        ax.set_title("Polytomy score", fontsize=10)
        ax.set_ylim(0, 1)
        ax.set_xlim(zoom_start/1_000_000, zoom_end/1_000_000)
        #ax.legend(fontsize=12, numpoints = 1)
        handles, labels = ax.get_legend_handles_labels()
        unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
        ax.legend(*zip(*unique))
        plt.show()
#main
