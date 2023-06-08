---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.5
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

```python
import sys
import importlib

import utils
import tszip
import tskit

importlib.reload(utils)

# ts = tskit.load("data/hgdp_tgp_sgdp_high_cov_ancients_chr20_q.dated.trees")
ts = tszip.decompress("data/hgdp_tgp_sgdp_chr2_q.dated.trees.tsz")
tinfo = utils.TreeInfo(ts, 6)
tinfo
```

```python
tinfo.plot_polytomy_fractions(window_size=500_000, overlap=0)
```

```python
tinfo.plot_mutations_per_site()
```

```python
tinfo.plot_mutations_per_site_along_seq(region_start=None, region_end=None)
```

```python
tinfo.plot_mutations_per_node(show_counts=True, max_num_muts=10)
```

```python
tinfo.plot_tree_spans(log_transform=True, region_end=200_000_000, region_start=190_000_000, show_counts=True)
```

```python
tinfo.plot_mean_node_arity(show_counts=True)
```

```python
tinfo.plot_mutations_per_tree(show_counts=True)
```

```python
tinfo.plot_mutations_per_tree_along_seq(hist_bins=500,region_start=150_000_000, region_end=190_000_000)
```

```python
tinfo.plot_sites_per_tree(show_counts=True)
```

```python
tinfo.plot_sites_per_tree_along_seq(hist_bins=500,region_start=150_000_000, region_end=190_000_000)
```
