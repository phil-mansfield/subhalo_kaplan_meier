# Subhalo Survival Analysis with Kaplan-Meier
A lightweight Python library for performing Kaplan-Meier analysis on the survival statistics of subhaloes. Mansfield et al. 2023 recommends several novel tests to assess the reliability of subhalo finders, many of which are implemented in this library.

Code for converting your merger tree into a grid format and computing some of the Mansfield et al. 2023 branch annotations for it.
```python
import symlib

# Read in merger tree quantities. You'll need to write these yourself. the arrays
# returned by read_my_merger_tree can be in any order (depth-first, depth-first, random,
# etc.) but needs to include the full evolution of all haloes.
mpeak_cutoff = my_mpeak_cutoff()
id, pid, leaf_id, snap, mvir, vmax = read_my_merger_tree()

# Convert into a numpy strucutred array.
tree = symlib.structured_array(
    [id, pid, leaf_id, snap, mvir, vmax],
    ["id", "pid", "leaf_id, "snap", "mvir", "vmax"],
    ["i8", "i8", "i8", "i8", "f4", "f4"]
)
# Convert to a (halo, snap) grid.
haloes, ok, _ = symlib.tree_to_grid(tree, haloes["leaf_id"], haloes["snap"])

# Make intial cut to reduce computation time.
mpeak_initial = symlib.mpeak_inital(halos["mvir"])
cut = mpeak_initial > mpeak_cutoff
haloes, ok = halos[cut,:], ok[cut,:]

# Compute annotated branch values.
infall_snap = symlib.infall_snapshot(haloes["id"], haloes["pid"], mpeak_initial)
mpeak_pre = symlib.mpeak_pre(haloes["mvir"], infall_snap)

# Make true cut on catalogue
cut = mpeak_pre > mpeak_cutoff
haloes, infall_snap, mpeak_pre = haloes[cut], infall_snap[cut], mpeak_pre[cut]
```

Compute subhalo-testing quantities
```python
# Select for subhaloes
is_sub = infall_snap > 0
haloes, ok = haloes[is_sub], ok[is_sub]
infall_snap, mpeak_pre = infall_snap[is_sub], mpeak_pre[is_sub]

# Survival curves
mu_eval = 10**np.linspace(-4, 0, 40)
S, S_err = symlib.survival_curve(haloes["mvir"], mpeak_pre)

# Mass-dependent mu_90 cutoffs

```
