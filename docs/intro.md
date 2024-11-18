(sec_intro)=

# Introduction to tsbrowse
Tsbrowse is an open-source Python web-app intended as a genome browser for ARGs. It provides interactive visual summaries of the fundamental components of an ARG, i.e. the mutations, nodes, edges and trees it encodes. Tsbrowse operates on ARGs encoded as succinct tree sequences (see [tskit](https://tskit.dev/tskit/docs/stable/introduction.html)). It can be used to inspect biobank-scale ARGs without requiring data uploads to a server.

## Data Model
The foundation of tsbrowse is the [tskit data model](https://tskit.dev/tskit/docs/stable/data-model.html), which is a generalised ARG representation. Tsbrowse can therefore be applied to ARGs inferred using widely-used methods including [tsinfer](https://tskit.dev/tsinfer/docs/stable/), [Relate](https://myersgroup.github.io/relate/), [ARG-Needle](https://palamaralab.github.io/software/argneedle/) and [SINGER](https://github.com/popgenmethods/SINGER). Tsbrowse augments the underlying tree sequence tables with additional columns to create a compressed `.tsbrowse` file, containing all necessary semantic metadata required for automatic visualisation with the [Holoviz](https://holoviz.org/) ecosystem. 

## How to
Consider an example tree sequence simulated with msprime:
```
def make_sweep_ts(n, Ne, L, rho=1e-8, mu=1e-8):
    sweep_model = msprime.SweepGenicSelection(
        position=L/2, start_frequency=0.0001, end_frequency=0.9999, s=0.25, dt=1e-6)
        models = [sweep_model, msprime.StandardCoalescent()]
    ts = msprime.sim_ancestry(
        n, model=models, population_size=Ne, sequence_length=L, recombination_rate=rho, random_seed=1234)
    return msprime.sim_mutations(ts, rate=mu, random_seed=4321)

sim_ts = make_sweep_ts(300, Ne=10_000, L=5_000_000, mu=1e-8)
sim_ts.dump("example.trees")
```
To run tsbrowse on this tree sequence, we first run the preprocessing step.
```
python -m tsbrowse preprocess example.trees
```
This creates a `.tsbrowse` file in the same directory as the input tree sequence.
We can then run the serve step to view the app in a web browser:
```
python -m tsbrowse serve example.tsbrowse --port 8080
```