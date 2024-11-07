(sec_nodes)=

# Nodes page

![Nodes Page](tsbrowse:example.tsbrowse:nodes)

A node defines a sampled or ancestral sequence represented in an ARG. It is identified by a numerical ID and can be present in many marginal trees. 

The interactive top plot visualises the total span on the sequence for each ancestral node over time. 

The histograms at the bottom show the distributions of node spans over different dimensions. The leftmost histogram summarises the span of nodes on the sequence; the middle plot summarises the span of nodes over time and the rightmost plot summarises the edge "area" defined as the product of sequence span and time span for each node.

## Plot controls (sidebar):
* The `Node flags`checkbox group 
    The tskit Nodes table includes a bitwise `flags` column used to store information about a node. For example. a value of 1 indicates that the node is a sample. By using these checkboxes it is possible to select which nodes to include in the node spans plot. based on their flag values.
* The `log y-axis` checkbox plots node time on log scale.