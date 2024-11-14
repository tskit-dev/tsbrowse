(sec_tables)=

# Tables page

![Tables Page](tsbrowse:example.tsbrowse:tables)

In tskit, tree sequences are represented as a collection of tables.
This page allows the querying and inspection of the raw data in these tables.
For more detail on what the columns mean see the [tskit data model](https://tskit.dev/tskit/docs/stable/data-model.html). Note that additional columns, used by tsbrowse, are added to the data model and displayed here. These include convenience columns such as the number of 
inheritors for a given mutation.


## Controls (sidebar):
There are two sidebar controls, the first allows the user to select which table to display, and the second allows the user to filter the rows of the table based on a Python-like expression. This expression is passed to the [pandas query method](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.query.html). Column names can be used as variables in the expression, and the expression should evaluate to a boolean value.
For example, to filter the mutations table to only show mutations with a derived state of 'C' that have more than 20,000 inheritors, you could use the expression 
`derived_state == 'C' and num_inheritors > 20000`.
