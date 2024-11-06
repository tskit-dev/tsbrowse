import panel as pn

WANTED_COLUMNS = {
    "edges": [
        "id",
        "left",
        "right",
        "parent",
        "child",
        "branch_length",
        "span",
        "parent_time",
        "child_time",
        "metadata",
    ],
    "trees": [
        "id",
        "left",
        "right",
        "num_mutations",
        "num_sites",
        "total_branch_length",
        "max_internal_arity",
        "mean_internal_arity",
    ],
    "mutations": [
        "id",
        "site",
        "node",
        "derived_state",
        "parent",
        "time",
        "position",
        "num_parents",
        "num_descendants",
        "num_inheritors",
        "inherited_state",
        "metadata",
    ],
    "nodes": [
        "id",
        "flags",
        "time",
        "population",
        "individual",
        "num_mutations",
        "ancestors_span",
        "metadata",
    ],
    "sites": ["id", "position", "ancestral_state", "num_mutations", "metadata"],
    "individuals": ["id", "flags", "parents", "location", "metadata"],
    "populations": ["id", "metadata"],
    "migrations": ["id", "left", "right", "node", "source", "dest", "time", "metadata"],
    "provenances": ["id", "timestamp", "record"],
}


class TablesPage:
    key = "tables"
    title = "Tables"

    def __init__(self, tsm):
        self.tsm = tsm

        self.table_selector = pn.widgets.Select(
            name="Select Table",
            options=[
                "",
                "edges",
                "trees",
                "mutations",
                "nodes",
                "sites",
                "individuals",
                "populations",
                "migrations",
                "provenances",
            ],
            value="",
        )

        self.filter_input = pn.widgets.TextInput(
            name="Filter Expression",
            placeholder="Enter query expression (e.g., num_mutations > 5)",
        )

        self.filter_status = pn.pane.Markdown("")

        self.content = pn.Column(
            pn.bind(self.get_filtered_table, self.table_selector, self.filter_input),
        )
        self.sidebar = pn.Column(
            pn.pane.Markdown("# Tables"),
            self.table_selector,
            self.filter_input,
            self.filter_status,
        )

    def get_filtered_table(self, table_name, filter_expr):
        if table_name == "":
            return pn.pane.Markdown("# Select a table to display")

        df = getattr(self.tsm, f"{table_name}_df")
        # Limit and order the columns
        df = df[WANTED_COLUMNS[table_name]]

        if filter_expr:
            try:
                filtered_df = df.query(filter_expr)
                self.filter_status.object = f"Filter applied: {len(filtered_df)} rows"
                df = filtered_df
            except Exception as e:
                self.filter_status.object = f"Error in filter: {str(e)}"
        else:
            self.filter_status.object = f"No filter: {len(df)} rows"

        return pn.widgets.Tabulator(
            df,
            hidden_columns=["index"],
            pagination="remote",
            page_size=25,
            height=800,
        )
