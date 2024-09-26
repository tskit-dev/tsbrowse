import panel as pn


class TablesPage:
    key = "tables"
    title = "Tables"

    def __init__(self, tsm):
        self.tsm = tsm

        self.table_selector = pn.widgets.Select(
            name="Select Table",
            options=["", "edges", "trees", "mutations", "nodes", "sites"],
            value="",
        )

        self.filter_input = pn.widgets.TextInput(
            name="Filter Expression",
            placeholder="Enter query expression (e.g., num_mutations > 5)",
        )

        self.filter_status = pn.pane.Markdown("")

        self.content = pn.Column(
            pn.bind(self.get_filtered_table, self.table_selector, self.filter_input),
            sizing_mode="stretch_both",
        )
        self.sidebar = pn.Column(
            pn.pane.Markdown("# Tables"),
            self.table_selector,
            self.filter_input,
            self.filter_status,
            sizing_mode="stretch_height",
        )

    def get_filtered_table(self, table_name, filter_expr):
        if table_name == "":
            return pn.pane.Markdown("# Select a table to display")

        df = getattr(self.tsm, f"{table_name}_df")

        if filter_expr:
            try:
                filtered_df = df.query(filter_expr)
                self.filter_status.object = f"Filter applied: {len(filtered_df)} rows"
                df = filtered_df
            except Exception as e:
                self.filter_status.object = f"Error in filter: {str(e)}"
        else:
            self.filter_status.object = ""

        return pn.widgets.Tabulator(
            df,
            pagination="remote",
            page_size=25,
            height=800,
        )
