import panel as pn


class OverviewPage:
    key = "overview"
    title = "Overview"

    def __init__(self, tsm):
        self.tsm = tsm
        self.content = pn.Column(
            pn.pane.Markdown(f"## {self.tsm.full_path}"), pn.pane.HTML(self.tsm.ts)
        )
        self.sidebar = pn.Column(pn.pane.Markdown("# Overview"))
