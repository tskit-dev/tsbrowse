import panel as pn


def page(tsm):
    return pn.Column(pn.pane.HTML(tsm.ts))
