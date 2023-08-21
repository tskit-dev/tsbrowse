import panel as pn


def page(tsm):
    return pn.pane.HTML(tsm.ts)
