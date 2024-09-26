import holoviews as hv
import numpy as np
import pandas as pd
import panel as pn

from .. import config


def make_afs_df(ts, span_normalise, afs_type, mode):
    polarised = afs_type == "folded"

    afs = ts.allele_frequency_spectrum(
        span_normalise=span_normalise, windows=None, polarised=polarised, mode=mode
    )
    num_samples = ts.num_samples
    afs_df = pd.DataFrame(
        {
            "allele_count": np.arange(1, num_samples + 1),
            "num_sites": afs[1:],
        }
    )
    return afs_df


def make_afs_panel(afs_df, log_bins, mode):
    len_df = len(afs_df)
    if log_bins:
        max_order = int(np.ceil(np.log10(len_df)))
        bin_edges = np.geomspace(1, 10**max_order, num=max_order + 1).astype(int)
        xrotation = 0
    else:
        num_bins = min(20, int(np.sqrt(len_df)))
        bin_edges = np.linspace(1, len_df, num_bins).astype(int)
        xrotation = 45

    labels = [f"{bin_edges[i]} - {bin_edges[i + 1]}" for i in range(len(bin_edges) - 1)]
    afs_df["bins"] = pd.cut(
        afs_df["allele_count"],
        bins=bin_edges,
        right=True,
        include_lowest=True,
        labels=labels,
    )

    summed_values = pd.DataFrame(
        afs_df.groupby("bins", observed=True)["num_sites"].sum().reset_index()
    )
    summed_values.columns = ["bins", "ac_sum"]

    afs_bar = hv.Bars(summed_values, kdims="bins", vdims="ac_sum").opts(
        height=config.PLOT_HEIGHT,
        responsive=True,
        xlabel="Allele count",
        ylabel="Number of sites",
        xrotation=xrotation,
        tools=["hover"],
        title=f"{mode} mode",
        bar_width=1,
        shared_axes=False,
    )
    return afs_bar


class FrequencySpectraPage:
    key = "sfs"
    title = "SFS"

    def __init__(self, tsm):
        log_bins_chk = pn.widgets.Checkbox(name="log-scale bins", value=True)
        span_normalise_chk = pn.widgets.Checkbox(name="span normalise", value=False)
        afs_type_radio = pn.widgets.RadioButtonGroup(
            name="SFS type", options=["folded", "unfolded"], value="folded"
        )
        plot_options = pn.Column(
            pn.pane.Markdown("# Frequency Spectra"),
            pn.Row(pn.Column(log_bins_chk, span_normalise_chk), afs_type_radio),
        )

        afs_site_df = pn.bind(
            make_afs_df, tsm.ts, span_normalise_chk, afs_type_radio, "site"
        )
        afs_site_bar = pn.bind(
            make_afs_panel, afs_df=afs_site_df, log_bins=log_bins_chk, mode="site"
        )
        afs_branch_df = pn.bind(
            make_afs_df, tsm.ts, span_normalise_chk, afs_type_radio, "branch"
        )
        afs_branch_bar = pn.bind(
            make_afs_panel, afs_df=afs_branch_df, log_bins=log_bins_chk, mode="branch"
        )

        site_count_table = pn.widgets.Tabulator(
            afs_site_df, show_index=False, disabled=True, layout="fit_data_stretch"
        )
        branch_count_table = pn.widgets.Tabulator(
            afs_branch_df, show_index=False, disabled=True, layout="fit_data_stretch"
        )
        layout = pn.Row(
            pn.Column(afs_site_bar, site_count_table),
            pn.Column(afs_branch_bar, branch_count_table),
        )

        self.content = pn.Column(layout, sizing_mode="stretch_width")
        self.sidebar = plot_options
