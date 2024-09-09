from tsbrowse import model


def raster_component(page, ts, png_filename, *, width=None, height=None, **kwargs):
    tsm = model.TSModel(ts)
    p = page(tsm, **kwargs)
    if width is not None:
        p.width = width
    if height is not None:
        p.height = height
    p.save(png_filename, as_png=True)
