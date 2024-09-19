def raster_component(page, png_filename):
    content = page.content
    content.save(png_filename, as_png=True)
