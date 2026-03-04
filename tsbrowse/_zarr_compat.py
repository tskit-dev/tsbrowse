import zarr

_ZARR_V3 = int(zarr.__version__.split(".")[0]) >= 3


def open_zip_store(path, mode):
    """Open a ZipStore compatible with zarr v2 and v3."""
    if _ZARR_V3:
        return zarr.storage.ZipStore(path, mode=mode)
    else:
        return zarr.storage.ZipStore(str(path), mode=mode)


def open_group_for_read(store):
    """Open a zarr group for reading (v2 or v3 format auto-detected)."""
    return zarr.open_group(store=store, mode="r")


def open_group_for_write(store):
    """Open a zarr group for writing in zarr v2 format."""
    if _ZARR_V3:
        return zarr.open_group(store=store, zarr_format=2)
    else:
        return zarr.open_group(store=store)
