import functools
import pathlib

import appdirs
import daiquiri
import diskcache

logger = daiquiri.getLogger("cache")


def get_cache_dir():
    cache_dir = pathlib.Path(appdirs.user_cache_dir("tsbrowse", "tsbrowse"))
    cache_dir.mkdir(exist_ok=True, parents=True)
    return cache_dir


cache = diskcache.Cache(get_cache_dir())


def disk_cache(version):
    def decorator(func):
        @functools.wraps(func)
        def wrapper(self, *args, **kwargs):
            uuid = self.file_uuid
            if uuid is None:
                logger.info(f"No uuid, not caching {func.__name__}")
                return func(self, *args, **kwargs)
            key = f"{self.file_uuid}-{func.__name__}-{version}"
            if key in cache:
                logger.info(f"Fetching {key} from cache")
                return cache[key]
            logger.info(f"Calculating {key} and caching")
            result = func(self, *args, **kwargs)
            cache[key] = result
            return result

        return wrapper

    return decorator
