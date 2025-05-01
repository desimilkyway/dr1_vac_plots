import os
import hashlib
import numpy as np
# Make sure this is the correct import path for your app
# Adjust this import path to match your project structure
from crossmatcher_cache import (_hash_call, _save_npz, _load_npz)
from config import cache_dir

CACHE_DIR = cache_dir


def _hash_query(query, kwargs):
    """
    Produce a unique hash based on both the query string
    and the relevant kwargs (so that changing asDict, etc.
    leads to a different cache file).
    """
    hasher = hashlib.sha256()
    hasher.update(query.encode("utf-8"))
    # Append certain kwargs to the hash (at minimum asDict).
    # You could include others if needed.
    as_dict_val = kwargs.get("asDict", False)
    hasher.update(f"asDict={as_dict_val}".encode("utf-8"))

    # If you want to incorporate all kwargs (besides asDict),
    # do something like this:
    # for k, v in sorted(kwargs.items()):
    #     if k != 'asDict':
    #         hasher.update(str(k).encode('utf-8'))
    #         hasher.update(str(v).encode('utf-8'))

    return hasher.hexdigest()


def get(query, **kwargs):
    """
    A caching wrapper around sqlutil.get().

    - If kwargs['asDict'] is True, we expect a dictionary of
      {colname -> np.array}.
    - Otherwise, we expect a tuple of np arrays.

    The results are stored in a .npz file in `cache/<hash>.npz`.
    If that file already exists, the data are loaded from it
    instead of running the query.
    """
    # Ensure the cache directory exists
    if not os.path.exists(CACHE_DIR):
        os.makedirs(CACHE_DIR)

    # Create a hash from the query string plus relevant kwargs
    query_hash = _hash_query(query, kwargs)
    cache_file = os.path.join(CACHE_DIR, f"{query_hash}.npz")

    # If a cached file exists, try to load it
    if os.path.isfile(cache_file):
        loaded = np.load(cache_file, allow_pickle=False)

        # Check if the data were written as a dictionary or a tuple
        is_dict = bool(loaded["IS_DICT"][0])  # 1 => dict, 0 => tuple
        if is_dict:
            # Reconstruct dictionary
            data = {}
            for k in loaded.files:
                if k == "IS_DICT" or k == "NUM_ARRAYS":
                    continue
                data[k] = loaded[k]
        else:
            # Reconstruct tuple
            num_arrays = loaded["NUM_ARRAYS"][0]
            data = tuple(loaded[f"arr_{i}"] for i in range(num_arrays))

        return data

    # Otherwise, run the query via sqlutil
    import sqlutilpy as sqlutil
    result = sqlutil.get(query, **kwargs)

    # Figure out whether the user requested a dictionary or a tuple
    as_dict = bool(kwargs.get("asDict", False))

    # Save data to .npz for future calls
    if as_dict:
        # Expecting result to be {colname -> np.array}
        arrays = {}
        arrays["IS_DICT"] = np.array([1], dtype=np.int8)

        # Each key in the dictionary becomes an array in the .npz
        for k, v in result.items():
            arrays[k] = v

        np.savez(cache_file, **arrays)
    else:
        # Expecting a tuple of arrays
        arrays = {}
        arrays["IS_DICT"] = np.array([0], dtype=np.int8)
        arrays["NUM_ARRAYS"] = np.array([len(result)], dtype=np.int32)

        # Store each array with a known name
        for i, arr in enumerate(result):
            arrays[f"arr_{i}"] = arr

        np.savez(cache_file, **arrays)

    return result


def local_join(query, table_name, arr_tuple, colnames_tuple, **kwargs):
    """
    A caching wrapper around sqlutil.local_join(...).

    Usage example:
        CA_dr3_id, = sqlutil_cache.local_join(
            '''
            WITH x AS (
              SELECT rid, dr3_source_id,
                     ROW_NUMBER() OVER (PARTITION BY rid ORDER BY
                    angular_distance ASC)
              FROM gaia_edr3.dr2_neighbourhood AS g, m
              WHERE m.source_id = g.dr2_source_id
            )
            SELECT dr3_source_id
            FROM x
            WHERE row_number = 1
            ORDER BY rid
            ''',
            'm',
            (rid_array, xid_array),
            ('rid','source_id')
        )

    This will return a tuple of numpy arrays, just as the real
    # sqlutil.local_join does.
    The result is cached into cache/<hash>.npz if not already present.
    """

    # Ensure our cache folder exists
    if not os.path.exists(CACHE_DIR):
        os.makedirs(CACHE_DIR)

    # Build a unique hash based on all arguments: query, table, arrays,
    # colnames, and kwargs
    call_hash = _hash_call("local_join", query, table_name, arr_tuple,
                           colnames_tuple, **kwargs)
    cache_file = os.path.join(CACHE_DIR, f"{call_hash}.npz")

    # If we already have a cached result, load and return
    if os.path.isfile(cache_file):
        return _load_npz(cache_file)

    # Otherwise, run the actual sqlutil.local_join
    import sqlutilpy as sqlutil
    result = sqlutil.local_join(query, table_name, arr_tuple, colnames_tuple,
                                **kwargs)

    # local_join typically returns a tuple of arrays, so we store as a "tuple"
    # (as_dict=False)
    _save_npz(cache_file, result, as_dict=False)

    return result
