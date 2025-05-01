import os
import hashlib
# Make sure this is the correct import path for your app
# Adjust this import path to match your project structure
from crossmatcher_cache import (_hash_call, _save_hdf5, _load_hdf5)
from config import cache_dir as CACHE_DIR


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

    The results are stored in a .hdf5 file in `cache/<hash>.hdf5`.
    If that file already exists, the data are loaded from it
    instead of running the query.
    """
    # Ensure the cache directory exists
    if not os.path.exists(CACHE_DIR):
        os.makedirs(CACHE_DIR)

    # Create a hash from the query string plus relevant kwargs
    query_hash = _hash_query(query, kwargs)
    cache_file = os.path.join(CACHE_DIR, f"{query_hash}.hdf5")

    # If a cached file exists, try to load it
    if os.path.isfile(cache_file):
        data = _load_hdf5(cache_file)

        return data

    # Otherwise, run the query via sqlutil
    import sqlutilpy as sqlutil
    result = sqlutil.get(query, **kwargs)

    # Figure out whether the user requested a dictionary or a tuple
    as_dict = bool(kwargs.get("asDict", False))

    # Save data to .hdf5 for future calls
    _save_hdf5(cache_file, result, as_dict=as_dict)
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
    The result is cached into cache/<hash>.hdf5 if not already present.
    """

    # Ensure our cache folder exists
    if not os.path.exists(CACHE_DIR):
        os.makedirs(CACHE_DIR)

    # Build a unique hash based on all arguments: query, table, arrays,
    # colnames, and kwargs
    call_hash = _hash_call("local_join", query, table_name, arr_tuple,
                           colnames_tuple, **kwargs)
    cache_file = os.path.join(CACHE_DIR, f"{call_hash}.hdf5")

    # If we already have a cached result, load and return
    if os.path.isfile(cache_file):
        return _load_hdf5(cache_file)

    # Otherwise, run the actual sqlutil.local_join
    import sqlutilpy as sqlutil
    result = sqlutil.local_join(query, table_name, arr_tuple, colnames_tuple,
                                **kwargs)

    # local_join typically returns a tuple of arrays, so we store as a "tuple"
    # (as_dict=False)
    _save_hdf5(cache_file, result, as_dict=False)

    return result
