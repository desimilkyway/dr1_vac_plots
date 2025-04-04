# crossmatcher_cache.py

import os
import hashlib
import numpy as np

CACHE_DIR = "cache"


def _encode_for_hash(item):
    """
    Convert 'item' into a bytes-like object for hashing.
    - If it's a numpy array, include shape, dtype, and raw bytes of data.
    - Otherwise, just use its string representation.
    """
    if isinstance(item, np.ndarray):
        # Combine shape, dtype, and raw bytes:
        shape_bytes = str(item.shape).encode('utf-8')
        dtype_bytes = str(item.dtype).encode('utf-8')
        data_bytes = item.tobytes()
        return shape_bytes + dtype_bytes + data_bytes
    elif isinstance(item, (list, tuple)):
        # Possibly you have lists or tuples of numeric data.
        # You could either turn them into np.array or hash them recursively.
        # Here, let's just convert them to np array for consistency:
        arr = np.array(item)
        return _encode_for_hash(arr)
    else:
        # Fallback: string representation
        return str(item).encode('utf-8')


def _hash_call(method_name, *args, **kwargs):
    """
    Produce a SHA256 hash capturing:
      - The method name (e.g. 'doit_by_key' or 'doit')
      - Each positional argument
      - Keyword arguments (including asDict, db, plus any arrays)
    """
    hasher = hashlib.sha256()

    # Include the method (so 'doit_by_key' vs 'doit' differ)
    hasher.update(method_name.encode("utf-8"))

    # Include all positional args
    for arg in args:
        hasher.update(_encode_for_hash(arg))

    # Include all kwargs (in sorted order for consistency)
    for k in sorted(kwargs.keys()):
        hasher.update(k.encode("utf-8"))
        hasher.update(_encode_for_hash(kwargs[k]))

    return hasher.hexdigest()


def _save_npz(filename, result, as_dict):
    """
    Save 'result' to an .npz file.
    If as_dict=True, result is a dict {colName -> np.array}.
    Otherwise result is a tuple of arrays.
    """
    if as_dict:
        # Expecting a dictionary of {colname -> np.array}
        arrays = {}
        arrays["IS_DICT"] = np.array([1], dtype=np.int8)

        # Each key in the dictionary becomes an array in the .npz
        for k, v in result.items():
            arrays[k] = v

        np.savez(filename, **arrays)
    else:
        # Expecting a tuple of arrays
        arrays = {}
        arrays["IS_DICT"] = np.array([0], dtype=np.int8)
        arrays["NUM_ARRAYS"] = np.array([len(result)], dtype=np.int32)

        # Store each array with a known name
        for i, arr in enumerate(result):
            arrays[f"arr_{i}"] = arr

        np.savez(filename, **arrays)


def _load_npz(filename):
    """
    Load from an .npz file created by _save_npz.
    Return either a dict or tuple depending on the 'IS_DICT' flag.
    """
    loaded = np.load(filename, allow_pickle=False)
    is_dict = bool(loaded["IS_DICT"][0])  # 1 => dict, 0 => tuple

    if is_dict:
        # Reconstruct dictionary
        data = {}
        for k in loaded.files:
            if k in ("IS_DICT", "NUM_ARRAYS"):
                continue
            data[k] = loaded[k]
    else:
        # Reconstruct tuple
        num_arrays = loaded["NUM_ARRAYS"][0]
        data = tuple(loaded[f"arr_{i}"] for i in range(num_arrays))

    return data


def doit_by_key(table, key_values, columns, **kwargs):
    """
    A caching wrapper around crossmatcher.doit_by_key(...).

    Example usage:

        D_AP = crossmatcher_cache.doit_by_key(
            'apogee_dr17.allstar',
            G_T['SOURCE_ID'],  # np array or list of keys
            'alpha_m,fe_h,c_fe,...',
            key_col='gaiaedr3_source_id',
            db='wsdb',
            asDict=True
        )
    """
    # Ensure the cache directory exists
    if not os.path.exists(CACHE_DIR):
        os.makedirs(CACHE_DIR)

    # Figure out whether user wants a dict back
    as_dict = bool(kwargs.get("asDict", False))

    # Build a unique hash for this call
    call_hash = _hash_call("doit_by_key", table, key_values, columns, **kwargs)
    cache_file = os.path.join(CACHE_DIR, f"{call_hash}.npz")

    # If it exists, load from cache
    if os.path.isfile(cache_file):
        return _load_npz(cache_file)

    # Otherwise, actually call crossmatcher

    import crossmatcher
    result = crossmatcher.doit_by_key(table, key_values, columns, **kwargs)

    # Save to cache
    _save_npz(cache_file, result, as_dict)

    return result


def doit(table, ra_values, dec_values, columns, **kwargs):
    """
    A caching wrapper around crossmatcher.doit(...).

    Example usage:

        D_AP = crossmatcher_cache.doit(
            'apogee_dr17.allstar',
            TT['TARGET_RA'],    # np array or list
            TT['TARGET_DEC'],   # np array or list
            'vhelio_avg, aspcapflag, ...',
            db='wsdb',
            asDict=True
        )
    """
    # Ensure the cache directory exists
    if not os.path.exists(CACHE_DIR):
        os.makedirs(CACHE_DIR)

    # Figure out whether user wants a dict back
    as_dict = bool(kwargs.get("asDict", False))

    # Build a unique hash for this call
    call_hash = _hash_call("doit", table, ra_values, dec_values, columns,
                           **kwargs)
    cache_file = os.path.join(CACHE_DIR, f"{call_hash}.npz")

    # If it exists, load from cache
    if os.path.isfile(cache_file):
        return _load_npz(cache_file)

    # Otherwise, call crossmatcher
    import crossmatcher
    result = crossmatcher.doit(table, ra_values, dec_values, columns, **kwargs)

    # Save to cache
    _save_npz(cache_file, result, as_dict)

    return result
