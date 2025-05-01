# crossmatcher_cache.py

import os
import hashlib
import numpy as np
import h5py
from config import cache_dir

CACHE_DIR = cache_dir


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


def _save_hdf5(filename, result, as_dict):
    """
    Save 'result' to an .hdf5 file.

    If as_dict=True, result is a dict {colName -> np.array}.
    Otherwise result is a tuple of arrays.
    """
    with h5py.File(filename, 'w') as f:
        if as_dict:
            # Expecting a dictionary of {colname -> np.array}
            f.attrs['IS_DICT'] = np.int8(1)
            for k, v in result.items():
                f.create_dataset(k, data=v)
        else:
            # Expecting a tuple of arrays
            f.attrs['IS_DICT'] = np.int8(0)
            f.attrs['NUM_ARRAYS'] = np.int32(len(result))
            # Store each array with a known name
            for i, arr in enumerate(result):
                f.create_dataset(f"arr_{i}", data=arr)


def _load_hdf5(filename):
    """
    Load from an .hdf5 file created by _save_hdf5.

    Return either a dict or tuple depending on the 'IS_DICT' attribute.
    """
    with h5py.File(filename, 'r') as f:
        is_dict = bool(f.attrs['IS_DICT'])  # 1 => dict, 0 => tuple

        if is_dict:
            # Reconstruct dictionary
            data = {}
            for k in f.keys():
                data[k] = np.array(f[k])
        else:
            # Reconstruct tuple
            num_arrays = f.attrs['NUM_ARRAYS']
            data = tuple(np.array(f[f"arr_{i}"]) for i in range(num_arrays))

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
    cache_file = os.path.join(CACHE_DIR, f"{call_hash}.hdf5")

    # If it exists, load from cache
    if os.path.isfile(cache_file):
        return _load_hdf5(cache_file)

    # Otherwise, actually call crossmatcher

    import crossmatcher
    result = crossmatcher.doit_by_key(table, key_values, columns, **kwargs)

    # Save to cache
    _save_hdf5(cache_file, result, as_dict)

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
    cache_file = os.path.join(CACHE_DIR, f"{call_hash}.hdf5")

    # If it exists, load from cache
    if os.path.isfile(cache_file):
        return _load_hdf5(cache_file)

    # Otherwise, call crossmatcher
    import crossmatcher
    result = crossmatcher.doit(table, ra_values, dec_values, columns, **kwargs)

    # Save to cache
    _save_hdf5(cache_file, result, as_dict)

    return result
