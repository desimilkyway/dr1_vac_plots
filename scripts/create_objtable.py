#!/usr/bin/env python3

import re
import math
import sys
from astropy.table import Table


def capitalize_words(s):
    """
    Capitalize the first letter of each word in the given string.
    Example: "canes venatici 1" -> "Canes Venatici 1"
    """
    return " ".join(w.capitalize() for w in s.split())


def transform_gc_name(raw_name):
    """
    If name is of the form NGC_XXXX_M_YY, return 'M YY'.
    Otherwise, preserve uppercase 'NGC' and underscores->spaces,
    but do not force-lower or force-capitalize the entire string.
    """
    pattern = r'^NGC_\d+_M_(\d+)$'
    match = re.match(pattern, raw_name)
    if match:
        # e.g., "NGC_5024_M_53" -> "M 53"
        return f"M {match.group(1)}"
    else:
        # Just underscores->spaces; keep the case as is
        return raw_name.replace('_', ' ')


def transform_dsph_name(raw_name):
    """
    Replace underscores with spaces and uppercase the first letter of each word.
    E.g., "sextans_1" -> "Sextans 1"
    """
    return capitalize_words(raw_name.replace('_', ' '))


def transform_stream_name(raw_name):
    """
    Replace underscores with spaces, special-case "Pal-5" -> "Pal 5",
    then uppercase the first letter of each word.
    """
    name = raw_name.replace('_', ' ')
    # Special replacement
    name = name.replace("Pal-5", "Pal 5")
    pattern = r'^M(\d+)$'
    match = re.match(pattern, name)
    if match:
        # e.g., "NGC_5024_M_53" -> "M 53"
        name = f"M {match.group(1)}"
    return name


def transform_default_name(raw_name):
    """Default transformation: underscores->spaces, then capitalize each word."""
    return (raw_name.replace('_', ' '))


def format_percentiles(p16_str, p50_str, p84_str):
    """
    Formats percentiles into LaTeX string: $Value_{-err_low}^{+err_up}$
    Returns an empty string "" if any input is 'nan', cannot be converted,
    or results in a nan float.
    """
    nan_strings = {'nan'}
    # Convert to string in case row has numeric/numpy types
    p16_str, p50_str, p84_str = str(p16_str), str(p50_str), str(p84_str)

    if (p16_str.lower() in nan_strings or p50_str.lower() in nan_strings
            or p84_str.lower() in nan_strings):
        return ""

    try:
        p16_f = float(p16_str)
        p50_f = float(p50_str)
        p84_f = float(p84_str)

        if math.isnan(p16_f) or math.isnan(p50_f) or math.isnan(p84_f):
            return ""

        val = p50_f
        err_low = val - p16_f
        err_up = p84_f - val
        return f"${val:.2f}_{{-{err_low:.2f}}}^{{+{err_up:.2f}}}$"
    except ValueError:
        return ""


def create_minimal_table(entries):
    """
    Minimal table: Type, Name, Count, in 3 column groups
    -> so each group is (Type & Name & Count).
    """
    num_output_columns = 3
    total = len(entries)
    if total == 0:
        print("No valid entries found for minimal table.")
        return []

    # Divide entries into column groups
    groups = []
    start_index = 0
    base_size = total // num_output_columns
    remainder = total % num_output_columns

    for i in range(num_output_columns):
        size = base_size + (1 if i < remainder else 0)
        groups.append(entries[start_index:start_index + size])
        start_index += size

    # Prepare LaTeX
    latex_lines = []
    latex_lines.append("\\begin{table*}[htbp]")
    latex_lines.append("\\centering")
    latex_lines.append("\\caption{Minimal table: Type, Name, Count}")
    latex_lines.append("\\label{tab:mytable_minimal}")

    # Each group has 3 columns
    tabular_cols = " ".join(["lll"] * num_output_columns)
    latex_lines.append(f"\\begin{{tabular}}{{{tabular_cols}}}")
    latex_lines.append("\\hline")

    # Build the header
    header_parts = []
    col_header = "\\textbf{Type} & \\textbf{Name} & \\textbf{Count}"
    for _ in range(num_output_columns):
        header_parts.append(col_header)
    latex_lines.append(" & ".join(header_parts) + " \\\\")
    latex_lines.append("\\hline")

    # Number of rows
    max_rows = max(len(g) for g in groups)
    for i in range(max_rows):
        row_parts = []
        for j in range(num_output_columns):
            if i < len(groups[j]):
                e = groups[j][i]
                row_parts.append(f"{e['type']} & {e['name']} & {e['count']}")
            else:
                row_parts.append("& & ")
        latex_lines.append(" & ".join(row_parts) + " \\\\")
    latex_lines.append("\\hline")
    latex_lines.append("\\end{tabular}")
    latex_lines.append("\\end{table*}")

    return latex_lines


def create_full_table(entries):
    """
    Full table: do NOT output Type or Count. Instead, each group has 5 columns:
      1) Name
      2) [Fe/H]_RVS
      3) Sigma_RVS
      4) [Fe/H]_SP
      5) Sigma_SP

    We keep 2 groups of columns. So 2 sets of these 5 columns -> total 10 columns.
    """
    num_output_columns = 2

    # Convert the percentile columns to latex strings
    for e in entries:
        e['feh_rvs_str'] = format_percentiles(e['feh_rvs_1'], e['feh_rvs_2'],
                                              e['feh_rvs_3'])
        e['sfeh_rvs_str'] = format_percentiles(e['sfeh_rvs_1'],
                                               e['sfeh_rvs_2'],
                                               e['sfeh_rvs_3'])
        e['feh_sp_str'] = format_percentiles(e['feh_sp_1'], e['feh_sp_2'],
                                             e['feh_sp_3'])
        e['sfeh_sp_str'] = format_percentiles(e['sfeh_sp_1'], e['sfeh_sp_2'],
                                              e['sfeh_sp_3'])

    total = len(entries)
    if total == 0:
        print("No valid entries found for full table.")
        return []

    # Divide entries into column groups
    groups = []
    start_index = 0
    base_size = total // num_output_columns
    remainder = total % num_output_columns

    for i in range(num_output_columns):
        size = base_size + (1 if i < remainder else 0)
        groups.append(entries[start_index:start_index + size])
        start_index += size

    # Prepare LaTeX
    latex_lines = []
    latex_lines.append("\\begin{table*}[htbp]")
    latex_lines.append("\\centering")
    latex_lines.append(
        "\\caption{Full table: Name, RVS [Fe/H], RVS $\\sigma$, SP [Fe/H], SP $\\sigma$.}"
    )
    latex_lines.append("\\label{tab:mytable_full}")

    # 5 columns per group => "lllll"
    tabular_cols = " ".join(["lllll"] * num_output_columns)
    latex_lines.append(f"\\begin{{tabular}}{{{tabular_cols}}}")
    latex_lines.append("\\hline")

    # Build the header
    # 5 columns: Name & [Fe/H]_RVS & sigma_RVS & [Fe/H]_SP & sigma_SP
    col_header = (
        "\\textbf{Name} & "
        "\\textbf{[Fe/H]$_\\mathrm{RVS}$} & \\textbf{$\\sigma$[Fe/H]$_\\mathrm{RVS}$} & "
        "\\textbf{[Fe/H]$_\\mathrm{SP}$} & \\textbf{$\\sigma$[Fe/H]$_\\mathrm{SP}$}"
    )

    header_parts = []
    for _ in range(num_output_columns):
        header_parts.append(col_header)
    latex_lines.append(" & ".join(header_parts) + " \\\\")
    latex_lines.append("\\hline")

    # Number of rows is max of group lengths
    max_rows = max(len(g) for g in groups)
    for i in range(max_rows):
        row_parts = []
        for j in range(num_output_columns):
            if i < len(groups[j]):
                e = groups[j][i]
                row_text = (f"{e['name']} & "
                            f"{e['feh_rvs_str']} & {e['sfeh_rvs_str']} & "
                            f"{e['feh_sp_str']} & {e['sfeh_sp_str']}")
                row_parts.append(row_text)
            else:
                # 5 columns => need 4 ampersands
                row_parts.append("& & & &")
        latex_lines.append(" & ".join(row_parts) + " \\\\")
    latex_lines.append("\\hline")
    latex_lines.append("\\end{tabular}")
    latex_lines.append("\\end{table*}")

    return latex_lines


def main():
    # Decide which table type to produce (default "minimal" or "full")
    if len(sys.argv) > 1:
        table_type = sys.argv[1].lower()
    else:
        table_type = "minimal"  # fallback

    # Location of input file
    input_file = "output/objs.txt"

    # Read the data with astropy.table
    try:
        t = Table.read(input_file, format='ascii')
    except Exception as e:
        print(f"Error reading {input_file} with astropy: {e}")
        sys.exit(1)

    entries = []
    for row in t:
        # 'type', 'name', 'count'
        # 'feh_rvs_1', 'feh_rvs_2', 'feh_rvs_3'
        # 'sfeh_rvs_1', 'sfeh_rvs_2', 'sfeh_rvs_3'
        # 'feh_sp_1', 'feh_sp_2', 'feh_sp_3'
        # 'sfeh_sp_1', 'sfeh_sp_2', 'sfeh_sp_3'
        obj_type = str(row['type'])
        raw_name = str(row['name'])
        count_str = str(row['count'])

        # Decide how to transform name
        if obj_type == "GC":
            name = transform_gc_name(raw_name)
        elif obj_type == "dSph":
            name = transform_dsph_name(raw_name)
        elif obj_type == "Stream":
            name = transform_stream_name(raw_name)
        else:
            name = transform_default_name(raw_name)

        data_dict = {
            'type': obj_type,
            'raw_name': raw_name,
            'name': name,
            'count': count_str,
            # RVS
            'feh_rvs_1': row['feh_rvs_1'],
            'feh_rvs_2': row['feh_rvs_2'],
            'feh_rvs_3': row['feh_rvs_3'],
            'sfeh_rvs_1': row['sfeh_rvs_1'],
            'sfeh_rvs_2': row['sfeh_rvs_2'],
            'sfeh_rvs_3': row['sfeh_rvs_3'],
            # SP
            'feh_sp_1': row['feh_sp_1'],
            'feh_sp_2': row['feh_sp_2'],
            'feh_sp_3': row['feh_sp_3'],
            'sfeh_sp_1': row['sfeh_sp_1'],
            'sfeh_sp_2': row['sfeh_sp_2'],
            'sfeh_sp_3': row['sfeh_sp_3'],
        }
        entries.append(data_dict)

    # Generate the requested table
    if table_type == "minimal":
        latex_lines = create_minimal_table(entries)
    else:
        latex_lines = create_full_table(entries)

    # Print LaTeX lines to stdout
    if latex_lines:
        for line in latex_lines:
            print(line)


if __name__ == "__main__":
    main()
