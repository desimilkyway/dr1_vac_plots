#!/usr/bin/env python3

import re


def transform_gc_name(name):
    """
    If name is of the form NGC_XXXX_M_YY, return M_YY.
    Otherwise, return the original name.
    """
    # Regex to match: starts with NGC_, then some digits, then _M_, then digits
    # Example: "NGC_5904_M_5" -> "M_5"
    pattern = r'^NGC_\d+_M_(\d+)$'
    match = re.match(pattern, name)
    if match:
        ret = f"M_{match.group(1)}"
    else:
        ret = name
    ret = ret.replace('_', ' ')
    print(name, 'x', ret)
    return ret


def transform_dsph_name(n):
    n = n.replace('_', ' ')
    return n


def main():
    input_file = "objs.txt"

    # 1) Read lines from file, store as (object_type, name, number)
    entries = []
    with open(input_file, "r") as f:
        for line in f:
            # Skip empty or comment lines
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) != 3:
                # If line doesn't have exactly 3 parts, skip or handle error
                continue

            obj_type, raw_name, number_str = parts

            # 2) Transform GC name if needed
            if obj_type == "GC":
                name = transform_gc_name(raw_name)
            elif obj_type == "dSph":
                name = transform_dsph_name(raw_name)
            else:
                name = raw_name

            # Convert number to int (optional) or keep as string
            entries.append((obj_type, name, number_str))
    print(entries)
    # 3) Split into three roughly equal groups
    total = len(entries)
    # Example way: group1 gets (total+2)//3, group2 gets (total+1)//3, group3 gets total//3
    # But you can also do a simpler approach:
    size1 = (total + 2) // 3  # a bit “round up” for group1
    size2 = (total + 1) // 3  # a bit “round up” for group2
    # Then size3 is remainder
    size3 = total - size1 - size2

    group1 = entries[:size1]
    group2 = entries[size1:size1 + size2]
    group3 = entries[size1 + size2:]

    # Sanity check
    assert len(group1) == size1
    assert len(group2) == size2
    assert len(group3) == size3

    # 4) Prepare lines of LaTeX
    #    We'll produce a table with at most `max_rows` rows
    max_rows = max(size1, size2, size3)

    latex_lines = []
    latex_lines.append("\\begin{table}[htbp]")
    latex_lines.append("\\centering")
    latex_lines.append(
        "\\caption{Split into three column groups, with NGC\\_XXXX\\_M\\_YY replaced by M\\_YY.}"
    )
    latex_lines.append("\\label{tab:mytable}")
    latex_lines.append("\\begin{tabular}{lll lll lll}")
    latex_lines.append("\\hline")
    latex_lines.append(
        "\\textbf{Object Type} & \\textbf{Name} & \\textbf{Number} &"
        " \\textbf{Object Type} & \\textbf{Name} & \\textbf{Number} &"
        " \\textbf{Object Type} & \\textbf{Name} & \\textbf{Number} \\\\")
    latex_lines.append("\\hline")

    # 5) Build each row using group1[i], group2[i], group3[i]
    for i in range(max_rows):
        # For each group, check if i < length, else use blanks
        if i < len(group1):
            obj1, name1, num1 = group1[i]
        else:
            obj1, name1, num1 = "", "", ""

        if i < len(group2):
            obj2, name2, num2 = group2[i]
        else:
            obj2, name2, num2 = "", "", ""

        if i < len(group3):
            obj3, name3, num3 = group3[i]
        else:
            obj3, name3, num3 = "", "", ""

        # Build LaTeX row
        row_text = (f"{obj1} & {name1} & {num1} & "
                    f"{obj2} & {name2} & {num2} & "
                    f"{obj3} & {name3} & {num3} \\\\")
        latex_lines.append(row_text)

    latex_lines.append("\\hline")
    latex_lines.append("\\end{tabular}")
    latex_lines.append("\\end{table}")

    # 6) Print or write the LaTeX output to stdout
    #    (If you want to write to a file, open a new file and write lines)
    for line in latex_lines:
        print(line)


if __name__ == "__main__":
    main()
