import os
import sys
import pandas as pd


if "snakemake" in globals():
    sys.stderr = open(snakemake.log[0], "w")


def _normalize_column_names(frame):
    frame = frame.copy()
    frame.columns = [str(col).strip().lower() for col in frame.columns]
    return frame


def _resolve_path(inputs, candidates):
    if isinstance(inputs, dict):
        for candidate in candidates:
            if candidate in inputs:
                return str(inputs[candidate])
    if isinstance(inputs, (list, tuple)):
        if inputs:
            return str(inputs[0])
    if isinstance(inputs, str):
        return str(inputs)
    return None


def _read_table(path):
    if not path or not os.path.exists(path):
        raise FileNotFoundError(f"Expected input file not found: {path}")
    return pd.read_csv(path, sep=None, engine="python")


def _check_duplicates(frame, label):
    errors = []
    if frame.empty:
        return errors

    first_column = frame.columns[0]
    values = frame[first_column].dropna().astype(str).str.strip()
    values = values[values != ""]
    duplicates = values[values.duplicated(keep=False)]

    if not duplicates.empty:
        duplicate_values = sorted(set(duplicates.tolist()))
        errors.append(
            f"Duplicate sample names were found in {label} in column '{first_column}': "
            + ", ".join(duplicate_values)
        )

    return errors


def _check_empty_values(frame, label):
    errors = []
    if frame.empty:
        return errors

    blank_mask = frame.apply(lambda col: col.astype(str).str.strip() == "")
    if blank_mask.any().any():
        bad_rows = frame[blank_mask.any(axis=1)]
        errors.append(
            f"{label} contains empty values in the following rows:\n"
            + bad_rows.to_string(index=False)
        )

    empty_columns = [
        col for col in frame.columns if frame[col].astype(str).str.strip().eq("").all()
    ]
    if empty_columns:
        errors.append(f"{label} contains empty columns: " + ", ".join(empty_columns))

    return errors


def main():
    inputs = snakemake.input if "snakemake" in globals() else {}
    sample_tsv_path = _resolve_path(
        inputs, ["sample_tsv", "sample.tsv", "metadata", "sample_metadata"]
    )
    sample_info_path = _resolve_path(
        inputs, ["sample_info", "sample_info.txt", "sample_info_file"]
    )

    if not sample_tsv_path:
        sample_tsv_path = "config/pep/sample.tsv"
    if not sample_info_path:
        sample_info_path = "config/pep/sample_info.txt"

    sample_tsv = _normalize_column_names(_read_table(sample_tsv_path))
    sample_info = _normalize_column_names(_read_table(sample_info_path))

    errors = []
    errors.extend(_check_duplicates(sample_tsv, "sample.tsv"))
    errors.extend(_check_duplicates(sample_info, "sample_info.txt"))
    errors.extend(_check_empty_values(sample_info, "sample_info.txt"))

    if errors:
        raise ValueError("\n".join(errors))


if __name__ == "__main__":
    main()
