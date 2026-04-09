#!/usr/bin/env python3
"""
Convert a Snakemake --containerize Dockerfile to a Singularity .def file.

Usage:
    python dockerfile_to_singularity.py Dockerfile [output.def]
"""

import argparse
import sys
from pathlib import Path


def convert(dockerfile_path: str, output_path: str = None):
    if dockerfile_path == "-":
        lines = sys.stdin.read().splitlines()
        if output_path is None:
            print("Error: -o/--output is required when reading from stdin")
            sys.exit(1)
    else:
        dockerfile = Path(dockerfile_path)
        if not dockerfile.exists():
            print(f"Error: {dockerfile_path} not found")
            sys.exit(1)
        if output_path is None:
            output_path = str(dockerfile.with_suffix(".def"))
        lines = dockerfile.read_text().splitlines()

    joined_lines = []
    buffer = ""

    for line in lines:
        stripped = line.rstrip()

        if stripped.endswith("\\"):
            buffer += stripped[:-1] + " "
        else:
            buffer += stripped
            joined_lines.append(buffer)
            buffer = ""

    lines = joined_lines

    bootstrap = "docker"
    base_image = None
    local_files = []  # (src, dst) from COPY
    remote_files = []  # (url, dst) from ADD
    mkdirs = []
    conda_envs = []

    for line in lines:
        line = line.strip()

        if line.startswith("FROM "):
            base_image = line[5:].strip()

        elif line.startswith("RUN mkdir -p "):
            mkdir_paths = line[13:].strip().split()
            mkdirs.extend(mkdir_paths)

        elif line.startswith("COPY "):
            parts = line[5:].strip().split()
            if len(parts) == 2:
                local_files.append((parts[0], parts[1]))

        elif line.startswith("ADD "):
            parts = line[4:].strip().split()
            if len(parts) == 2:
                remote_files.append((parts[0], parts[1]))

        elif line.startswith("RUN "):
            cmd = line[4:].strip()  # alles nach RUN

            # Split chained commands
            parts = [c.strip() for c in cmd.split("&&")]

            for part in parts:
                if part.startswith("conda env create"):
                    conda_envs.append(part)

    if base_image is None:
        print("Error: could not find FROM in Dockerfile")
        sys.exit(1)

    out = []
    out.append(f"Bootstrap: {bootstrap}")
    out.append(f"From: {base_image}")
    out.append("")

    # %files section (local files only)
    if local_files:
        out.append("%files")
        for src, dst in local_files:
            out.append(f"    {src} {dst}")
        out.append("")

    # %post section
    out.append("%post")

    # mkdirs
    for cmd in mkdirs:
        out.append(f"    mkdir -p {cmd}")
    out.append("")

    # install build tools (needed for source-compiled packages like whatshap)
    out.append("    apt-get update && apt-get install -y g++ gcc")
    out.append("")

    # download remote files with wget
    if remote_files:
        for url, dst in remote_files:
            out.append(f"    wget -q {url} -O {dst}")
        out.append("")

    # conda env creates
    for cmd in conda_envs:
        out.append(f"    {cmd}")
    out.append("    conda clean --all -y")
    out.append("")

    result = "\n".join(out)
    Path(output_path).write_text(result)
    print(f"Written to {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert a Snakemake --containerize Dockerfile to a Singularity .def file."
    )
    parser.add_argument("dockerfile", help="Path to the input Dockerfile")
    parser.add_argument(
        "-o",
        "--output",
        help="Path to the output .def file (default: input path with .def extension)",
        default=None,
    )
    args = parser.parse_args()
    convert(args.dockerfile, args.output)
