#!/usr/bin/env python3
"""Prepare a four-node FRP-grid interface UEL input for Abaqus.

The script replaces a user-selected placeholder *ELEMENT block with U1 elements,
computes the node-influence weight omega_n, assigns an orthonormal local frame,
and writes element-specific *UEL PROPERTY blocks.

The input placeholder may have four or more nodes. For an 8-node cohesive
placeholder, define ``node_order`` in the JSON configuration to select and order
the four nodes used by the line-interface UEL.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


FLOAT_KEYS = [
    "branch_width",
    "k_uu0", "t_u_max", "delta_u_p", "delta_u_f", "t_u_res",
    "k_vc", "k_vt", "delta_v0", "k_uv0", "alpha_c", "alpha_t",
    "beta_k", "beta_tau", "lambda_r", "k_ww",
    "eta_k_u", "eta_k_v", "eta_k_uv", "eta_t_u", "eta_delta_v",
    "eta_beta_k", "eta_beta_tau",
]


@dataclass(frozen=True)
class Element:
    label: int
    connectivity: tuple[int, ...]


def parse_keyword(line: str) -> tuple[str, dict[str, str]]:
    parts = [part.strip() for part in line.strip().split(",")]
    keyword = parts[0].upper()
    options: dict[str, str] = {}
    for item in parts[1:]:
        if "=" in item:
            key, value = item.split("=", 1)
            options[key.strip().upper()] = value.strip()
        elif item:
            options[item.upper()] = ""
    return keyword, options


def split_numbers(line: str) -> list[str]:
    return [x.strip() for x in line.split(",") if x.strip()]


def read_nodes(lines: Sequence[str]) -> dict[int, tuple[float, float, float]]:
    nodes: dict[int, tuple[float, float, float]] = {}
    i = 0
    while i < len(lines):
        if lines[i].lstrip().startswith("*"):
            keyword, _ = parse_keyword(lines[i])
            if keyword == "*NODE":
                i += 1
                while i < len(lines) and not lines[i].lstrip().startswith("*"):
                    if lines[i].strip() and not lines[i].lstrip().startswith("**"):
                        fields = split_numbers(lines[i])
                        if len(fields) >= 4:
                            nodes[int(fields[0])] = (
                                float(fields[1]), float(fields[2]), float(fields[3])
                            )
                    i += 1
                continue
        i += 1
    if not nodes:
        raise ValueError("No *NODE data were found in the input file.")
    return nodes


def find_target_blocks(lines: Sequence[str], target_elset: str) -> list[tuple[int, int, list[Element]]]:
    target = target_elset.upper()
    blocks: list[tuple[int, int, list[Element]]] = []
    i = 0
    while i < len(lines):
        if lines[i].lstrip().startswith("*"):
            keyword, options = parse_keyword(lines[i])
            if keyword == "*ELEMENT" and options.get("ELSET", "").upper() == target:
                start = i
                i += 1
                elems: list[Element] = []
                while i < len(lines) and not lines[i].lstrip().startswith("*"):
                    if lines[i].strip() and not lines[i].lstrip().startswith("**"):
                        fields = split_numbers(lines[i])
                        elems.append(Element(int(fields[0]), tuple(map(int, fields[1:]))))
                    i += 1
                blocks.append((start, i, elems))
                continue
        i += 1
    if not blocks:
        raise ValueError(f"No *ELEMENT block with ELSET={target_elset!r} was found.")
    return blocks


def vec_add(a: Sequence[float], b: Sequence[float]) -> tuple[float, float, float]:
    return (a[0]+b[0], a[1]+b[1], a[2]+b[2])


def vec_sub(a: Sequence[float], b: Sequence[float]) -> tuple[float, float, float]:
    return (a[0]-b[0], a[1]-b[1], a[2]-b[2])


def vec_scale(a: Sequence[float], s: float) -> tuple[float, float, float]:
    return (a[0]*s, a[1]*s, a[2]*s)


def dot(a: Sequence[float], b: Sequence[float]) -> float:
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]


def cross(a: Sequence[float], b: Sequence[float]) -> tuple[float, float, float]:
    return (
        a[1]*b[2]-a[2]*b[1],
        a[2]*b[0]-a[0]*b[2],
        a[0]*b[1]-a[1]*b[0],
    )


def norm(a: Sequence[float]) -> float:
    return math.sqrt(dot(a, a))


def unit(a: Sequence[float], name: str) -> tuple[float, float, float]:
    n = norm(a)
    if n <= 1.0e-12:
        raise ValueError(f"Cannot normalize the {name}; its magnitude is zero.")
    return vec_scale(a, 1.0/n)


def orthonormal_frame(branch: Sequence[float], normal: Sequence[float]) -> tuple[tuple[float,float,float], ...]:
    eu = unit(branch, "branch direction")
    ev0 = vec_sub(normal, vec_scale(eu, dot(normal, eu)))
    ev = unit(ev0, "normal direction after orthogonalization")
    ew = unit(cross(eu, ev), "transverse direction")
    ev = unit(cross(ew, eu), "corrected normal direction")
    return eu, ev, ew


def midpoint(a: Sequence[float], b: Sequence[float]) -> tuple[float, float, float]:
    return vec_scale(vec_add(a, b), 0.5)


def overlap_length(x_minus: float, x_plus: float, center: float, influence_length: float) -> float:
    left = center - influence_length/2.0
    right = center + influence_length/2.0
    return max(0.0, min(x_plus, right) - max(x_minus, left))


def format_property_lines(values: Sequence[float], per_line: int = 8) -> list[str]:
    out: list[str] = []
    for i in range(0, len(values), per_line):
        chunk = values[i:i+per_line]
        out.append(", ".join(f"{v:.12g}" for v in chunk) + "\n")
    return out


def load_config(path: Path) -> dict:
    cfg = json.loads(path.read_text(encoding="utf-8"))
    required = ["target_elset", "node_order", "branch_direction", "normal_direction",
                "origin", "node_centers", "node_influence_length", "parameters"]
    missing = [key for key in required if key not in cfg]
    if missing:
        raise ValueError(f"Missing configuration entries: {', '.join(missing)}")
    params = cfg["parameters"]
    missing_params = [key for key in FLOAT_KEYS if key not in params]
    if missing_params:
        raise ValueError(f"Missing parameter entries: {', '.join(missing_params)}")
    return cfg


def selected_connectivity(element: Element, node_order: Sequence[int]) -> tuple[int, int, int, int]:
    try:
        selected = tuple(element.connectivity[index-1] for index in node_order)
    except IndexError as exc:
        raise ValueError(
            f"Element {element.label} has {len(element.connectivity)} nodes, "
            f"but node_order={node_order} is incompatible."
        ) from exc
    if len(selected) != 4:
        raise ValueError("node_order must contain exactly four local node positions.")
    return selected  # type: ignore[return-value]


def process(input_path: Path, config_path: Path, output_path: Path) -> None:
    lines = input_path.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)
    nodes = read_nodes(lines)
    cfg = load_config(config_path)
    blocks = find_target_blocks(lines, cfg["target_elset"])

    node_order = [int(x) for x in cfg["node_order"]]
    eu_ref, ev, ew = orthonormal_frame(cfg["branch_direction"], cfg["normal_direction"])
    origin = tuple(map(float, cfg["origin"]))
    node_centers = [float(x) for x in cfg["node_centers"]]
    ln = float(cfg["node_influence_length"])
    params = {key: float(cfg["parameters"][key]) for key in FLOAT_KEYS}

    transformed: list[tuple[Element, tuple[int,int,int,int], float, tuple[float,float,float], tuple[float,float,float], tuple[float,float,float]]] = []
    for _, _, elems in blocks:
        for elem in elems:
            conn = selected_connectivity(elem, node_order)
            xyz = [nodes[n] for n in conn]
            start = midpoint(xyz[0], xyz[1])
            end = midpoint(xyz[3], xyz[2])
            e_vec = vec_sub(end, start)
            eu = unit(e_vec, f"element {elem.label} branch axis")
            if dot(eu, eu_ref) < 0.0:
                eu = vec_scale(eu, -1.0)
            # Rebuild an element-specific orthonormal frame using the requested normal.
            eu, ev_e, ew_e = orthonormal_frame(eu, ev)
            x_start = dot(vec_sub(start, origin), eu_ref)
            x_end = dot(vec_sub(end, origin), eu_ref)
            x_minus, x_plus = sorted((x_start, x_end))
            le = x_plus - x_minus
            if le <= 1.0e-12:
                raise ValueError(f"Element {elem.label} has zero projected length.")
            max_overlap = max((overlap_length(x_minus, x_plus, c, ln) for c in node_centers), default=0.0)
            omega = max(0.0, min(1.0, max_overlap/le))
            transformed.append((elem, conn, omega, eu, ev_e, ew_e))

    skip_ranges = {(start, end) for start, end, _ in blocks}
    out: list[str] = []
    i = 0
    inserted = False
    first_start = min(start for start, _, _ in blocks)
    while i < len(lines):
        if i == first_start and not inserted:
            out.extend([
                "** ----------------------------------------------------------------\n",
                "** Four-node FRP-grid interface user element generated by\n",
                "** preprocessing/frp_grid_preprocess.py\n",
                "*USER ELEMENT, TYPE=U1, NODES=4, COORDINATES=3, PROPERTIES=33, VARIABLES=18\n",
                "1, 2, 3\n",
            ])
            for elem, conn, omega, eu, ev_e, ew_e in transformed:
                elset = f"UEL_E_{elem.label}"
                out.append(f"*ELEMENT, TYPE=U1, ELSET={elset}\n")
                out.append(f"{elem.label}, {conn[0]}, {conn[1]}, {conn[2]}, {conn[3]}\n")
                out.append(f"*UEL PROPERTY, ELSET={elset}\n")
                p = params
                values = [
                    p["branch_width"],
                    p["k_uu0"], p["t_u_max"], p["delta_u_p"], p["delta_u_f"], p["t_u_res"],
                    p["k_vc"], p["k_vt"], p["delta_v0"], p["k_uv0"], p["alpha_c"], p["alpha_t"],
                    p["beta_k"], p["beta_tau"], p["lambda_r"], p["k_ww"],
                    p["eta_k_u"], p["eta_k_v"], p["eta_k_uv"], p["eta_t_u"], p["eta_delta_v"],
                    p["eta_beta_k"], p["eta_beta_tau"], omega,
                    *eu, *ev_e, *ew_e,
                ]
                out.extend(format_property_lines(values))
            labels = [str(item[0].label) for item in transformed]
            out.append(f"*ELSET, ELSET={cfg['target_elset']}\n")
            for j in range(0, len(labels), 16):
                out.append(", ".join(labels[j:j+16]) + "\n")
            out.append("** ----------------------------------------------------------------\n")
            inserted = True

        matched_range = next(((s,e) for s,e in skip_ranges if i == s), None)
        if matched_range:
            i = matched_range[1]
            continue
        out.append(lines[i])
        i += 1

    output_path.write_text("".join(out), encoding="utf-8")
    print(f"Wrote {output_path}")
    print(f"Converted {len(transformed)} interface elements.")
    omegas = [item[2] for item in transformed]
    print(f"omega_n range: {min(omegas):.6g} to {max(omegas):.6g}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Original Abaqus .inp file")
    parser.add_argument("config", type=Path, help="JSON configuration file")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Processed Abaqus .inp file")
    args = parser.parse_args()
    try:
        process(args.input, args.config, args.output)
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
