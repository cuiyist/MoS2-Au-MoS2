#!/usr/bin/env python3
"""Pressure-scale models for MoS2-Au-MoS2 confinement.

All public functions accept SI units. The command-line interface accepts nm and
prints human-readable pressures or generates a CSV parameter sweep.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

EV_TO_J = 1.602_176_634e-19
GPA = 1e9
NM = 1e-9


@dataclass(frozen=True)
class ModelParameters:
    """Material parameters used by the two comparison models."""

    hamaker_empty_j: float = 0.70 * EV_TO_J
    c11_pa: float = 238.0 * GPA
    c33_pa: float = 52.0 * GPA
    c13_pa: float = 23.0 * GPA
    c44_pa: float = 18.6 * GPA


DEFAULTS = ModelParameters()


def indentation_modulus(params: ModelParameters = DEFAULTS) -> float:
    """Return the longitudinal indentation modulus of transversely isotropic MoS2."""

    first = (params.c11_pa * params.c33_pa - params.c13_pa**2) / params.c11_pa
    second_inverse = (
        1.0 / params.c44_pa
        + 2.0 / (math.sqrt(params.c11_pa * params.c33_pa) + params.c13_pa)
    )
    return 2.0 * math.sqrt(first / second_inverse)


def finite_slab_kernel(gap_m: float, slab_thickness_m: float) -> float:
    """Return the finite-slab geometric kernel in m^-3."""

    if gap_m <= 0 or slab_thickness_m <= 0:
        raise ValueError("gap and slab thickness must be positive")
    return (
        gap_m**-3
        - 2.0 * (gap_m + slab_thickness_m) ** -3
        + (gap_m + 2.0 * slab_thickness_m) ** -3
    )


def hamaker_pressure(
    gap_m: float,
    slab_thickness_m: float,
    hamaker_j: float,
) -> float:
    """Return the magnitude of the attractive finite-slab Hamaker pressure in Pa."""

    if hamaker_j <= 0:
        raise ValueError("Hamaker constant must be positive")
    return hamaker_j * finite_slab_kernel(gap_m, slab_thickness_m) / (6.0 * math.pi)


def half_space_hamaker_pressure(gap_m: float, hamaker_j: float) -> float:
    """Return the magnitude of the two-half-space Hamaker pressure in Pa."""

    if gap_m <= 0 or hamaker_j <= 0:
        raise ValueError("gap and Hamaker constant must be positive")
    return hamaker_j / (6.0 * math.pi * gap_m**3)


def free_standing_flat_punch_pressure(
    imposed_displacement_m: float,
    disc_diameter_m: float,
    modulus_pa: float,
) -> float:
    """Return the mean pressure for two compliant, free-standing MoS2 flakes.

    The two identical flakes share the imposed relative displacement equally,
    so each flake accommodates half of the Au-disc thickness.
    """

    if imposed_displacement_m <= 0 or disc_diameter_m <= 0 or modulus_pa <= 0:
        raise ValueError("displacement, diameter, and modulus must be positive")
    radius_m = disc_diameter_m / 2.0
    return modulus_pa * imposed_displacement_m / (math.pi * radius_m)


def evaluate_geometry(
    mos2_thickness_nm: float,
    gap_nm: float,
    diameter_nm: float,
    params: ModelParameters = DEFAULTS,
) -> dict[str, float | bool]:
    """Evaluate all comparison models for one geometry."""

    t = mos2_thickness_nm * NM
    d = gap_nm * NM
    diameter = diameter_nm * NM
    modulus = indentation_modulus(params)
    free_standing = free_standing_flat_punch_pressure(d, diameter, modulus)
    empty = hamaker_pressure(d, t, params.hamaker_empty_j)
    half_space = half_space_hamaker_pressure(d, params.hamaker_empty_j)
    radius_nm = diameter_nm / 2.0

    return {
        "mos2_thickness_nm": mos2_thickness_nm,
        "gap_nm": gap_nm,
        "diameter_nm": diameter_nm,
        "indentation_modulus_gpa": modulus / GPA,
        "hamaker_empty_mpa": empty / 1e6,
        "hamaker_fraction_of_half_space": empty / half_space,
        "flat_punch_free_standing_gpa": free_standing / GPA,
        "t_over_a": mos2_thickness_nm / radius_nm,
        "half_space_valid": mos2_thickness_nm / radius_nm >= 2.0,
    }


def inclusive_values(start: float, stop: float, step: float) -> Iterable[float]:
    """Yield decimal grid values including the endpoint within rounding tolerance."""

    if step <= 0 or stop < start:
        raise ValueError("sweep bounds require stop >= start and step > 0")
    count = int(round((stop - start) / step))
    for index in range(count + 1):
        yield round(start + index * step, 10)


def write_sweep(path: Path, params: ModelParameters = DEFAULTS) -> int:
    """Write the default thickness-gap-diameter sweep and return its row count."""

    fieldnames = list(evaluate_geometry(50.0, 2.5, 50.0, params).keys())
    rows = 0
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for thickness in inclusive_values(10.0, 100.0, 5.0):
            for gap in inclusive_values(1.0, 8.0, 0.1):
                for diameter in (30.0, 50.0, 75.0, 100.0):
                    writer.writerow(evaluate_geometry(thickness, gap, diameter, params))
                    rows += 1
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mos2-thickness", type=float, default=50.0, help="MoS2 thickness in nm")
    parser.add_argument("--gap", type=float, default=2.5, help="gap/Au thickness in nm")
    parser.add_argument("--diameter", type=float, default=50.0, help="Au-disc diameter in nm")
    parser.add_argument("--sweep", action="store_true", help="write the default parameter sweep to CSV")
    parser.add_argument("--output", type=Path, default=Path("pressure_sweep.csv"), help="CSV output path")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.sweep:
        rows = write_sweep(args.output)
        print(f"Wrote {rows} rows to {args.output}")
        return

    result = evaluate_geometry(args.mos2_thickness, args.gap, args.diameter)
    print(
        f"Geometry: t={args.mos2_thickness:g} nm, d={args.gap:g} nm, "
        f"D={args.diameter:g} nm"
    )
    print(f"MoS2 indentation modulus: {result['indentation_modulus_gpa']:.3f} GPa")
    print(f"Empty-gap Hamaker pressure: {result['hamaker_empty_mpa']:.6f} MPa")
    print(f"Free-standing finite-disc elastic contact: {result['flat_punch_free_standing_gpa']:.6f} GPa")
    validity = "within" if result["half_space_valid"] else "outside"
    print(f"t/a={result['t_over_a']:.3f}: {validity} the t/a >= 2 half-space criterion")


if __name__ == "__main__":
    main()
