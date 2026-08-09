# Confinement-pressure estimates for MoS₂–Au–MoS₂ nanogaps

This directory contains a reproducible estimate of the normal-pressure scales relevant to an Au nanodisc confined between two MoS₂ flakes during annealing. The calculation uses the experimentally relevant **50 nm mean Au-disc diameter**, approximately **50 nm MoS₂ thickness per flake**, and a **2.4–5.3 nm disc/gap thickness**.

The models answer different physical questions and should not be added together:

1. **Finite-slab Hamaker pressure, empty gap** — direct MoS₂–MoS₂ attraction across vacuum/Ar.
2. **Finite-slab Hamaker pressure, Au-filled gap** — an illustrative effective Lifshitz/Hamaker term across Au.
3. **Finite-disc elastic contact** — the reaction stress needed for transversely anisotropic MoS₂ to conform to a 50 nm Au nanodisc.

![Schematics of the three pressure models](figures/model-schematics.svg)

## Main result

For two 50 nm MoS₂ flakes and a 50 nm diameter Au disc:

| Gap / Au thickness | Empty-gap Hamaker | Au-filled Hamaker | Adhesion scale, Γ/d | Flat punch, two compliant slabs | Flat punch, supported lower slab |
|---:|---:|---:|---:|---:|---:|
| 2.4 nm | 0.430 MPa | 0.357 MPa | 0.201 GPa | 1.646 GPa | 3.291 GPa |
| 2.5 nm | 0.381 MPa | 0.315 MPa | 0.193 GPa | 1.714 GPa | 3.428 GPa |
| 5.3 nm | 0.040 MPa | 0.033 MPa | 0.091 GPa | 3.634 GPa | 7.268 GPa |

The direct Hamaker attraction is therefore **sub-MPa** over the observed thickness range. The cohesion-energy scale is approximately **0.09–0.20 GPa**, while a fully conformal, linear-elastic flat-punch model gives **1.65–3.63 GPa** when the two flakes share the deformation and **3.29–7.27 GPa** when the lower flake is treated as supported.

These values are not a single predicted hydrostatic pressure. The Hamaker result is a direct long-range traction; Γ/d is an energy-density scale; and the flat-punch result is an elastic reaction-stress scale. At the annealing temperature, interfacial sliding, local debonding, MoS₂ bending, Au plasticity, diffusion, creep, and imperfect contact can all relax the elastic stress. The defensible conclusion is that the architecture can generate an **effective confinement scale from roughly 0.1 GPa to several GPa**, with the GPa values interpreted as ideal-elastic upper scales rather than a uniquely determined sustained pressure. This range is consistent with the broader pressure scales reported for nanoscale confinement in van der Waals heterostructures.<sup>1</sup>

## Interactive parameter sweep

Open [`docs/index.html`](docs/index.html) in a browser to explore:

- MoS₂ thickness: 10–100 nm
- gap/Au thickness: 1–8 nm
- Au-disc diameter: 30–100 nm

![Static preview of the interactive pressure explorer](figures/interactive-preview.png)

From this directory, a simple local server avoids browser restrictions on local files:

```bash
python3 -m http.server 8000 --directory docs
```

Then visit `http://localhost:8000`.

## Model 1: finite MoS₂ slabs across an empty gap

For two identical slabs of thickness \(t\), separated by distance \(d\), the non-retarded pairwise-additive finite-slab result is

$$
P_{\mathrm{H}}(d,t)=\frac{A}{6\pi}
\left[
\frac{1}{d^3}-\frac{2}{(d+t)^3}+\frac{1}{(d+2t)^3}
\right],
$$

where the table uses \(A_{\mathrm{empty}}=0.70\ \mathrm{eV}\), within the range measured for few-layer and multilayer MoS₂ interfaces.<sup>2</sup> The sign is attractive; the repository reports its magnitude.

The pressure approaches the half-space result \(A/(6\pi d^3)\) once each flake is several gap widths thick. At \(d=5.3\) nm, the finite-thickness factor relative to two half-spaces is:

| MoS₂ thickness | Fraction of half-space pressure |
|---:|---:|
| 10 nm | 92.606% |
| 20 nm | 98.322% |
| 30 nm | 99.377% |
| 50 nm | 99.837% |
| 100 nm | 99.976% |

At \(d=2.4\) nm and \(t=50\) nm, the factor is 99.982%. Thus the Hamaker result is already effectively thickness-saturated for the experimental geometry and is dominated by its \(d^{-3}\) dependence.

## Model 2: finite MoS₂ slabs with Au in the gap

The same geometry is evaluated with an illustrative effective constant
\(A_{\mathrm{MoS_2|Au|MoS_2}}=0.58\ \mathrm{eV}\), informed by elemental-metal Hamaker calculations and first-principles descriptions of the Au–MoS₂ interface:<sup>3,4</sup>

$$
P_{\mathrm{H,Au}}(d,t)=\frac{A_{\mathrm{eff}}}{6\pi}
\left[
\frac{1}{d^3}-\frac{2}{(d+t)^3}+\frac{1}{(d+2t)^3}
\right].
$$

This is a screening estimate, not a full multilayer Lifshitz calculation. A quantitatively rigorous value requires the frequency-dependent dielectric response \(\varepsilon(i\xi)\) of MoS₂ and Au at the annealing temperature, including MoS₂ anisotropy and the actual interfacial medium. The chosen effective constant should therefore be varied; it is exposed by the Python API.

## Model 3: a finite Au disc as an anisotropic flat punch

The Au feature is not a 200 nm plate. It is treated as a circular disc with mean diameter \(D=50\) nm and radius \(a=25\) nm. For loading normal to the basal plane, a transversely anisotropic indentation modulus is calculated using the anisotropic half-space contact treatment of Vlassak and Nix:<sup>5</sup>

$$
M_L=2\sqrt{
\frac{C_{11}C_{33}-C_{13}^2}{C_{11}}
\left[
\frac{1}{C_{44}}+
\frac{2}{\sqrt{C_{11}C_{33}}+C_{13}}
\right]^{-1}}
$$

using reported elastic constants \(C_{11}=238\), \(C_{33}=52\), \(C_{13}=23\), and \(C_{44}=18.6\) GPa.<sup>6</sup> This gives \(M_L=53.85\) GPa.

For a flat circular contact, the mean pressure scales as

$$
\bar P_{\mathrm{two}}=\frac{M_L d}{\pi a},
\qquad
\bar P_{\mathrm{supported}}=\frac{2M_L d}{\pi a}.
$$

The first expression lets two identical MoS₂ bodies share the total relative displacement \(d\); the second places the full displacement in one MoS₂ body above a supported lower side. Treating the full nanodisc thickness as the imposed relative displacement is an idealized, fully conformal limit, so these are reaction-stress scales rather than direct van der Waals tractions.

The diameter sensitivity at \(d=2.5\) nm is:

| Au diameter | Two compliant slabs | Supported lower slab |
|---:|---:|---:|
| 30 nm | 2.86 GPa | 5.71 GPa |
| 50 nm | 1.71 GPa | 3.43 GPa |
| 100 nm | 0.86 GPa | 1.71 GPa |

The half-space approximation is most credible when \(t/a\gtrsim2\), equivalently \(t\gtrsim D\). The experimental point \(t\approx D\approx50\) nm is therefore near the lower edge of that regime. For thinner MoS₂, the script flags the result as an extrapolation; a finite-layer anisotropic contact calculation or finite-element model is needed to predict the thickness correction. The script deliberately does not hide this limitation behind an unsupported empirical correction.

## Adhesion-energy cross-check

A separate energy-density estimate is

$$
P_\Gamma=\frac{\Gamma}{d},
$$

with the experimentally measured MoS₂ cohesion energy \(\Gamma=0.482\pm0.032\ \mathrm{J\,m^{-2}}\).<sup>7</sup> Using the central value produces 0.201 GPa at 2.4 nm and 0.091 GPa at 5.3 nm. It is useful as a work-per-volume scale, but it is not the local derivative of the interaction potential and should not be added to the Hamaker or elastic-contact pressure.

## Reproduce the calculations

The code uses only the Python standard library.

```bash
# One geometry (defaults: t=50 nm, d=2.5 nm, D=50 nm)
python3 pressure_models.py

# A selected geometry
python3 pressure_models.py --mos2-thickness 50 --gap 2.4 --diameter 50

# Write a parameter sweep to CSV
python3 pressure_models.py --sweep --output pressure_sweep.csv

# Run regression tests
python3 -m unittest -v
```

The CSV sweep includes all three models, the Γ/d scale, and a Boolean flag for the half-space validity criterion.

## Reviewer-response text

> To provide a quantitative estimate of the confinement pressure, we compared three limiting descriptions for the experimental geometry: two approximately 50 nm thick MoS₂ flakes confining an Au nanodisc with a mean lateral diameter of approximately 50 nm and a thickness of 2.4–5.3 nm. A finite-slab Hamaker model gives only 0.43–0.040 MPa for direct MoS₂–MoS₂ attraction across an empty gap and approximately 0.36–0.033 MPa when an illustrative effective Au-filled Hamaker constant is used. An independent cohesion-energy scale, Γ/d, is 0.20–0.091 GPa. To account for the finite Au nanodisc, we also treated the disc as a circular flat punch acting on transversely anisotropic MoS₂. Using an out-of-plane indentation modulus of 53.9 GPa gives a mean elastic reaction stress of 1.65–3.63 GPa when the two flakes share the deformation, or 3.29–7.27 GPa when the lower flake is treated as supported. The approximately 50 nm flake thickness lies at the lower edge of the half-space validity criterion for a 50 nm disc; thinner flakes require a finite-layer anisotropic calculation. These estimates refer to distinct physical quantities and are not summed. Because sliding, bending, imperfect contact, Au plasticity and diffusion can relax stress during annealing, we interpret the result as evidence for an effective confinement scale of roughly 0.1 GPa to several GPa, rather than as a unique hydrostatic pressure. This quantitative range connects the observed Au restructuring to confined nanoscale synthesis while making the model assumptions explicit.

## Assumptions and limitations

- The Hamaker models are non-retarded, pairwise-additive continuum approximations.
- The Au-filled effective Hamaker constant is illustrative. Full Lifshitz theory with \(\varepsilon(i\xi)\) is needed for a high-accuracy value.
- The flat-punch model assumes linear elasticity, full contact, a circular disc, and loading normal to the MoS₂ basal plane.
- MoS₂ is represented by a transversely anisotropic half-space where \(t/a\gtrsim2\). The elastic pressure shown at smaller \(t/a\) is an extrapolation.
- The calculation gives a mean contact pressure; edge stresses are nonuniform and singular in an ideal flat-punch solution.
- The analysis does not model thermally activated Au flow, interfacial slip, defects, or time-dependent relaxation.

## References

1. Khestanova, E., Guinea, F., Fumagalli, L., Geim, A. K. & Grigorieva, I. V. Universal shape and pressure inside bubbles appearing in van der Waals heterostructures. *Nat. Commun.* **7**, 12587 (2016). [https://doi.org/10.1038/ncomms12587](https://doi.org/10.1038/ncomms12587)
2. Gisbert, V. G. & Garcia, R. Fast and high-resolution mapping of van der Waals forces of 2D materials interfaces with bimodal AFM. *Nanoscale* **15**, 19196–19202 (2023). [https://doi.org/10.1039/D3NR05274E](https://doi.org/10.1039/D3NR05274E)
3. Tolias, P. Non-retarded room temperature Hamaker constants between elemental metals. *Surf. Sci.* **700**, 121652 (2020). [https://doi.org/10.1016/j.susc.2020.121652](https://doi.org/10.1016/j.susc.2020.121652)
4. Boschetto, G., Carapezzi, S. & Todri-Sanial, A. Non-volatile resistive switching mechanism in single-layer MoS₂ memristors: insights from *ab initio* modelling of Au and MoS₂ interfaces. *Nanoscale Adv.* **5**, 4203–4212 (2023). [https://doi.org/10.1039/D3NA00045A](https://doi.org/10.1039/D3NA00045A)
5. Vlassak, J. J. & Nix, W. D. Indentation modulus of elastically anisotropic half spaces. *Philos. Mag. A* **67**, 1045–1056 (1993). [https://doi.org/10.1080/01418619308224756](https://doi.org/10.1080/01418619308224756)
6. Yengejeh, S. I., Liu, J., Kazemi, S. A., Wen, W. & Wang, Y. Effect of structural phases on mechanical properties of molybdenum disulfide. *ACS Omega* **5**, 5994–6002 (2020). [https://doi.org/10.1021/acsomega.9b04360](https://doi.org/10.1021/acsomega.9b04360)
7. Rokni, H. & Lu, W. Direct measurements of interfacial adhesion in 2D materials and van der Waals heterostructures in ambient air. *Nat. Commun.* **11**, 5607 (2020). [https://doi.org/10.1038/s41467-020-19411-7](https://doi.org/10.1038/s41467-020-19411-7)
