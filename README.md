# compute_poc_flux MATLAB Package

This repository provides a production-ready MATLAB function `compute_poc_flux.m` for computing particulate organic carbon (POC) vertical export fluxes at Zeu and 100 m based on cp660-derived POC profiles, NPP, mixed-layer depth (MLD), and the fraction of microplankton.

It follows GitHub publication standards with function encapsulation, parameter validation, documentation, examples, and optional parallelization.

---

## 🚀 Features

* Full function encapsulation: `results = compute_poc_flux(params)`
* Supports Zeu-depth and 100 m flux estimation
* Robust handling of cp660 → POC conversion
* Preallocated outputs for performance
* Optional parallel computing using `parfor`
* Error handling to ensure full-grid completion
* Fully documented with inline comments

---

## 📁 Function Location

```
compute_poc_flux.m
```

The function is placed in the main directory. All inputs must be provided externally by the user before calling.

---

## 📥 Required Inputs

Inputs are grouped in a structure `params`:

| Field                      | Description                                         |
| -------------------------- | --------------------------------------------------- |
| `Clim_NPP`                 | 122×199×12 monthly NPP (mg C m⁻² d⁻¹ or compatible) |
| `SCSdepth9km`              | bathymetry depth (m)                                |
| `Clim_Zeu2`, `Clim_Zeu3`   | photic depth fields                                 |
| `Clidata_NN_cp660`         | cp660 profiles: 122×199×depth×12                    |
| `depth`                    | depth vector (length = profile layers)              |
| `Clim_MLD`                 | monthly mixed layer depth                           |
| `FM_bbp_Micro`             | microplankton fraction (%)                          |
| `useParallel` *(optional)* | use parfor (default false)                          |

---

## 🧠 Usage Example

```matlab
params.Clim_NPP       = Clim_NPP;
params.SCSdepth9km    = SCSdepth9km;
params.Clim_Zeu2      = Clim_Zeu2;
params.Clim_Zeu3      = Clim_Zeu3;
params.Clidata_NN_cp660 = Clidata_NN_cp660;
params.depth          = depth;
params.Clim_MLD       = Clim_MLD;
params.FM_bbp_Micro   = FM_bbp_Micro;
params.useParallel    = false;   % optional

results = compute_poc_flux(params);
```

Output fields include:

* `Flux_Zeu_nn`
* `Flux_100_nn`
* `DPOC`
* `DPOC100`
* `IPOC_Zeu_NN`, `IPOC_100_NN`
* corresponding `_1` fields for previous-month comparisons

---

## 📤 Output Structure

Example:

```matlab
results.Flux_Zeu_nn   % 122×199×12
results.DPOC100       % 122×199×12
```

All outputs use NaN for invalid ocean/non-ocean cells.

---

## 🧪 Unit Testing

A test file `test_compute_poc_flux.m` is provided to verify:

* dimension consistency
* NaN handling
* no-crash execution

---

## 📄 License

MIT License — free for academic and commercial use.

---

## 🙋 Questions?

If you encounter performance issues or want a vectorized version, feel free to open an Issue.
