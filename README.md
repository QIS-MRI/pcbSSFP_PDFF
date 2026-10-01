# Phase-cycled bSSFP PDFF mapping

MATLAB code for **proton-density fat fraction (PDFF)** mapping from **phase-cycled balanced SSFP**, with **ΔB0 estimation** so that fat–water swaps from field inhomogeneity can be corrected before SPARCQ reconstruction.

The method is described in:

> Acikgoz BC, Mackowiak ALC, Bongiolatti-Rossi GMC, et al. *Fat fraction mapping in the presence of magnetic field inhomogeneities with phase-cycled bSSFP at 3T.*

Typical acquisition (as in the manuscript): TR/TE ≈ 3.4/1.7 ms, FA 35°, 3T, complex profiles `[x y nPC]`.

## Pipeline

1. Phase-correct the measured profiles (remove coil / eddy / systematic phase).
2. Fit residuals over a discretized ΔB0 range, then pick a spatially smooth field map with **graph cuts**.
3. Correct profiles for ΔB0 (Fourier circular shift + accumulated phase).
4. Identify very low / very high PDFF with a phase lookup table.
5. Run **SPARCQ** dictionary matching on the remaining voxels.

## How to run

MATLAB R2023b or later. Toolboxes: Optimization (`lsqnonneg`, `lsqlin`), Image Processing (`imresize`, `bwareaopen`), Parallel Computing (`parfor`).

**One-time:** copy a working `max_flow_mex` (Windows: `max_flow_mex.mexw64`) into

```
third_party/matlab_bgl/private/
```

Do not add that `private` folder to the MATLAB path.

1. Put a `.mat` with complex `profiles` `[x y nPC]` into `example data/`.
2. Edit **USER SETTINGS** in `prepare_data.m` (filename, `type` = `Liver` | `Knee` | `Calimetrix` | `Butters`, TR, FA, …) and run it. This writes a `dataset` struct.
3. Set `prepared_name` in `run_sparcq.m` to that filename and run it.
4. Draw a background mask (slider + **Continue**).
5. Optional: in `run_sparcq.m`, set `plot_vials = true` to click phantom vials against the commercial (Calimetrix) or custom (Butters) ground-truth lists.

Prepared `dataset` structs and optional `*_results.mat` stay in `example data/`. Sample commercial-phantom data (shim settings and RF excitation angle) are included there. Code: https://github.com/QIS-MRI/pcbSSFP_PDFF

```
prepare_data.m          pack scan + algorithm parameters
run_sparcq.m            mask GUI + reconstruction + figures
lib/                    SPARCQ pipeline
third_party/hernando/   graph-cut (Hernando)
third_party/matlab_bgl/ max_flow (paste MEX in private/)
example data/           input and output .mat files
```

## Citations

If you use this code, please cite the SPARCQ PDFF papers and, kindly, the graph-cut method and the ISMRM Fat-Water Toolbox that this reconstruction relies on (as in the manuscript, refs 10 and 28):

**Graph-cut field mapping**

Hernando D, Kellman P, Haldar JP, Liang ZP. Robust water/fat separation in the presence of large field inhomogeneities using a graph cut algorithm. *Magn Reson Med.* 2010;63(1):79-90. doi:10.1002/mrm.22177 [PMID: 19859956](https://pubmed.ncbi.nlm.nih.gov/19859956/)

**ISMRM Fat-Water Separation Toolbox**

Hu HH, Börnert P, Hernando D, Kellman P, Ma J, Reeder S, Sirlin C. ISMRM workshop on fat-water separation: Insights, applications and progress in MRI. *Magn Reson Med.* 2012;68(2):378-388. doi:10.1002/mrm.24369

Toolbox: <https://www.ismrm.org/workshops/FatWater12/data.htm>

**SPARCQ**

Rossi GMC, Mackowiak ALC, Açikgöz BC, Pierzchała K, Kober T, Hilbert T, Bastiaansen JAM. SPARCQ: A new approach for fat fraction mapping using asymmetries in the phase-cycled balanced SSFP signal profile. *Magn Reson Med.* 2023;90(6):2348-2361. doi:10.1002/mrm.29813

Please also cite the manuscript above when using this ΔB0-corrected implementation.

## License

Academic, non-commercial research only. The bundled graph-cut code is Diego Hernando / University of Wisconsin academic software; see `LICENSE`.
