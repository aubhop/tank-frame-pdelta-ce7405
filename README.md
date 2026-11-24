# Tank Frame P-Delta Analyzer (CE 7405, Fall 2026)

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=aubhop/tank-frame-pdelta&file=CE7405_FinalProj2.mlapp)

This MATLAB app analyzes a 2D braced steel frame with fixed bases that holds up an elevated cylindrical tank (R = 10 m). It runs both
**first-order** and **second-order (P-Delta, Newton-Raphson)** direct-stiffness analysis. The tank gravity load,
the wave drag F_D, and the uplift F_z come from 2D OpenFOAM storm simulations. Morison loads act on the legs.

## Run it in your browser (no install)
1. Click the **Open in MATLAB Online** button above.
2. Sign in with a MathWorks account. 
3. When asked, click **Clone** to copy the repository into your MATLAB Drive.
4. `CE7405_FinalProj2.mlapp` opens in App Designer. Click the green **Run** button.
   You can also type `CE7405_FinalProj2` in the Command Window.

## Try this example
The app opens with these preset values. Just press **RUN**:

| Input | Preset |
|---|---|
| Number of columns | 3 |
| Storm | Extreme (CFD sim #24) |
| Tank contents | Water |
| Fill level | 100 % |
| Analysis | 1st order + 2nd order |

Next, change **Number of columns to 2** and press RUN again. The P-Delta amplification factor goes up to about 1.08.

## Run it on your own computer instead
You need MATLAB R2021a or newer. No toolboxes are needed.
```
git clone https://github.com/aubhop/tank-frame-pdelta.git
```
Open the folder in MATLAB and run `CE7405_FinalProj2`.

## Files
| File | Purpose |
|---|---|
| `CE7405_FinalProj2.mlapp` | App Designer GUI |
| `runAnalysis.m` | Called by the RUN button. Reads the inputs, runs the solvers, fills the plots and the table. |
| `buildFrame.m` | Builds the nodes, elements, sections, and DOF map |
| `buildLoads.m` | Builds the tank load, the CFD storm loads (F_D, F_z), and the Morison loads |
| `assembleK.m` | Assembles the elastic and geometric stiffness matrices |
| `solveLinear.m` | First-order solution |
| `solvePDelta.m` | Second-order incremental Newton-Raphson solution and critical load factor |
| `plotModel.m`, `plotLoads.m`, `plotDeflected.m` | Plotting |

## Assumptions
- Fixed bases, 2D frame, linear-elastic steel.
- Member self-weight and tank shell weight are neglected.
- Section properties are representative values.

