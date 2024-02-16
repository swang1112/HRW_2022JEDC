# Replication Files for: "Proxy SVAR identification of monetary policy shocks - Monte Carlo evidence and insights for the US", Journal of Economic Dynamics and Control, 139, 104457
by Herwartz, H., H. Rohloff, and S. Wang

### Code:

1. **Monte Carlo Simulation:**
    - For replication of the results, please use `Code/Simulation/Replica.R` to extract the files from the folder `Code/Simulation/Server`. These results were produced on a server with the following specifications: 128 CPUs, each with an Intel(R) Xeon(R) CPU E5-4660 v4 @ 2.20GHz. They can be reproduced by calling `Calc.R` in the corresponding working directory indicated by the sample size.

2. **Empirical Application:**
    - The main file is `Code/Main.R`, which contains most of our calculations such as specifications and diagnostics of reduced-form VAR, identifications of structural VAR, computation of IRFs, and bootstrap inference.
    - **Supplement files:**
        - `Code/Fa.R`: Performs factor-augmented VAR analysis.
        - `Code/Hist_decomp.R`: Performs historical decompositions.
        - `Code/barplot.R`: Plots the cumulative contributions of Volcker's disinflation as a bar plot.
    - **Please always run `Code/Main.R` first to run other supplement files.**

3. **Functions and Procedures:**
    - Functions and procedures are collected in the local package `Code/Functions` as well as in `aux_funs`.

### Data:

1. `data/USA_Tri.csv`: Contains variables in the VAR system.
2. Folder `data/Instruments`: Contains all employed monetary policy proxies and purged instruments as implied by the CV model (`RR_s.csv` and `SZ_s.csv`).
3. `data/Factors/fred-database_code/current.csv`: Contains informational variables from which latent factors are extracted.
