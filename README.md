# Two-Agent New Keynesian Model: Heterogeneous Response to Monetary Tightening

Companion code to:

**Gokhale, A.A. (2026). "Who Leaves When Rates Rise? Balance Sheet Heterogeneity and the Emigration Response to Monetary Tightening." SSRN Working Paper. [ssrn.com/abstract=7604979](https://ssrn.com/abstract=7604979)**

## What this does

Simulates a two-agent New Keynesian model with mortgaged homeowners and renters facing an identical interest rate shock. Because owners have direct mortgage exposure (θ_O = 0.6) and renters do not (θ_R = 0.0), the model predicts a differential consumption contraction — and a corresponding differential in emigration incentives.

This theoretical mechanism rationalises the empirical finding in the companion paper: a 2.14 percentage point increase in emigration intentions among mortgage-exposed households following the 2022–23 UK tightening cycle.

## Run it

```julia
julia two_agent_nk.jl
```

Generates two figures:
- `irf_clean.png` — consumption IRFs by household type
- `irf_emigration.png` — consumption contraction and emigration intent gap

## Key parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| θ_O | 0.6 | Mortgage exposure: owners |
| θ_R | 0.0 | Mortgage exposure: renters |
| β | 0.99 | Discount factor |
| ρ_r | 0.7 | Rate shock persistence |

## Author

Anuja A. Gokhale · NUS Merit Scholar · [ssrn.com/author=10973290](https://ssrn.com/author=10973290)
