# Delta Hedging Midterm — README

**Author:** Arjo Bhattacharya  
**Course:** ISYE6767
**Submission date:** 10/15/2025

---

## Project summary

This project implements and evaluates a **daily delta-hedging** strategy for European call options under the Black–Scholes model. It has two parts:

- **Task 1 (simulation)** — Monte Carlo simulation of GBM stock paths and evaluation of discrete daily delta hedging across simulated paths.
- **Task 2 (market-data)** — Apply the same hedging algorithm to real market data (GOOG) using option quotes, underlying prices, and short rates; compute implied vol by inverting Black–Scholes, calculate delta, perform daily rebalancing, and record hedging errors and P&L.

The code is written in C++ (numerical core) and Python (visualization & analysis). Unit tests validate key numerical routines (implied-vol inversion and delta).

---

## Repository contents

