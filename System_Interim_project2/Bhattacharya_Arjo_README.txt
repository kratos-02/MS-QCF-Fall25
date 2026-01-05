================================================================================
ISYE 6767 - Project 2: Data Analysis and Machine Learning
README
================================================================================

PYTHON VERSION AND DEPENDENCIES
-------------------------------
Python Version: 3.7.x

Required Packages:
- pandas >= 1.3.0
- numpy >= 1.21.0
- scikit-learn >= 0.24.0
- yfinance >= 0.1.70
- backtrader >= 1.9.76
- quantstats >= 0.0.59
- ta >= 0.10.0

Install dependencies:
    pip install pandas numpy scikit-learn yfinance backtrader quantstats ta

EXECUTION INSTRUCTIONS
----------------------
Run the main script:
    python Bhattacharya_Arjo_project2.py

OUTPUT FILES
------------
1. small_universe_results.csv - ML metrics for small universe
2. small_universe_backtest.csv - Backtest metrics for small universe
3. large_universe_results.csv - ML metrics for large universe
4. large_universe_top10_by_accuracy.csv - Top 10 stocks with backtest results

All outputs are saved to the 'final_outputs' directory.
Trading reports (HTML) are saved to the 'trading_reports' directory.

NOTES
-----
- HIBB and WBA tickers were excluded due to data availability issues
- Models are trained separately for each stock
- 60/40 train/test split is used (chronological)
- TimeSeriesSplit cross-validation respects temporal order
================================================================================
