# FX-Multi-FGARCH
This repository provide codes and dataset to replicate the procedure and results in the paper "Forecasting intraday foreign exchange volatility with functional GARCH approaches", by Fearghal Kearney, HanLin Shang, Yuqian Zhao

For relevant dataset, one could find from the shared folder:
https://www.dropbox.com/scl/fo/qfv0ufla9bnzeb5k0hefk/ACX6jtiDlSGdR2d6ymUoP1k?rlkey=dxyun3exl5vh1fsbvhoolg19q&st=5l1f7sz4&dl=0

There are 8 files in the package:

(1) Pre_tests.R : preliminary tests for intraday data, which one could get the results in Table 3.1 and Figure 4.1

(2) FPCA_basis.R : four types of functional principal component analysis (FPCA), including truncated FPCA, dynamic FPCA, long-range dependent FPCA, and multi-level FPCA. One could use this code to get the results in Table 5.1

(3) FGARCH-FGARCH-X.R : the main estimation functions of FGARCH(1,1) model and FGARCH-X model

(4) Forecasting_fgarch(x).R : functions for forecasting h-steps of intraday conditional volatility using FGARCH(1,1) or FGARCH-X model, the case of USD/EUR is taken as an example

(5) Forecasting_classical_vol.R : functions for forecasting h-steps of interdaily conditional volatility using 8 scalar models including: GARCH(1,1), GARCH-X, GJR-GARCH, FIGARCH, HAR, HAR-X, HAR-ARFIMA, RGARCH, and GARCH-MIDAS

(6) Forecasting_implementation.R : code to deply the forecasting and compare the forecasting errors. One could use the code from files (2)-(5) to get the results in Figures 5.1-5.3, Tables 5.2-5.6 

(7) VaR_backtesting.R : functions for backtesting the intraday VaR curves. One could use this code to get results in Table 6.1

(8) Application_VaRcorrect_strategy.R : functions for implementing VaR-corrected trading strategy. One could use the code from files (7)-(8) to get the results in Table 6.2 Figure 6.1
