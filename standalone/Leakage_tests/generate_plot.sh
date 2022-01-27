#!/bin/bash

python plots.py
convert -trim -colors 32 NRPyLeakage_semi_analytic_results_beta_only.png NRPyLeakage_semi_analytic_results_beta_only.png
convert -trim -colors 32 NRPyLeakage_semi_analytic_results_all.png NRPyLeakage_semi_analytic_results_all.png
pngcrush NRPyLeakage_semi_analytic_results_beta_only.png pngout.png && mv pngout.png NRPyLeakage_semi_analytic_results_beta_only.png
pngcrush NRPyLeakage_semi_analytic_results_all.png pngout.png && mv pngout.png NRPyLeakage_semi_analytic_results_all.png
