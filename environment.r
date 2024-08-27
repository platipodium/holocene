#! /bin/env R
# SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de
# SPDX-License-Identifier: CC0-1.0
# SPDX-FileCopyrightText: 2024 Helmholtz-Zentrum hereon GmbH
#
# Script to setup your environment, i.e. installing necessary R libraries
# Run this with, e.g., the command `Rscript environment.r`

if (!require("pacman")) install.packages("pacman")
pacman::p_load(sp, rworldmap, stats, "R.matlab", "RColorBrewer", sf, rcarbon)
