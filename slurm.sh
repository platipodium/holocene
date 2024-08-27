#!/bin/bash -x
# SPDX-FileContributor: Kai Wirtz <kai.wirtz@hereon.de>
# SPDX-License-Identifier: CC0-1.0
# SPDX-FileCopyrightText: 2023-2024 Helmholtz-Zentrum hereon GmbH
#SBATCH --ntasks=1
#SBATCH --output=gg-%j-%a.stdout
#SBATCH --error=gg-%j-%a.stderr
#SBATCH --time=0:59:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=grid_growth
#SBATCH --account=???

srun Rscript grid_growth.r $SLURM_ARRAY_TASK_ID
