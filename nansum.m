% SPDX-FileCopyrightText: 2021-2024 Helmholtz-Zentrum hereon GmbH
% SPDX-FileContributor: Carsten Lemmen  <carsten.lemmen@hereon.de>
% SPDX-License-Identifier: CC0-1.0 

function nsum = nansum(arg)
  arg(isnan(arg)) = 0;
  nsum = sum(arg);