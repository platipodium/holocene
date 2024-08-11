# SPDX-FileCopyrightText: 2023-204 Helmholtz-Zentrum hereon GmbH
# SPDX-FileCopyrightText: 2024 2023-2024 Helmholtz-Zentrum hereon GmbH
# SPDX-FileContributor: Kai W. Wirtz <kai.wirtz@hereon.de>
#
# SPDX-License-Identifier: GPL-3.0-or-later

#
# moving average filter
#
movavg <- function(x,n) {
    stopifnot(is.numeric(x), is.numeric(n))
    if (length(n) != 1 || ceiling(n != floor(n)) || n <= 1)
        stop("Time window must be a single integer greater 1.")
    nx <- length(x)
    if (n >= nx)
        stop("Time window length 'n' must be greater then length of time series.")
    y <- numeric(nx)

    n2 <- ceiling((n + 0)/2)
        for (k in 1:nx)  y[k] <- mean(x[max(k-n2,1):min(k+n2,nx)])
    return(y)
}
