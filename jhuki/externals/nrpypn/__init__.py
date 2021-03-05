# BSD 2-Clause License

# Copyright (c) 2018, Zachariah Etienne All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.

# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from sympy import N

from jhuki.externals.nrpypn.NRPyPN_shortcuts import eval__P_t__and__P_r


def compute_quasicircular_momenta(
    mass_ratio, coordinate_distance, chi_plus, chi_minus
):
    """Wrapper around :py:func:`~eval__P_t__and__P_r does`. It is guaranteed that a number is returned

        # TODO: Clean up this docstring

    :param mass_ratio: Mass ratio of the binary (>= 1).
    :type mass_ratio: float
    :param coordinate_distance: Initial coordinate separation.
    :type coordinate_distance: float

    :param chi_plus: Dimensionless spin of the black hole on the positive side of the x
                     (or z) axis along the three directions.
    :type chi_plus: tuple/list with three numbers

    :param chi_minus: Dimensionless spin of the black hole on the negative side of the x
                     (or z) axis along the three directions.
    :type chi_minus: tuple/list with three numbers


    # The numerical values for
    #  * P_t: the tangential momentum
    #  * P_r: the radial momentum
    """
    (nchi1x, nchi1y, nchi1z), (nchi2x, nchi2y, nchi2z) = chi_plus, chi_minus

    return tuple(
        map(
            N,
            eval__P_t__and__P_r(
                mass_ratio,
                coordinate_distance,
                nchi1x,
                nchi1y,
                nchi1z,
                nchi2x,
                nchi2y,
                nchi2z,
            ),
        )
    )
