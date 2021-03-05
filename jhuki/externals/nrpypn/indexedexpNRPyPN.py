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
# indexedexpNRPyPN.py: functions related to indexed expressions,
# including e.g., tensors and pseudotensors.
# *** This is a stripped-down version of the indexedexp.py module
#     in the NRPy+ root directory.

# Step 1: Load needed modules
import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy+ depends


def zerorank1(DIM=-1):
    if DIM == -1:
        DIM = 3  # default to 3D
    return [sp.sympify(0) for i in range(DIM)]


def zerorank3(DIM=-1):
    if DIM == -1:
        DIM = 3  # default to 3D
    return [
        [[sp.sympify(0) for i in range(DIM)] for j in range(DIM)]
        for k in range(DIM)
    ]


def declarerank1(objname, DIM=-1):
    if DIM == -1:
        DIM = 3  # default to 3D
    return [sp.sympify(objname + str(i)) for i in range(DIM)]


class NonInvertibleMatrixError(ZeroDivisionError):
    """ Matrix Not Invertible; Division By Zero """


# Define the rank-3 version of the Levi-Civita symbol.
def LeviCivitaSymbol_dim3_rank3():
    LeviCivitaSymbol = zerorank3(DIM=3)

    for i in range(3):
        for j in range(3):
            for k in range(3):
                # From https://codegolf.stackexchange.com/questions/160359/levi-civita-symbol :
                LeviCivitaSymbol[i][j][k] = (
                    (i - j) * (j - k) * (k - i) * sp.Rational(1, 2)
                )
    return LeviCivitaSymbol
