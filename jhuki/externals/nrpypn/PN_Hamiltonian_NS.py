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
# As documented in the NRPyPN notebook
# PN-Hamiltonian-Nonspinning.ipynb, this Python script
# generates nonspinning pieces of the post-Newtonian (PN)
# Hamiltonian, up to and including third PN order.

# Basic functions:
# f_H_Newt__H_NS_1PN__H_NS_2PN(m1,m2, pU, nU, q): Compute H_Newt,
#                                                 H_NS_1PN, and H_NS_2PN
#                                                 and store to global
#                                                 variables of the same
#                                                 names.
# f_H_NS_3PN: Compute H_NS_3PN, and store to global variable of same name

# Author:  Zach Etienne
#          zachetie **at** gmail **dot* com

# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy+ depends
import jhuki.externals.nrpypn.indexedexpNRPyPN as ixp  # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
from jhuki.externals.nrpypn.NRPyPN_shortcuts import (
    div,
    dot,
)  # NRPyPN: shortcuts for e.g., vector operations


def f_H_Newt__H_NS_1PN__H_NS_2PN(m1, m2, PU, nU, q):
    mu = m1 * m2 / (m1 + m2)
    eta = m1 * m2 / (m1 + m2) ** 2
    pU = ixp.zerorank1()
    for i in range(3):
        pU[i] = PU[i] / mu

    global H_Newt, H_NS_1PN, H_NS_2PN
    H_Newt = mu * (+div(1, 2) * dot(pU, pU) - 1 / q)

    H_NS_1PN = mu * (
        +div(1, 8) * (3 * eta - 1) * dot(pU, pU) ** 2
        - div(1, 2) * ((3 + eta) * dot(pU, pU) + eta * dot(nU, pU) ** 2) / q
        + div(1, 2) / q ** 2
    )

    H_NS_2PN = mu * (
        +div(1, 16) * (1 - 5 * eta + 5 * eta ** 2) * dot(pU, pU) ** 3
        + div(1, 8)
        * (
            +(5 - 20 * eta - 3 * eta ** 2) * dot(pU, pU) ** 2
            - 2 * eta ** 2 * dot(nU, pU) ** 2 * dot(pU, pU)
            - 3 * eta ** 2 * dot(nU, pU) ** 4
        )
        / q
        + div(1, 2)
        * ((5 + 8 * eta) * dot(pU, pU) + 3 * eta * dot(nU, pU) ** 2)
        / q ** 2
        - div(1, 4) * (1 + 3 * eta) / q ** 3
    )


def f_H_NS_3PN(m1, m2, PU, nU, q):
    mu = m1 * m2 / (m1 + m2)
    eta = m1 * m2 / (m1 + m2) ** 2
    pU = ixp.zerorank1()
    for i in range(3):
        pU[i] = PU[i] / mu

    global H_NS_3PN
    H_NS_3PN = mu * (
        +div(1, 128)
        * (-5 + 35 * eta - 70 * eta ** 2 + 35 * eta ** 3)
        * dot(pU, pU) ** 4
        + div(1, 16)
        * (
            +(-7 + 42 * eta - 53 * eta ** 2 - 5 * eta ** 3) * dot(pU, pU) ** 3
            + (2 - 3 * eta) * eta ** 2 * dot(nU, pU) ** 2 * dot(pU, pU) ** 2
            + 3 * (1 - eta) * eta ** 2 * dot(nU, pU) ** 4 * dot(pU, pU)
            - 5 * eta ** 3 * dot(nU, pU) ** 6
        )
        / q
        + (
            +div(1, 16) * (-27 + 136 * eta + 109 * eta ** 2) * dot(pU, pU) ** 2
            + div(1, 16)
            * (+17 + 30 * eta)
            * eta
            * dot(nU, pU) ** 2
            * dot(pU, pU)
            + div(1, 12) * (+5 + 43 * eta) * eta * dot(nU, pU) ** 4
        )
        / q ** 2
        + (
            +(
                -div(25, 8)
                + (div(1, 64) * sp.pi ** 2 - div(335, 48)) * eta
                - div(23, 8) * eta ** 2
            )
            * dot(pU, pU)
            + (-div(85, 16) - div(3, 64) * sp.pi ** 2 - div(7, 4) * eta)
            * eta
            * dot(nU, pU) ** 2
        )
        / q ** 3
        + (+div(1, 8) + (div(109, 12) - div(21, 32) * sp.pi ** 2) * eta)
        / q ** 4
    )
