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
# PN-Hamiltonian-SSS.ipynb, this Python script
# generates spin-spin-spin coupling pieces of the
# post-Newtonian (PN) Hamiltonian, up to and
# including 3PN order.

# Core functions:
# f_H_SSS_3PN(m1,m2, n12U,n21U, S1U,S2U, p1U,p2U, q)
#       Compute the complete H_SSS_3PN term and store to
#                     global variable of the same name.

# Author:  Zach Etienne
#          zachetie **at** gmail **dot* com

# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
from . import (
    indexedexpNRPyPN as ixp,  # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
)
from .NRPyPN_shortcuts import (  # NRPyPN: shortcuts for e.g., vector operations
    cross,
    div,
    dot,
)


#################################
#################################
# Step 1: 3PN spin-spin-spin term, from Eq. 3.12 of
#        Levi and Steinhoff (2015):
#     https://arxiv.org/abs/1410.2601
def f_H_SSS_3PN(m1, m2, n12U, n21U, S1U, S2U, p1U, p2U, r12):
    def f_H_SSS_3PN_pt(m1, m2, nU, S1U, S2U, p1U, p2U, r):
        p2_minus_m2_over_4m1_p1 = ixp.zerorank1()
        for i in range(3):
            p2_minus_m2_over_4m1_p1[i] = p2U[i] - m2 / (4 * m1) * p1U[i]
        H_SSS_3PN_pt = (
            +div(3, 2)
            * (
                +dot(S1U, S1U) * dot(S2U, cross(nU, p1U))
                + dot(S1U, nU) * dot(S2U, cross(S1U, p1U))
                - 5 * dot(S1U, nU) ** 2 * dot(S2U, cross(nU, p1U))
                + dot(nU, cross(S1U, S2U))
                * (+dot(S1U, p1U) - 5 * dot(S1U, nU) * dot(p1U, nU))
            )
            - 3
            * m1
            / (2 * m2)
            * (
                +dot(S1U, S1U) * dot(S2U, cross(nU, p2U))
                + 2 * dot(S1U, nU) * dot(S2U, cross(S1U, p2U))
                - 5 * dot(S1U, nU) ** 2 * dot(S2U, cross(nU, p2U))
            )
            - dot(cross(S1U, nU), p2_minus_m2_over_4m1_p1)
            * (dot(S1U, S1U) - 5 * dot(S1U, nU) ** 2)
        ) / (m1 ** 2 * r ** 4)
        return H_SSS_3PN_pt

    global H_SSS_3PN
    H_SSS_3PN = +f_H_SSS_3PN_pt(
        m1, m2, n12U, S1U, S2U, p1U, p2U, r12
    ) + f_H_SSS_3PN_pt(m2, m1, n21U, S2U, S1U, p2U, p1U, r12)
