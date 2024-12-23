function d = DAR(cu, cu_p, WL, WR, DPART, Nz, C1, C2, coeff_4_d_3_m_SQRDT_m_SQRTKAPPA, coeff_1i_m_SQRDZ_m_kappa) %#codegen

Nzm1 = Nz - 1;

d1  = DD( 1, WL, DPART(1),  C1, coeff_4_d_3_m_SQRDT_m_SQRTKAPPA);
d2_Nzm1 = -coeff_1i_m_SQRDZ_m_kappa * (cu_p  + cu) + DPART(2:Nzm1);
dNz = DD(-1, WR, DPART(Nz), C2, coeff_4_d_3_m_SQRDT_m_SQRTKAPPA);

d = [d1; d2_Nzm1; dNz];

% D(1)  =   coeff_C1_m_coeff_4_d_3_m_SQRDT * WL(IDX(step)) + D_0_PART;
% D(Nz) = - coeff_C2_m_coeff_4_d_3_m_SQRDT * WR(IDX(step)) + D_END_PART;

end

function dd = DD(sgn, W, DPART, C, coeff_4_d_3_m_SQRDT_m_SQRTKAPPA)
dd  = sgn * C * coeff_4_d_3_m_SQRDT_m_SQRTKAPPA * W + DPART;
end