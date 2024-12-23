function ilr = ILR(step, u, coeff_4_d_3_m_SQRDT, SQRT2D2, SQRT2M2minus2p5)

if step == 1
    ilr = 0;
elseif step == 2
    ilr = coeff_4_d_3_m_SQRDT * (u(0)*(1 - SQRT2D2) + u(1)*(SQRT2M2minus2p5));
else
    j = 1:step-2;
    ilr = coeff_4_d_3_m_SQRDT * (u(0)*((step - 1).^(1.5) - (step - 1.5)*sqrt(step))...
        + sum(u(j).*((step - j - 1).^(1.5) - 2*(step - j).^(1.5) + (step - j + 1).^(1.5)))...
        + u(step - 1)*(SQRT2M2minus2p5));
    % ilr = coeff_4_d_3_m_SQRDT * (u(0)*((step - 1.0D0).^(1.5) - (step - 1.5D0)*sqrt(step)) + u(step - 1)*(SQRT2M2minus2p5));
    % for j = 1:step-2
    %     ilr = ilr + coeff_4_d_3_m_SQRDT * (u(j).*((step - j - 1.0D0).^(1.5) - 2.0D0*(step - j).^(1.5) + (step - j + 1.0D0).^(1.5)));
    % end
end

% if step == 1
%     IL = 0;
% elseif step == 2
%     IL = coeff_4_d_3_m_SQRDT * (ul(0)*(1 - SQRT2D2) + ul(1)*(SQRT2M2minus2p5));
% else
%     j = 1:step-2;
%     IL = coeff_4_d_3_m_SQRDT * (ul(0)*((step - 1).^(1.5) - (step - 1.5)*sqrt(step))...
%         + sum(ul(j).*((step - j - 1).^(1.5) - 2*(step - j).^(1.5) + (step - j + 1).^(1.5)))...
%         + ul(step - 1)*(SQRT2M2minus2p5));
%     % IL = coeff_4_d_3_m_SQRDT * (ul(0)*((step - 1.0D0).^(1.5) - (step - 1.5D0)*sqrt(step)) + ul(step - 1)*(SQRT2M2minus2p5));
%     % for j = 1:step-2
%     %     IL = IL + coeff_4_d_3_m_SQRDT * (ul(j).*((step - j - 1.0D0).^(1.5) - 2.0D0*(step - j).^(1.5) + (step - j + 1.0D0).^(1.5)));
%     % end
% end
% 
% if step == 1
%     IR = 0;
% elseif step == 2
%     IR = coeff_4_d_3_m_SQRDT * (ur(0)*(1 - SQRT2D2) + ur(1)*(SQRT2M2minus2p5));
% else
%     j = 1:step-2;
%     IR = coeff_4_d_3_m_SQRDT * (ur(0)*((step - 1).^(1.5) - (step - 1.5)*sqrt(step))...
%         + sum(ur(j).*((step - j - 1).^(1.5) - 2*(step - j).^(1.5) + (step - j + 1).^(1.5)))...
%         + ur(step - 1)*(SQRT2M2minus2p5));
%     % IR = coeff_4_d_3_m_SQRDT * (ur(0)*((step - 1.0D0).^(1.5) - (step - 1.5D0)*sqrt(step)) + ur(step - 1)*(SQRT2M2minus2p5));
%     % for j = 1:step-2
%     %     IR = IR + coeff_4_d_3_m_SQRDT * (ur(j).*((step - j - 1.0D0).^(1.5) - 2.0D0*(step - j).^(1.5) + (step - j + 1.0D0).^(1.5)));
%     % end
% end
end