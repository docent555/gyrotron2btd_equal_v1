function [OUTF, OUTJ, p, Eff, Omega, ConLow, jout] = gyroscr(Nz, Nzi, Nt, Ne, ZAxis, ZAxisi, TAxis, Delta, Ic, dt, dz, dzi, tol, kpar2, INTT, INTZ, OUTNz, OUTNt, InitialField) %#codegen

WR = complex(zeros(Nt,1));
FNz = complex(zeros(Nt,1));
FNzm1 = complex(zeros(Nt,1));
BNz = complex(zeros(Nt,1));
BNzm1 = complex(zeros(Nt,1));
SigmaNz = complex(zeros(Nt,1));
SigmaNzm1 = complex(zeros(Nt,1));
WL = complex(zeros(Nt,1));
F0 = complex(zeros(Nt,1));
F1 = complex(zeros(Nt,1));
B0 = complex(zeros(Nt,1));
B1 = complex(zeros(Nt,1));
Sigma0 = complex(zeros(Nt,1));
Sigma1 = complex(zeros(Nt,1));
fmax = zeros(Nt, 1);
jmax = zeros(Nt, 1);
field = complex(zeros(Nz,1));
testfield = complex(zeros(Nz,1));
field_p = complex(zeros(Nz,1));
rfield_p = complex(zeros(Nz,1));
lfield_p = complex(zeros(Nz,1));
B_p = complex(zeros(Nz,1));
B = complex(zeros(Nz,1));
p_p = complex(zeros(Nzi,Ne));
Jr = complex(zeros(Nzi,1));
J_p = complex(zeros(Nz,1));
J = complex(zeros(Nz,1));
OUTF = complex(zeros(OUTNz, OUTNt));
OUTJ = complex(zeros(OUTNz, OUTNt));
Eff = zeros(Nt,1);
Omega = zeros(Nt,1);
ConLow = zeros(Nt,1);
ZAxis_ipart = ZAxis(ZAxis <= ZAxisi(end));
ZAxis_ipart_end = length(ZAxis_ipart);
% theta = zeros(Nz, Ne);
% p = zeros(Nz, Ne);
% pv = zeros(Nz, 2*Ne);
% p0 = zeros(Ne,1);
% p0v = zeros(2*Ne,1);
% reidx = zeros(1,Ne);
% imidx = zeros(1,Ne);

if INTZ > 1
    IZ = 0:INTZ:length(ZAxis);
    IZ(1) = 1;
    SIZEZ = length(IZ);
else
    IZ = 1:INTZ:length(ZAxis);
end

% kpar2 = zeros(length(ZAxis),1);
% N = length(ZAxis);
AA = complex(zeros(Nz,1));
BB = complex(zeros(Nz-1,1));
CC = complex(zeros(Nz-1,1));
DD = complex(zeros(Nz,1));

% SQR2 = sqrt(2.0D0);
SQRT2M2 = 2.828427124746190;
SQRT2D2 = 0.707106781186548;
SQRTDT = sqrt(dt);
SQRDZ = dz*dz;
SQRT2M2minus2p5 = SQRT2M2 - 2.5;

Nzm1 = Nz - 1;
C0 = 1.0D0;
CL = -1i*kpar2(1);
CR = -1i*kpar2(Nz);
% C1 = 1.0D0/sqrt(1i*pi);
C1 = 0;
C2 = 1.0D0/sqrt(1i*pi);
W0 = ((-1i*2.0D0/3.0D0*C0*dz/dt) - 1.0D0/dz);
W1 = ((-1i*C0/3.0D0*dz/dt) + 1.0D0/dz);
WNz = -((-1i*2.0D0/3.0D0*C0*dz/dt) - 1.0D0/dz);
WNzm1 = -((-1i*C0/3.0D0*dz/dt) + 1.0D0/dz);

AA(1) =  1.0D0 - 4.0D0/3.0D0*C1*W0*SQRTDT;
AA(2:Nzm1) = -2.0D0*(1.0D0 + 1i * SQRDZ/dt*C0);
AA(Nz) = 1.0D0 + 4.0D0/3.0D0*C2*WNz*SQRTDT;
BB(1) =   -4.0D0/3.0D0*C1*W1*SQRTDT;
BB(2:Nzm1) = 1.0D0;
CC(1:Nzm1 - 1) = 1.0D0;
CC(Nzm1) = 4.0D0/3.0D0*C2*WNzm1*SQRTDT;

M = spdiags([[CC; 0] AA [0 ;BB]], -1:1, Nz, Nz);

% Initial values
jout = 1;
field(:,1) = InitialField;
% FforP = SF(ZAxisi(1:Nzi-1));
OUTF(:, jout) = field(IZ,1);
th0 = 2.0D0*pi*(0:Ne-1)/Ne;
p0 = exp(1i*th0).';
p0v = [real(p0); imag(p0)];
reidx = 1:Ne;
imidx = Ne+1:2*Ne;

% PforF   = complex(zeros(length(ZAxis_ipart), Ne)); % то что попадает на точки
% PforF_p = complex(zeros(length(ZAxis_ipart), Ne)); % уравнения возбуждения

optsp = odeset('RelTol',1e-8,'AbsTol',1e-10);
% tmp = cell(Ne,1);
% Sp = coder.nullcopy(tmp);
% Sp = cell(Ne,1);
SF = griddedInterpolant(ZAxis, field,'spline');
p = oscill_ode(SF, Nzi, ZAxisi, Delta, p0v, reidx, imidx, optsp);
% p = oscill_reim(field, Nzi, ZAxis, ZAxisi, Delta, p0v, reidx, imidx, opts);
% p = oscill_cmplx(field(1:Nzi), ZAxis(1:Nzi), Delta, p0);
% for i=1:Ne
%     Sp{i} = griddedInterpolant(ZAxisi, p(:,i),'spline');
%     % PforF(:,i) = Sp{i}(ZAxis_ipart);
%     % plot(ZAxisi,abs(Sp{i}(ZAxisi)))
%     % hold on
% end
% for i=1:Ne
%     PforF(:,i) = Sp{i}(ZAxis_ipart);   
% end
% J(1:ZAxis_ipart_end,1) = Ic * sum(PforF, 2)/Ne;
Jr(:) = Ic * sum(p, 2)/Ne;
SJ = griddedInterpolant(ZAxisi, Jr,'spline');
J(1:ZAxis_ipart_end,1) = SJ(ZAxis_ipart);
B(:,1) = J(:) - 1i*kpar2(:).*field(:);
OUTJ(:,jout) = J(IZ,1);

% IDX = @(j) (j + 1);

fmax(IDX(0)) = max(abs(field(:,1)));
jmax(IDX(0)) = max(abs(B(:,1)));
F0(IDX(0)) = field(1);
F1(IDX(0)) = field(2);
B0(IDX(0)) = B(1);
B1(IDX(0)) = B(2);
Sigma0(IDX(0)) = 0;
Sigma1(IDX(0)) = 0;
FNz(IDX(0)) = field(Nz);
FNzm1(IDX(0)) = field(Nzm1);
BNz(IDX(0)) = B(Nz);
BNzm1(IDX(0)) = B(Nzm1);
SigmaNz(IDX(0)) = 0;
SigmaNzm1(IDX(0)) = 0;
Eff(IDX(0)) = 1 - sum(abs(p(Nzi,:)).^2)/Ne;
Omega(IDX(0)) = 0;

WL(IDX(0)) = -dz * (-1i/6 * (2 * B0(IDX(0)) + B1(IDX(0))));
WR(IDX(0)) =  dz * (-1i/6 * (2 * BNz(IDX(0)) + BNzm1(IDX(0))));

SHOW = 1;
if SHOW == 1
    [lhfmax, lhfabs, lhjmax, lhjabs, hFig] = makeFig(ZAxis, ZAxis_ipart , TAxis);
end

%Coefficients
coeff_1i_m_C0_m_2_d_3_d_dt = 1i*C0*2.0D0/3.0D0/dt;
coeff_1i_m_C0_d_3_d_dt = 1i*C0/3.0D0/dt;
coeff_1i_d_6 = 1i/6.0D0;
coeff_4_d_3_m_SQRDT = 4.0D0/3.0D0*SQRTDT;
coeff_2_d_3_m_SQRDT = 2.0D0/3.0D0*SQRTDT;
coeff_1i_m_SQRDZ = 1i*SQRDZ;
coeff_1i_m_C0_m_SQRDZ_d_dt = 1i*C0*SQRDZ/dt;
% coeff_C1_m_coeff_4_d_3_m_SQRDT = C1*coeff_4_d_3_m_SQRDT;
% coeff_C2_m_coeff_4_d_3_m_SQRDT = C2*coeff_4_d_3_m_SQRDT;
coeff_exp_CL_m_dt = exp(CL*dt);
coeff_CL_m_dt = CL*dt;
coeff_exp_CR_m_dt = exp(CR*dt);
coeff_CR_m_dt = CR*dt;
coeff_dz_m_coeff_1i_d_6 = dz*coeff_1i_d_6;

wp  = @(sgn, fEND, fENDm1, SigmaEND, SigmaENDm1) W_PART(sgn, fEND, fENDm1, ...
    SigmaEND, SigmaENDm1, dz, coeff_1i_m_C0_m_2_d_3_d_dt, coeff_1i_m_C0_d_3_d_dt);
w   = @(sgn, cuEND, cuENDm1, cu_pEND, cu_pENDm1, W_PART) WLR(sgn, cuEND, cuENDm1, cu_pEND, cu_pENDm1, W_PART, coeff_dz_m_coeff_1i_d_6);
ilr = @(step, u) ILR(step, u, coeff_4_d_3_m_SQRDT, SQRT2D2, SQRT2M2minus2p5);
dp  = @(IL, IR, field, WL, WR) DAR_PART(IL, IR, field, WL, WR, W0, W1, WNzm1, WNz, ...
    coeff_2_d_3_m_SQRDT, coeff_1i_m_C0_m_SQRDZ_d_dt, ...
    coeff_exp_CL_m_dt, coeff_exp_CR_m_dt, C1, C2, Nz);
d   = @(B, B_p, WL, WR, DPART) DAR(B, B_p, WL, WR, DPART, Nz, C1, C2, coeff_4_d_3_m_SQRDT, coeff_1i_m_SQRDZ);
sgm = @(f1, f2, j1, j2, Sold) Sigma(f1, f2, j1, j2, Sold, coeff_1i_m_C0_d_3_d_dt, coeff_1i_d_6);

num_st_test_iter = 0;
fmax_glob_old = max(abs(field(:,1)));

fprintf('\n');
timerVal = tic;
for step=1:Nt-1

    if SHOW == 1
        lhfmax.YData(1:step) = Omega(1:step);
        lhfmax.XData(1:step) = TAxis(1:step);
        lhfabs.YData = abs(field);

        lhjmax.YData(1:step) = Eff(1:step);
        lhjmax.XData(1:step) = TAxis(1:step);
        lhjabs.YData = abs(J);

        drawnow
    end   

    WL_PART = wp(-1, field(1),  field(2),    Sigma0(IDX(step-1)),  Sigma1(IDX(step-1)));
    WR_PART = wp( 1, field(Nz), field(Nzm1), SigmaNz(IDX(step-1)), SigmaNzm1(IDX(step-1)));

    % WL_PART = -dz * ((-coeff_1i_m_C0_m_2_d_3_d_dt) * field(1)...
    %     + (-coeff_1i_m_C0_d_3_d_dt) * field(2)...
    %     - (2.0D0 * Sigma0(IDX(step-1)) + Sigma1(IDX(step-1))));    
    % WR_PART =  dz * ((-coeff_1i_m_C0_m_2_d_3_d_dt) * field(Nz)...
    %     + (-coeff_1i_m_C0_d_3_d_dt) * field(Nzm1)...
    %     - (2.0D0 * SigmaNz(IDX(step-1)) + SigmaNzm1(IDX(step-1))));

    WL(IDX(step)) = w( 1, B(1),  B(2),    B(1),  B(2),    WL_PART);    
    WR(IDX(step)) = w(-1, B(Nz), B(Nzm1), B(Nz), B(Nzm1), WR_PART);

    % WL(IDX(step)) = coeff_dz_m_coeff_1i_d_6*(4.0D0 * B(1) + 2.0D0 * B(2)) + WL_PART;
    % WR(IDX(step)) = - coeff_dz_m_coeff_1i_d_6*(4.0D0 * B(Nz) + 2.0D0 * B(Nzm1)) + WR_PART;

    %     u = @(j) (WNzm1 * FNzm1(IDX(j)) + WNz * FNz(IDX(j)) + WR(IDX(j))).' .* exp(coeff_CR_m_dt * (step - j));

    IL = ilr(step, @ul);
    IR = ilr(step, @ur);

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

    % DD(1) = 0;
    %         DD(1) = IN.TimeAmp * exp(1i * IN.TimeFreq * AxisTau(step));

 
    D_PART = dp(IL, IR, field, WL(IDX(step-1)), WR(IDX(step-1)));
    % D_0_PART = D_PART(1);
    % D_MIDDLE_PART = D_PART(2:Nzm1);
    % D_END_PART = D_PART(Nz);

    % D_0_PART =   C1 * (IL + coeff_2_d_3_m_SQRDT * (W0 * field(1)...
    %     + W1 * field(2) + WL(IDX(step-1))) * coeff_exp_CL_m_dt);  
    % D_MIDDLE_PART = 2.0D0 * (1.0D0 - coeff_1i_m_C0_m_SQRDZ_d_dt) .* field(2:Nzm1)...
    %     - (field(1:Nz - 2) + field(3:Nz));    
    % D_END_PART = - C2 * (IR + coeff_2_d_3_m_SQRDT * (WNzm1 * field(Nzm1)...
    %     + WNz * field(Nz) + WR(IDX(step-1))) * coeff_exp_CR_m_dt);

    DD = d(B(2:Nzm1), B(2:Nzm1), WL(IDX(step)), WR(IDX(step)), D_PART);    

    % DD(1)  =   coeff_C1_m_coeff_4_d_3_m_SQRDT * WL(IDX(step)) + D_0_PART;
    % DD(2:Nzm1) = -coeff_1i_m_SQRDZ * (2.0D0*B(2:Nzm1)) + D_MIDDLE_PART;
    % DD(Nz) = - coeff_C2_m_coeff_4_d_3_m_SQRDT * WR(IDX(step)) + D_END_PART; 

    % nesamosoglasovannoe pole
    field_p = M \ DD;
    % rfield_p = rtridag(CC,AA,BB,DD);
    % lfield_p = ltridag(CC,AA,BB,DD);
    % field_p = (rfield_p + lfield_p)/2.0D0;

    num_insteps = 0;
    maxfield = max(abs(field_p(:,1)));
    testfield = field_p;
    while 1
        num_insteps = num_insteps + 1;
        SF = griddedInterpolant(ZAxis, field_p,'spline');
        p_p = oscill_ode(SF, Nzi, ZAxisi, Delta, p0v, reidx, imidx, optsp);
        % p = oscill_reim(field, Nzi, ZAxis, ZAxisi, Delta, p0v, reidx, imidx, opts);
        % p = oscill_cmplx(field(1:Nzi), ZAxis(1:Nzi), Delta, p0);
        % for i=1:Ne
        %     Sp{i} = griddedInterpolant(ZAxisi, p_p(:,i),'spline');
        %     PforF_p(:,i) = Sp{i}(ZAxis_ipart);
        %     % plot(ZAxisi,abs(Sp{i}(ZAxisi)))
        %     % hold on
        % end
        % J_p(1:ZAxis_ipart_end,1) = Ic * sum(PforF, 2)/Ne;
        Jr(:) = Ic * sum(p_p, 2)/Ne;
        SJ = griddedInterpolant(ZAxisi, Jr,'spline');
        J_p(1:ZAxis_ipart_end,1) = SJ(ZAxis_ipart);
        B_p(:,1) = J_p(:) - 1i*kpar2(:).*field_p(:);

        WL(IDX(step)) = w( 1, B(1),  B(2),    B_p(1),  B_p(2),    WL_PART);    
        WR(IDX(step)) = w(-1, B(Nz), B(Nzm1), B_p(Nz), B_p(Nzm1), WR_PART);

        % WL(IDX(step)) =  coeff_dz_m_coeff_1i_d_6 * (2.0D0 * B_p(1)  + 2.0D0 * B(1)  + B_p(2)    + B(2))    + WL_PART;
        % WR(IDX(step)) = -coeff_dz_m_coeff_1i_d_6 * (2.0D0 * B_p(Nz) + 2.0D0 * B(Nz) + B_p(Nzm1) + B(Nzm1)) + WR_PART;

        % DD(1) =   coeff_C1_m_coeff_4_d_3_m_SQRDT * WL(IDX(step)) + D_0_PART;
        % DD(2:Nzm1) = -coeff_1i_m_SQRDZ * (B_p(2:Nzm1) + B(2:Nzm1)) + D_MIDDLE_PART;
        % DD(Nz) = -coeff_C2_m_coeff_4_d_3_m_SQRDT * WR(IDX(step)) + D_END_PART;

        DD = d(B(2:Nzm1), B_p(2:Nzm1), WL(IDX(step)), WR(IDX(step)), D_PART);

        % samosoglasovannoe pole
        field_p(:,1) = M \ DD;
        % rfield_p(:,1) = rtridag(CC,AA,BB,DD);
        % lfield_p(:,1) = ltridag(CC,AA,BB,DD);
        % field_p = (rfield_p + lfield_p)/2.0D0;

        maxdiff = max(abs(testfield - field_p));
        err = maxdiff/maxfield;
        if err < tol
            break
        end
        testfield = field_p;
        if num_insteps > 1000
            fprintf('\nToo many inner steps!\n');
            pause;
        end
    end

    field(:,1) = field_p(:,1);
    SF = griddedInterpolant(ZAxis, field,'spline');
    p = oscill_ode(SF, Nzi, ZAxisi, Delta, p0v, reidx, imidx, optsp);
    % p = oscill_reim(field, Nzi, ZAxis, ZAxisi, Delta, p0v, reidx, imidx, opts);
    % p = oscill_cmplx(field(1:Nzi), ZAxis(1:Nzi), Delta, p0);
    % for i=1:Ne
    %     Sp{i} = griddedInterpolant(ZAxisi, p(:,i),'spline');
    %     PforF(:,i) = Sp{i}(ZAxis_ipart);
    %     % plot(ZAxisi,abs(Sp{i}(ZAxisi)))
    %     % hold on
    % end
    % J(1:ZAxis_ipart_end,1) = Ic * sum(PforF, 2)/Ne;
    Jr(:) = Ic * sum(p, 2)/Ne;
    SJ = griddedInterpolant(ZAxisi, Jr,'spline');
    J(1:ZAxis_ipart_end,1) = SJ(ZAxis_ipart);
    B(:,1) = J(:) - 1i*kpar2(:).*field(:);
    fmax(IDX(step)) = max(abs(field(:,1)));
    jmax(IDX(step)) = max(abs(B(:,1)));

    F0(IDX(step)) =  field(1);
    F1(IDX(step)) = field(2);
    FNz(IDX(step)) =  field(Nz);
    FNzm1(IDX(step)) = field(Nzm1);
    B0(IDX(step)) = B(1);
    B1(IDX(step)) = B(2);
    BNz(IDX(step)) = B(Nz);
    BNzm1(IDX(step)) = B(Nzm1);

    %     Omega(IDX(step)) = (angle(field(Nz)) - angle(FNz(IDX(step-1))))/dt;
    Omega(IDX(step)) = imag(log(field(Nz)/FNz(IDX(step-1))))/dt;
    Eff(IDX(step)) = 1 -  sum(abs(p(Nzi,:)).^2)/Ne;

    if (mod(num_st_test_iter,1000))
        fmax_glob_new = max(abs(field(:,1)));
        if abs(fmax_glob_new - fmax_glob_old)/fmax_glob_old < tol
            jout = jout + 1;
            OUTF(:, jout) = field(IZ,1);
            OUTJ(:, jout) = J(IZ,1);
            fprintf('Emergency exit!\n');
            return;
        end
        num_st_test_iter = num_st_test_iter + 1;
        fmax_glob_old = fmax_glob_new;
    end

    if mod(step,INTT) == 0
        jout = jout + 1;
        OUTF(:, jout) = field(IZ,1);
        OUTJ(:, jout) = J(IZ,1);
    end

    % Sigma0(IDX(step)) = -(-coeff_1i_m_C0_d_3_d_dt) * field(1) ...
    %     + (-coeff_1i_m_C0_d_3_d_dt) * F0(IDX(step - 1)) ...
    %     -coeff_1i_d_6*(B(1) + B0(IDX(step - 1))) - Sigma0(IDX(step - 1));
    % Sigma1(IDX(step)) = -(-coeff_1i_m_C0_d_3_d_dt) * field(2) ...
    %     + (-coeff_1i_m_C0_d_3_d_dt) * F1(IDX(step - 1)) ...
    %     -coeff_1i_d_6*(B(2) + B1(IDX(step - 1))) - Sigma1(IDX(step - 1));
    % 
    % SigmaNz(IDX(step)) = -(-coeff_1i_m_C0_d_3_d_dt) * field(Nz) ...
    %     + (-coeff_1i_m_C0_d_3_d_dt) * FNz(IDX(step - 1)) ...
    %     -coeff_1i_d_6*(B(Nz) + BNz(IDX(step - 1))) - SigmaNz(IDX(step - 1));
    % SigmaNzm1(IDX(step)) = -(-coeff_1i_m_C0_d_3_d_dt) * field(Nzm1) ...
    %     + (-coeff_1i_m_C0_d_3_d_dt) * FNzm1(IDX(step - 1)) ...
    %     -coeff_1i_d_6*(B(Nzm1) + BNzm1(IDX(step - 1))) - SigmaNzm1(IDX(step - 1));

    Sigma0(IDX(step))    = sgm(field(1),    F0(IDX(step - 1)),    B(1),     B0(IDX(step - 1)),    Sigma0(IDX(step - 1)));    
    Sigma1(IDX(step))    = sgm(field(2),    F1(IDX(step - 1)),    B(2),     B1(IDX(step - 1)),    Sigma1(IDX(step - 1)));    
    SigmaNz(IDX(step))   = sgm(field(Nz),   FNz(IDX(step - 1)),   B(Nz),    BNz(IDX(step - 1)),   SigmaNz(IDX(step - 1)));
    SigmaNzm1(IDX(step)) = sgm(field(Nzm1), FNzm1(IDX(step - 1)), B(Nzm1),  BNzm1(IDX(step - 1)), SigmaNzm1(IDX(step - 1)));    

    k = step + 1;

    ConLow(IDX(step)) = (2*imag(field(Nz)*conj(dfdzNz(step)) - field(1)*conj(dfdz0(step))) - Ic*Eff(IDX(step)))/(Ic*Eff(IDX(step)))*100;

    fprintf(['\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'...
        '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'...
        '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'...
        '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'...
        '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'...
        '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'...
        'Step = %8i   Time = %10.4f   Bmax = %+15.10e   Jmax = %+15.10e   W = %+15.10e   E = %+15.10e   CL = %+10.5f %%'],...
        int64(step), TAxis(k), fmax(k), max(abs(B(:,1))), Omega(IDX(step)), Eff(IDX(step)), ConLow(IDX(step)));
end

OUTJ(:,jout) = J(IZ,1);

fprintf("\n\n\n");

ExecutionTime = toc(timerVal);

hours = fix(ExecutionTime/3600);
minutes = fix((ExecutionTime - fix(hours*3600))/60);
seconds = ExecutionTime - hours*3600 - minutes*60;

fprintf("ExecitionTime = %8.4f [h]   %8.4f [m]   %8.4f [s]\n", hours, minutes, seconds);

fprintf(" \n\n");

    function  f = ur(j)
        coder.inline("always");
        f = (WNzm1 * FNzm1(IDX(j)) + WNz * FNz(IDX(j)) + WR(IDX(j))).' .* exp(coeff_CR_m_dt * (step - j));
    end

    function  f = ul(j)
        coder.inline("always");
        f = (W1 * F1(IDX(j)) + W0 * F0(IDX(j)) + WL(IDX(j))).' .* exp(coeff_CL_m_dt * (step - j));
    end

    function j = IDX(j)
        coder.inline("always");
        j = j + 1;
    end

    function  dfdz = dfdzNz(j)
        coder.inline("always");
        dfdz = (WNzm1 * FNzm1(IDX(j)) + WNz * FNz(IDX(j)) + WR(IDX(j)));
    end

    function  dfdz = dfdz0(j)
        coder.inline("always");
        dfdz = (W1 * F1(IDX(j)) + W0 * F0(IDX(j)) + WL(IDX(j)));
    end
end




