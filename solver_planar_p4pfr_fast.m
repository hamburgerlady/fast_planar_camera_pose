function [R1sol,R2sol,R3sol,tsol,fsol,ksol] = solver_planar_p4pfr_fast(X,U)
% function [Rsol,tsol,fsol,ksol] = solver_planar_p4pfr_fast(X,U)
%
% Input:
% X: 2x4 3D-coordinates (z=0)
% U: 2x4 image coordinates
% Output:
% returns n solutions, 
% Rsol: contains cell 1xn of rotations
% tsol: 3xn translations
% fsol: 1xn focal lengths
% ksol: 1xn radial distortion 
%
% (c) 2018 Magnus Oskarsson




% 3D coordinates Xi = (xi,yi,0) i = 1,...,4
% image coordinates Ui = (ui,vi) i = 1,...,4
% 1/f = w
% radial dist k


% v = [p11 p12 p14 p21 p22 p24]'

%X = randn(2,4);
%U = randn(2,4);



r = sum(U.^2);


M = [-U(2,1)*X(1,1) -U(2,1)*X(2,1) -U(2,1) U(1,1)*X(1,1) U(1,1)*X(2,1) U(1,1);...
     -U(2,2)*X(1,2) -U(2,2)*X(2,2) -U(2,2) U(1,2)*X(1,2) U(1,2)*X(2,2) U(1,2);...
     -U(2,3)*X(1,3) -U(2,3)*X(2,3) -U(2,3) U(1,3)*X(1,3) U(1,3)*X(2,3) U(1,3);...
     -U(2,4)*X(1,4) -U(2,4)*X(2,4) -U(2,4) U(1,4)*X(1,4) U(1,4)*X(2,4) U(1,4)];
 
    
N = null(M);

C = [U(1,1)*X(1,1) U(1,1)*X(2,1) U(1,1);...
     U(1,2)*X(1,2) U(1,2)*X(2,2) U(1,2);...
     U(1,3)*X(1,3) U(1,3)*X(2,3) U(1,3)];
 


D = [X(1,1)*N(1,1)+X(2,1)*N(2,1)+N(3,1) r(1)*(X(1,1)*N(1,1)+X(2,1)*N(2,1)+N(3,1)) r(1)*(X(1,1)*N(1,2)+X(2,1)*N(2,2)+N(3,2)) X(1,1)*N(1,2)+X(2,1)*N(2,2)+N(3,2);...
    X(1,2)*N(1,1)+X(2,2)*N(2,1)+N(3,1) r(2)*(X(1,2)*N(1,1)+X(2,2)*N(2,1)+N(3,1)) r(2)*(X(1,2)*N(1,2)+X(2,2)*N(2,2)+N(3,2)) X(1,2)*N(1,2)+X(2,2)*N(2,2)+N(3,2);...
    X(1,3)*N(1,1)+X(2,3)*N(2,1)+N(3,1) r(3)*(X(1,3)*N(1,1)+X(2,3)*N(2,1)+N(3,1)) r(3)*(X(1,3)*N(1,2)+X(2,3)*N(2,2)+N(3,2)) X(1,3)*N(1,2)+X(2,3)*N(2,2)+N(3,2)];
    

CiD = C\D;

%nomy = u4*(d34 + b*d31 + x4*(d14 + b*d11) + y4*(d24 + b*d21)) - b*n31 - n32 - x4*(n12 + b*n11) - y4*(n22 + b*n21);
%denny = r4*(n32 + b*n31 + x4*(n12 + b*n11) + y4*(n22 + b*n21)) - u4*(d33 + b*d32 + x4*(d13 + b*d12) + y4*(d23 + b*d22));

d11 = CiD(1,1);
d12 = CiD(1,2);
d13 = CiD(1,3);
d14 = CiD(1,4);
d21 = CiD(2,1);
d22 = CiD(2,2);
d23 = CiD(2,3);
d24 = CiD(2,4);
d31 = CiD(3,1);
d32 = CiD(3,2);
d33 = CiD(3,3);
d34 = CiD(3,4);
n11 = N(1,1);
n12 = N(1,2);
n21 = N(2,1);
n22 = N(2,2);
n31 = N(3,1);
n32 = N(3,2);
n41 = N(4,1);
n42 = N(4,2);
n51 = N(5,1);
n52 = N(5,2);
%n61 = N(6,1);
%n62 = N(6,2);

u4 = U(1,4);
r4 = r(4);
x4 = X(1,4);
y4 = X(2,4);

knomy_b = n31 - d31*u4 + n11*x4 + n21*y4 - d11*u4*x4 - d21*u4*y4;
knomy_1 = n32 - d34*u4 + n12*x4 + n22*y4 - d14*u4*x4 - d24*u4*y4;
kdenny_b = d32*u4 - n31*r4 + d12*u4*x4 + d22*u4*y4 - n11*r4*x4 - n21*r4*y4;
kdenny_1 = d33*u4 - n32*r4 + d13*u4*x4 + d23*u4*y4 - n12*r4*x4 - n22*r4*y4;

c11_0 = n12*n22 + n42*n52;
c11_1 = n11*n22 + n12*n21 + n41*n52 + n42*n51;
c11_2 = n11*n21 + n41*n51;

c21_0 = n12^2 - n22^2 + n42^2 - n52^2;
c21_1 = 2*n11*n12 - 2*n21*n22 + 2*n41*n42 - 2*n51*n52;
c21_2 = n11^2 - n21^2 + n41^2 - n51^2;

c12_0 = (d14*kdenny_1 + d13*knomy_1)*(d24*kdenny_1 + d23*knomy_1);
c12_1 = (d24*kdenny_1 + d23*knomy_1)*(d11*kdenny_1 + d14*kdenny_b + d12*knomy_1 + d13*knomy_b) + (d14*kdenny_1 + d13*knomy_1)*(d21*kdenny_1 + d24*kdenny_b + d22*knomy_1 + d23*knomy_b);
c12_2 = (d11*kdenny_1 + d14*kdenny_b + d12*knomy_1 + d13*knomy_b)*(d21*kdenny_1 + d24*kdenny_b + d22*knomy_1 + d23*knomy_b) + (d14*kdenny_1 + d13*knomy_1)*(d21*kdenny_b + d22*knomy_b) + (d24*kdenny_1 + d23*knomy_1)*(d11*kdenny_b + d12*knomy_b);
c12_3 = (d21*kdenny_b + d22*knomy_b)*(d11*kdenny_1 + d14*kdenny_b + d12*knomy_1 + d13*knomy_b) + (d11*kdenny_b + d12*knomy_b)*(d21*kdenny_1 + d24*kdenny_b + d22*knomy_1 + d23*knomy_b);
c12_4 = (d11*kdenny_b + d12*knomy_b)*(d21*kdenny_b + d22*knomy_b);

c22_0 = (d14*kdenny_1 + d24*kdenny_1 + d13*knomy_1 + d23*knomy_1)*(d14*kdenny_1 - d24*kdenny_1 + d13*knomy_1 - d23*knomy_1);
c22_1 = (d14*kdenny_1 - d24*kdenny_1 + d13*knomy_1 - d23*knomy_1)*(d11*kdenny_1 + d21*kdenny_1 + d14*kdenny_b + d24*kdenny_b + d12*knomy_1 + d22*knomy_1 + d13*knomy_b + d23*knomy_b) + (d14*kdenny_1 + d24*kdenny_1 + d13*knomy_1 + d23*knomy_1)*(d11*kdenny_1 - d21*kdenny_1 + d14*kdenny_b - d24*kdenny_b + d12*knomy_1 - d22*knomy_1 + d13*knomy_b - d23*knomy_b);
c22_2 = (d11*kdenny_1 + d21*kdenny_1 + d14*kdenny_b + d24*kdenny_b + d12*knomy_1 + d22*knomy_1 + d13*knomy_b + d23*knomy_b)*(d11*kdenny_1 - d21*kdenny_1 + d14*kdenny_b - d24*kdenny_b + d12*knomy_1 - d22*knomy_1 + d13*knomy_b - d23*knomy_b) + (d14*kdenny_1 + d24*kdenny_1 + d13*knomy_1 + d23*knomy_1)*(d11*kdenny_b - d21*kdenny_b + d12*knomy_b - d22*knomy_b) + (d14*kdenny_1 - d24*kdenny_1 + d13*knomy_1 - d23*knomy_1)*(d11*kdenny_b + d21*kdenny_b + d12*knomy_b + d22*knomy_b);
c22_3 = (d11*kdenny_b - d21*kdenny_b + d12*knomy_b - d22*knomy_b)*(d11*kdenny_1 + d21*kdenny_1 + d14*kdenny_b + d24*kdenny_b + d12*knomy_1 + d22*knomy_1 + d13*knomy_b + d23*knomy_b) + (d11*kdenny_b + d21*kdenny_b + d12*knomy_b + d22*knomy_b)*(d11*kdenny_1 - d21*kdenny_1 + d14*kdenny_b - d24*kdenny_b + d12*knomy_1 - d22*knomy_1 + d13*knomy_b - d23*knomy_b);
c22_4 = (d11*kdenny_b + d21*kdenny_b + d12*knomy_b + d22*knomy_b)*(d11*kdenny_b - d21*kdenny_b + d12*knomy_b - d22*knomy_b);


poly = zeros(1,7);
poly(1) = c11_2*c22_4 - c12_4*c21_2;
poly(2) = c11_1*c22_4 + c11_2*c22_3 - c12_3*c21_2 - c12_4*c21_1;
poly(3) = c11_0*c22_4 + c11_1*c22_3 + c11_2*c22_2 - c12_2*c21_2 - c12_3*c21_1 - c12_4*c21_0;
poly(4) = c11_0*c22_3 + c11_1*c22_2 + c11_2*c22_1 - c12_1*c21_2 - c12_2*c21_1 - c12_3*c21_0;
poly(5) = c11_0*c22_2 + c11_1*c22_1 + c11_2*c22_0 - c12_0*c21_2 - c12_1*c21_1 - c12_2*c21_0;
poly(6) = c11_0*c22_1 + c11_1*c22_0 - c12_0*c21_1 - c12_1*c21_0;
poly(7) = c11_0*c22_0 - c12_0*c21_0;

bsol = roots(poly)';
lille = 1e-7;

bsol = bsol(abs(imag(bsol))<lille);

fsol2 = -(c11_0+c11_1*bsol+c11_2*bsol.^2)./((c12_0+c12_1*bsol+c12_2*bsol.^2+c12_3*bsol.^3+c12_4*bsol.^4)).*(kdenny_b*bsol+kdenny_1).^2;
okids = fsol2>0;


bsol = bsol(okids);
fsol = sqrt(fsol2(okids));
nn = sum(okids);


ksol = (knomy_b*bsol+knomy_1)./(kdenny_b*bsol+kdenny_1);

vsol = N*[bsol;ones(1,nn)];

v2sol = CiD*[bsol; ksol.*bsol; ksol; ones(1,nn)];

nr = repmat(sqrt(vsol(1,:).^2+vsol(4,:).^2+fsol2(okids).*v2sol(1,:).^2),3,1);
R1sol = [vsol(1,:);vsol(4,:);fsol.*v2sol(1,:)]./nr;
R2sol = [vsol(2,:);vsol(5,:);fsol.*v2sol(2,:)]./nr;

%sg = repmat(sign(X(:,4)'*[R1sol(3,:); R2sol(3,:)]+v2sol(3,:)./nr(1,:)),3,1);

sg = repmat(mode(sign(X(:,:)'*[R1sol(3,:); R2sol(3,:)]+repmat(v2sol(3,:)./nr(1,:),4,1))),3,1);


R1sol = R1sol.*sg;
R2sol = R2sol.*sg;


R3sol = [R1sol(2,:).*R2sol(3,:) - R1sol(3,:).*R2sol(2,:);...
         R1sol(3,:).*R2sol(1,:) - R1sol(1,:).*R2sol(3,:);...
         R1sol(1,:).*R2sol(2,:) - R1sol(2,:).*R2sol(1,:)];

tsol = sg.*[vsol(3,:);vsol(6,:);fsol.*v2sol(3,:)]./nr;




