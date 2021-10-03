function [out,varargout] = residuals_induction(x,rotor,tsr,idx,p,p1,p2,t_max,Rnew,c)
% Calculate the residuals for axial and tangential induction factor

% unpack inputs
R = Rnew; B = rotor.B;
% spanwise position
r = rotor.r(idx);

% unpack x
%c  = x(1);
a  = x(1);
ap = x(2);

% calculate that, cl and cl/cd from design polynomials
t    = thickness(r,p,t_max,rotor.R);
that = t ./ c;
cl   = x_des(that*100,p1(1,:),p2(1,:));
alpha = x_des(that*100,p1(2,:),p2(2,:));
clcd = x_des(that*100,p1(3,:),p2(3,:));

% calculate intermediate variables
phi   = atan(((1-a)/(1+ap)) * (R/(r*tsr)));
cd    = cl/clcd;
cy    = cl * cos(phi) + cd * sin(phi);
cx    = cl * sin(phi) - cd * cos(phi);
f     = (B/2) * ((R-r)/(r * sin(phi)));
F     = (2/pi) * acos(exp(-f));
sigma = (c*B) ./ (2*pi*r);

% calculate residuals
res_c  = 4*pi*r * sin(phi)^2 * F * ((2*a) / ...
    (cy * B * (1-a))) - c;

% res_a = 1/((4*F*sin(phi)*cos(phi))/(sigma*cx)-1) - a;
res_ap = 1 / ((4 * F * sin(phi) * cos(phi) / ...
    (sigma * cx)) - 1) - ap;

% pack output
out = [res_c, res_ap];

% optional outputs
if nargout == 2
    beta = phi - deg2rad(alpha);
    cp = ((1-a)^2 + (tsr*(r/R))^2 * (1+ap)^2) * tsr * (r/R) * sigma * cx;
    ct = ((1-a)^2 + (tsr*(r/R))^2 * (1+ap)^2) * sigma * cy;
%     result.a = 1 / (((4*F*np.sin(phi)**2) / (sigma * cy)) + 1)
    varargout{1} = [t,a,phi,alpha,beta,cl,cd,ap,cp,ct];
%     varargout{1} = [t,c,phi,alpha,beta,cl,cd,a,ap,cp,ct];
end
end