function [F, Fcomponents, phantom] = evaluateFarfieldNystrom(objects, params, q, kappa, n2, plotflag)
% [F,geometry] = computefarfield(object,kappa,n,plotflag)
%
%  COMPUTEFARFIELD: Solves the direct scattering problem for the
%       Helmholtz equation with sound soft scatterers, sound hard 
%       scatterers, and scatterers carrying an impedance boundary condition
%       and plane wave incident fields. 
%       It's the Nystroem method described in Colton&Kress (1998) and 
%       Kress (1995).
%       Based on codes by Andreas Kirsch.
%
%  INPUT: kappa -- wave number
%         n -- no. of discretization points
%         plotflag -- binary (0 or 1), if plotflag = 1, then geometry and
%                     fields are plotted
%  OUTPUT: F -- farfield matrix
%          Fcomponents -- cell array containing far field radiated by each
%                         scatterer (taking into account multiple scattering)
%          geometry -- geometry data
%
%  written by R. Griesmaier (2011), 
%  further development of a code for a single sound-soft scatterer by A. Kirsch


% kappa = 1;   % wavenumber
%n2 = 128; % discretization of unit sphere S (i.e., number of directions)
%object = 'extended_obstacles';  % select of obstacles
eta = kappa;  % coupling parmeter in the integral equation

% the different obstacles are parameterized by the first two entries of
% the variable par. par(5) adds a rotation (mathematically positive direction),
% par(6) is a scaling factor for the obstacle, and (par(3),par(4)) adds a 
% final shift of the obstacle.
% the generation of the obstacles requires file "kurve.m"

% switch object
%     case 'circle'
%         objects = {'circle'};
%         nr_of_scatterers = 1;
%         params = [1 1 -2.5 -2.5 0 1]; % Ellipse
%         q=[-0.5];
%     case 'nut'
%         objects = {'nut'};
%         nr_of_scatterers = 1;
%         params = [1 4 3 3 20 1]; % Erdnuss
%         q=[1];
%     case 'kite'
%         objects = {'kite'};
%         nr_of_scatterers = 1;
%         params = [1.5 .65 -6 -6 0 1]; % Drachen
%         q=[1];
%     case 'kite2'
%         objects = {'kite'};
%         nr_of_scatterers = 1;
%         params = [1.5 0.65 4.5 4.5 45 1 ]; % Drachen
%         q=[-0.5];
%     case 'circle2'
%         objects = {'circle'};
%         nr_of_scatterers = 1;
%         params = [0.65 1.5 -2 -2 45 1]; % Ellipse
%         q=[0.5];
%     case '2circles'
%         objects = {'circle', 'circle'};
%         nr_of_scatterers = 2;
%         params = [1 2 -10 -10 0 0.5 ; 2 1 12 16 0 0.5];
%         q=[-0.5,0.5];
%     case 'kiteandcircle'
%         objects = {'kite', 'circle'};
%         nr_of_scatterers = 2;
%         params = [1.5 0.65 4.5 4.5 45 1 ; 0.65 1.5 -2 -2 45 1];
%         q=[-0.5,0.5];
%     case '3scatterers'
%         objects = {'kite', 'nut', 'circle'};
%         nr_of_scatterers = 3;
%         params = [1.5 .65 -6 12 70 1; 1 4 12 -8 -28 1; 2 0.5 -9 -9 135 1];
%         q=[-0.9,0.5,0.5];
%     case '3obstacles_wide'
%         objects = {'kite', 'nut', 'circle'};
%         nr_of_scatterers = 3;
%         params = [1.5 .65 26 -3 45 3; 1 4 -20 23 30 3; 2.8 1 -15 -21 135 3];
%         q=[1, 1, 1];
%     case '2obstacles_wide'
%         objects = {'kite', 'circle'};
%         nr_of_scatterers = 2;
%         params = [1 .5 26 -3 45 3; 1 .5 -15 -21 135 3];
%         q=[.25, .25];
%     case '4obstacles_wide'
%         objects = {'kite', 'nut', 'circle','circle'};
%         nr_of_scatterers = 4;
%         params = [1.5 .65 36 -3 45 2; 1 4 -15 18 30 2; 2.8 1 -15 -21 135 2; 1 2 25 35 -30 1.5];
%         q=[1,1,1,1];
%     case 'kiteandellipse'
%         objects = {'kite', 'circle'};
%         nr_of_scatterers = 2;
%         params = [1.5 .65 26 -3 45 3; 2.8 1 -15 -21 135 3];
%         R = [5, 4];  % radii of circles containing source components
%         z = [[24; -4], [-15; -21]];
%         % params = [1.5 0.65 4.5 4.5 45 1 ; 0.65 1.5 -2 -2 45 1];
%         q = [.25,.25];
%     case 'nut_paper'
%         objects = {'nut'};
%         nr_of_scatterers = 1;
%         params = [1 .5 26 -3 45 3];
%         R = [4.5];
%         z = [[26; -3]];
%         q = [.25];
%     case 'kite_paper'
%         objects = {'kite'};
%         nr_of_scatterers = 1;
%         params = [1 .5 -15 -21 135 3];
%         R = [3.5];
%         z = [[-15; -21]];
%         q = [.25];
%     case 'nut_and_kite_paper'
%         objects = {'nut', 'kite'};
%         nr_of_scatterers = 2;
%         params = [1 .5 26 -3 45 3; 1 .5 -15 -21 135 3];
%         R = [4.5, 3.5];
%         z = [[26; -3], [-15; -21]];
%         q = [.25,.25];
% end
nr_of_scatterers = length(objects);

n2m = n2 - 1;
n = n2/2;
pin = pi/n;

scatterers = cell(1, nr_of_scatterers);
for iteri = 1:nr_of_scatterers
    object = objects{iteri};
    param = params(iteri,:);
    [x, dx, ddx, dddx] = kurve(n2, object, param);
    scatterers{iteri} = scatterer(n, n2, n2m, pin, x, dx, ddx, dddx, q(iteri));
end

IOP = zeros(2*nr_of_scatterers*n2);
for iteri = 1 : nr_of_scatterers
        obji = scatterers{iteri};
        
        for iterj = 1 : nr_of_scatterers
            objj = scatterers{iterj};
            kappa0 = sqrt(1+objj.q)*kappa;
            if iterj == iteri
                IOPij11 = DirichletOnScattererPlus(kappa, eta, obji);
                IOPij12 = -DirichletOnScattererMinus(kappa0, eta, obji);
                IOPij21 = NeumannOnScattererPlus(kappa, eta, obji);
                IOPij22 = -NeumannOnScattererMinus(kappa0, eta, obji);                
            else
                IOPij11 = DirichletOffScatterer(kappa, eta, objj, obji);
                IOPij12 = zeros(n2);
                IOPij21 = NeumannOffScatterer(kappa, eta, objj, obji);
                IOPij22 = zeros(n2);
            end
            
            IOPij = [IOPij11 IOPij12; IOPij21 IOPij22];
            
            IOP(2*(iteri-1)*n2+1 : 2*iteri*n2, 2*(iterj-1)*n2+1 : 2*iterj*n2) = IOPij;
        end
end
clear IOPij IOPij11 IOPij12 IOPij21 IOPij22

cs = cos((0:n2m)*pin);
sn = sin((0:n2m)*pin);

RHS = zeros(2*nr_of_scatterers*n2, n2);
for iteri = 1 : nr_of_scatterers
    
        RHSi1 = EvalRHSDirichlet(kappa, cs, sn, scatterers{iteri});
        RHSi2 = EvalRHSNeumann(kappa, cs, sn, scatterers{iteri});
        RHS(2*(iteri-1)*n2+1 : (2*iteri-1)*n2, :) = RHSi1;
        RHS((2*iteri-1)*n2+1 : 2*iteri*n2, :) = RHSi2;
        
end
clear RHSi1 RHSi2;

PHI = IOP \ RHS;

for iteri = 1 : nr_of_scatterers
    scatterers{iteri}.PHI = PHI(2*(iteri-1)*n2+1 : (2*iteri-1)*n2,:);
end
% clear PHI;

sigma=-1i;  % sigma = (1-1i)/(4*sqrt(kappa*pi));
F = zeros(n2);
Fcomponents = cell(1, nr_of_scatterers);
for iteri = 1 : nr_of_scatterers
    
    Fhelp = EvalFarField(kappa, eta, cs, sn, sigma, scatterers{iteri});
    Fcomponents{iteri} = [Fhelp, Fhelp(:,1)];
    F = F + Fhelp;
    
end
clear Fhelp

phantom = zeros(n2, nr_of_scatterers*2);
for iteri = 1 : nr_of_scatterers
   phantom(:, (iteri-1)*2+1 : iteri*2) = scatterers{iteri}.x; 
end

phantom(n2+1,:) = phantom(1,:);

if plotflag == 1
    figure()
    hold on
    for iter = 1 : nr_of_scatterers
        plot(phantom(:,2*iter-1), phantom(:,2*iter))
    end
    hold off
    axis equal
    grid on
    title('Geometry of the scatterer')
    title('phantom')
end

FFloc = [ F F(:,1) ];
FFloc = [ FFloc ; FFloc(1,:) ];
xord = linspace(0,360,max(size(FFloc)));
yord = linspace(0,360,max(size(FFloc)));

if plotflag == 1
    figure()
    subplot(1,2,1)
    imagesc(xord,yord,real(FFloc))
    axis xy
    axis square
    title('farfield (real)')
    subplot(1,2,2)
    imagesc(xord,yord,imag(FFloc))
    axis xy
    axis square
    title('farfield (imag)')
end


end


function obj = scatterer(n, n2, n2m, pin, x, dx, ddx, dddx, q)
%
obj.n = n;
obj.n2 = n2;
obj.n2m = n2m;
obj.pin = pin;
obj.x = x;
obj.dx = dx;
obj.ddx = ddx;
obj.dddx = dddx;
obj.q = q;
end


function IOP = DirichletOnScattererPlus(kappa, eta, obj)
%
wR = (-1).^(0:obj.n2m)' / obj.n2 + ...
     cos((0:obj.n2m)'*(1:obj.n-1)*obj.pin)*(1./(1:obj.n-1)');
wR = -2*pi*wR/obj.n;
WR = toeplitz(wR);

normdx = sqrt( obj.dx(:,1).^2 + obj.dx(:,2).^2 );
r =     (obj.x(:,1)*ones(1,obj.n2)-ones(obj.n2,1)*obj.x(:,1)').^2;
r = r + (obj.x(:,2)*ones(1,obj.n2)-ones(obj.n2,1)*obj.x(:,2)').^2;
r = sqrt(r) + eye(obj.n2);

J0 = besselj(0,kappa*(r-eye(obj.n2)));
Y0 = bessely(0,kappa*r);
J1 = besselj(1,kappa*r);
Y1 = bessely(1,kappa*r);

xcrossdx = obj.x(:,1).*obj.dx(:,2)  - obj.x(:,2).*obj.dx(:,1);
Hhelp = ones(obj.n2,1)*xcrossdx' - obj.x(:,1)*obj.dx(:,2)' + ...
        obj.x(:,2)*obj.dx(:,1)';
Hhelp = kappa*Hhelp./r;

IJ = toeplitz((0:obj.n2m));
SN = sin(obj.pin*.5*IJ);

L = 0.5 * (1i*J1 - Y1) .* Hhelp;
L1 = -0.5 * J1/pi .* Hhelp;
L2 = L - log(4*SN.*SN+eye(obj.n2)) .* L1;
L = WR.*L1 + obj.pin*L2;
diagL = (0.5/pi)*( obj.dx(:,1).*obj.ddx(:,2) - ...
        obj.dx(:,2).*obj.ddx(:,1) )./(normdx.^2);
L = L + diag(obj.pin*diagL - diag(L));

Normdx = repmat(normdx, 1, obj.n2);
M = 0.5 * (1i*J0 - Y0) .* Normdx';
M1 = -0.5 * J0/pi .* Normdx';
M2 = M - log(4*SN.*SN+eye(obj.n2)) .* M1;
diagM2 = (1i/2 + psi(1)/pi - (0.5/pi)*log(0.25*kappa^2*normdx.^2)) .* normdx;
M2 = M2 + diag(diagM2 - diag(M2));
M = WR.*M1 + obj.pin*M2;

IOP = eye(obj.n2) - L - 1i*eta*M;
end


function IOP = DirichletOnScattererMinus(kappa, eta, obj)
%
wR = (-1).^(0:obj.n2m)' / obj.n2 + ...
     cos((0:obj.n2m)'*(1:obj.n-1)*obj.pin)*(1./(1:obj.n-1)');
wR = -2*pi*wR/obj.n;
WR = toeplitz(wR);

normdx = sqrt( obj.dx(:,1).^2 + obj.dx(:,2).^2 );
r =     (obj.x(:,1)*ones(1,obj.n2)-ones(obj.n2,1)*obj.x(:,1)').^2;
r = r + (obj.x(:,2)*ones(1,obj.n2)-ones(obj.n2,1)*obj.x(:,2)').^2;
r = sqrt(r) + eye(obj.n2);

J0 = besselj(0,kappa*(r-eye(obj.n2)));
Y0 = bessely(0,kappa*r);
J1 = besselj(1,kappa*r);
Y1 = bessely(1,kappa*r);

xcrossdx = obj.x(:,1).*obj.dx(:,2)  - obj.x(:,2).*obj.dx(:,1);
Hhelp = ones(obj.n2,1)*xcrossdx' - obj.x(:,1)*obj.dx(:,2)' + ...
        obj.x(:,2)*obj.dx(:,1)';
Hhelp = kappa*Hhelp./r;

IJ = toeplitz((0:obj.n2m));
SN = sin(obj.pin*.5*IJ);

L = 0.5 * (1i*J1 - Y1) .* Hhelp;
L1 = -0.5 * J1/pi .* Hhelp;
L2 = L - log(4*SN.*SN+eye(obj.n2)) .* L1;
L = WR.*L1 + obj.pin*L2;
diagL = (0.5/pi)*( obj.dx(:,1).*obj.ddx(:,2) - ...
        obj.dx(:,2).*obj.ddx(:,1) )./(normdx.^2);
L = L + diag(obj.pin*diagL - diag(L));

Normdx = repmat(normdx, 1, obj.n2);
M = 0.5 * (1i*J0 - Y0) .* Normdx';
M1 = -0.5 * J0/pi .* Normdx';
M2 = M - log(4*SN.*SN+eye(obj.n2)) .* M1;
diagM2 = (1i/2 + psi(1)/pi - (0.5/pi)*log(0.25*kappa^2*normdx.^2)) .* normdx;
M2 = M2 + diag(diagM2 - diag(M2));
M = WR.*M1 + obj.pin*M2;

IOP = -eye(obj.n2) - L - 1i*eta*M;
end


function IOP = DirichletOffScatterer(kappa, eta, obj1, obj2)
% Feld on obj1 (integral ueber obj1) auf obj2 ausgeewertet
normdx = sqrt( obj1.dx(:,1).^2 + obj1.dx(:,2).^2 );
r = (repmat(obj2.x(:,1), 1, obj1.n2) - repmat(obj1.x(:,1)', obj2.n2, 1)).^2;
r = r + (repmat(obj2.x(:,2), 1, obj1.n2) - ...
         repmat(obj1.x(:,2)', obj2.n2, 1)).^2;
r = sqrt(r);

J0 = besselj(0,kappa*r);
Y0 = bessely(0,kappa*r);
J1 = besselj(1,kappa*r);
Y1 = bessely(1,kappa*r);

xcrossdx = obj1.x(:,1).*obj1.dx(:,2)  - obj1.x(:,2).*obj1.dx(:,1);
Hhelp = repmat(xcrossdx', obj2.n2, 1) - obj2.x(:,1)*obj1.dx(:,2)' + ...
           obj2.x(:,2)*obj1.dx(:,1)';
Hhelp = kappa * Hhelp ./ r;

Normdx = repmat(normdx', obj2.n2, 1);
IOP = obj1.pin * (-0.5 * (1i*J1 - Y1) .* Hhelp - ...
        1i*eta * 0.5 * (1i*J0 - Y0) .* Normdx);
end


function IOP = NeumannOnScattererPlus(kappa, eta, obj)
%
IJ = toeplitz((0:obj.n2m));
SN = sin(obj.pin*.5*IJ);
clear IJ

wt = zeros(obj.n2,1);
wt(1) = -obj.n/2;
wt(2:2:end) = 1 ./ (obj.n2 * (SN(2:2:end,1)).^2);
WT = toeplitz(wt);

wr = cos((0:obj.n2m)'*(1:obj.n-1)*obj.pin)*(1./(1:obj.n-1)') + ...
        (-1).^(0:obj.n2m)' / obj.n2;
wr = -2*pi*wr/obj.n;
WR = toeplitz(wr);

normdx = sqrt( obj.dx(:,1).^2 + obj.dx(:,2).^2 );
r =     (repmat(obj.x(:,1),1,obj.n2)-repmat(obj.x(:,1)',obj.n2,1)).^2;
r = r + (repmat(obj.x(:,2),1,obj.n2)-repmat(obj.x(:,2)',obj.n2,1)).^2;
r = sqrt(r) + eye(obj.n2);

J0 = besselj(0,kappa*r);
Y0 = bessely(0,kappa*r);
J1 = besselj(1,kappa*r);
Y1 = bessely(1,kappa*r);

xcrossdx = obj.x(:,1).*obj.dx(:,2) - obj.x(:,2).*obj.dx(:,1);
Hhelp = -repmat(xcrossdx,1,obj.n2) + obj.dx(:,2)*(obj.x(:,1)') - ...
         obj.dx(:,1)*(obj.x(:,2)');
Hhelp = kappa * Hhelp ./ r;
Hhelp = Hhelp .* repmat(normdx',obj.n2,1);

H = 0.5 * Hhelp .* (1i*J1 - Y1);
H1 = -(0.5/pi) * Hhelp .* J1;
H1 = H1 - diag(diag(H1));
H2 = H - log(4*SN.*SN+eye(obj.n2)) .* H1;
diagH2 = (0.5/pi) * (obj.dx(:,2).*obj.ddx(:,1) - ...
                     obj.dx(:,1).*obj.ddx(:,2)) ./ normdx;
H2 = H2 + diag(diagH2 - diag(H2));
H = WR.*H1 + obj.pin*H2;
clear Hhelp H1 H2 diagH2

M = 0.5 * (1i*J0 - Y0);
M1 = -(0.5/pi) * J0;
diagM1 = -(0.5/pi) * ones(obj.n2,1);
M1 = M1 + diag(diagM1 - diag(M1));
M2 = M - log(4*SN.*SN+eye(obj.n2)) .* M1;
diagM2 = 1i/2 + psi(1)/pi - log(0.5*kappa*normdx)/pi;
M2 = M2 + diag(diagM2 - diag(M2));
M = WR.*M1 + obj.pin*M2;
clear M1 M2 diagM1 diagM2

xdotdx = obj.x(:,1).*obj.dx(:,1) + obj.x(:,2).*obj.dx(:,2);
Ntilde = repmat(xdotdx,1,obj.n2) - obj.dx(:,1)*(obj.x(:,1)') - ...
                obj.dx(:,2)*(obj.x(:,2)');
Ntilde = Ntilde .* (-Ntilde') ./ (r.^2);

Dxdotdx = obj.dx(:,1)*(obj.dx(:,1)') + obj.dx(:,2)*(obj.dx(:,2)');
N = kappa^2 * (1i*J0 - Y0) - 2*kappa * (1i*J1 - Y1) ./ r;
N = 0.5 * Ntilde .* N;
N = N + 0.5*kappa * Dxdotdx .* (1i*J1 - Y1) ./ r;
N = N + 0.25 ./ (pi*SN.*SN+eye(obj.n2));

N1 = kappa^2 * J0 - 2*kappa * J1 ./ r;
N1 = -(0.5/pi) * Ntilde .* N1;
N1 = N1 - (0.5*kappa/pi) * Dxdotdx .* J1 ./ r;
diagN1 = -(0.25*kappa^2/pi) * normdx.^2;
N1 = N1 + diag(diagN1 - diag(N1));

N2 = N - log(4*SN.*SN+eye(obj.n2)) .* N1;
dxdotddx = obj.dx(:,1).*obj.ddx(:,1) + obj.dx(:,2).*obj.ddx(:,2);
normddx = sqrt( obj.ddx(:,1).^2 + obj.ddx(:,2).^2 );
dxdotdddx = obj.dx(:,1).*obj.dddx(:,1) + obj.dx(:,2).*obj.dddx(:,2);
diagN2 = (pi*1i - 1 + 2*psi(1) - 2*log(0.5*kappa*normdx));
diagN2 = (0.25*kappa^2/pi) * diagN2 .*  normdx.^2 + 1/(12*pi);
diagN2 = diagN2 + dxdotddx.^2 ./ (2*pi*normdx.^4);
diagN2 = diagN2 - normddx.^2 ./ (4*pi*normdx.^2);
diagN2 = diagN2 - dxdotdddx ./ (6*pi*normdx.^2);
N2 = N2 + diag(diagN2 - diag(N2));

N = WR.*N1 + obj.pin*N2;
clear Ntilde N1 N2 diagN1 diagN2 WR

K = kappa^2 * M .* Dxdotdx - N - 1i*eta*H;
a = 1i * eta * normdx;
IOP = diag(a) + WT + K;
end


function IOP = NeumannOnScattererMinus(kappa, eta, obj)
%
IJ = toeplitz((0:obj.n2m));
SN = sin(obj.pin*.5*IJ);
clear IJ

wt = zeros(obj.n2,1);
wt(1) = -obj.n/2;
wt(2:2:end) = 1 ./ (obj.n2 * (SN(2:2:end,1)).^2);
WT = toeplitz(wt);

wr = cos((0:obj.n2m)'*(1:obj.n-1)*obj.pin)*(1./(1:obj.n-1)') + ...
        (-1).^(0:obj.n2m)' / obj.n2;
wr = -2*pi*wr/obj.n;
WR = toeplitz(wr);

normdx = sqrt( obj.dx(:,1).^2 + obj.dx(:,2).^2 );
r =     (repmat(obj.x(:,1),1,obj.n2)-repmat(obj.x(:,1)',obj.n2,1)).^2;
r = r + (repmat(obj.x(:,2),1,obj.n2)-repmat(obj.x(:,2)',obj.n2,1)).^2;
r = sqrt(r) + eye(obj.n2);

J0 = besselj(0,kappa*r);
Y0 = bessely(0,kappa*r);
J1 = besselj(1,kappa*r);
Y1 = bessely(1,kappa*r);

xcrossdx = obj.x(:,1).*obj.dx(:,2) - obj.x(:,2).*obj.dx(:,1);
Hhelp = -repmat(xcrossdx,1,obj.n2) + obj.dx(:,2)*(obj.x(:,1)') - ...
         obj.dx(:,1)*(obj.x(:,2)');
Hhelp = kappa * Hhelp ./ r;
Hhelp = Hhelp .* repmat(normdx',obj.n2,1);

H = 0.5 * Hhelp .* (1i*J1 - Y1);
H1 = -(0.5/pi) * Hhelp .* J1;
H1 = H1 - diag(diag(H1));
H2 = H - log(4*SN.*SN+eye(obj.n2)) .* H1;
diagH2 = (0.5/pi) * (obj.dx(:,2).*obj.ddx(:,1) - ...
                     obj.dx(:,1).*obj.ddx(:,2)) ./ normdx;
H2 = H2 + diag(diagH2 - diag(H2));
H = WR.*H1 + obj.pin*H2;
clear Hhelp H1 H2 diagH2

M = 0.5 * (1i*J0 - Y0);
M1 = -(0.5/pi) * J0;
diagM1 = -(0.5/pi) * ones(obj.n2,1);
M1 = M1 + diag(diagM1 - diag(M1));
M2 = M - log(4*SN.*SN+eye(obj.n2)) .* M1;
diagM2 = 1i/2 + psi(1)/pi - log(0.5*kappa*normdx)/pi;
M2 = M2 + diag(diagM2 - diag(M2));
M = WR.*M1 + obj.pin*M2;
clear M1 M2 diagM1 diagM2

xdotdx = obj.x(:,1).*obj.dx(:,1) + obj.x(:,2).*obj.dx(:,2);
Ntilde = repmat(xdotdx,1,obj.n2) - obj.dx(:,1)*(obj.x(:,1)') - ...
                obj.dx(:,2)*(obj.x(:,2)');
Ntilde = Ntilde .* (-Ntilde') ./ (r.^2);

Dxdotdx = obj.dx(:,1)*(obj.dx(:,1)') + obj.dx(:,2)*(obj.dx(:,2)');
N = kappa^2 * (1i*J0 - Y0) - 2*kappa * (1i*J1 - Y1) ./ r;
N = 0.5 * Ntilde .* N;
N = N + 0.5*kappa * Dxdotdx .* (1i*J1 - Y1) ./ r;
N = N + 0.25 ./ (pi*SN.*SN+eye(obj.n2));

N1 = kappa^2 * J0 - 2*kappa * J1 ./ r;
N1 = -(0.5/pi) * Ntilde .* N1;
N1 = N1 - (0.5*kappa/pi) * Dxdotdx .* J1 ./ r;
diagN1 = -(0.25*kappa^2/pi) * normdx.^2;
N1 = N1 + diag(diagN1 - diag(N1));

N2 = N - log(4*SN.*SN+eye(obj.n2)) .* N1;
dxdotddx = obj.dx(:,1).*obj.ddx(:,1) + obj.dx(:,2).*obj.ddx(:,2);
normddx = sqrt( obj.ddx(:,1).^2 + obj.ddx(:,2).^2 );
dxdotdddx = obj.dx(:,1).*obj.dddx(:,1) + obj.dx(:,2).*obj.dddx(:,2);
diagN2 = (pi*1i - 1 + 2*psi(1) - 2*log(0.5*kappa*normdx));
diagN2 = (0.25*kappa^2/pi) * diagN2 .*  normdx.^2 + 1/(12*pi);
diagN2 = diagN2 + dxdotddx.^2 ./ (2*pi*normdx.^4);
diagN2 = diagN2 - normddx.^2 ./ (4*pi*normdx.^2);
diagN2 = diagN2 - dxdotdddx ./ (6*pi*normdx.^2);
N2 = N2 + diag(diagN2 - diag(N2));

N = WR.*N1 + obj.pin*N2;
clear Ntilde N1 N2 diagN1 diagN2 WR

K = kappa^2 * M .* Dxdotdx - N - 1i*eta*H;
a = 1i * eta * normdx;
IOP = -diag(a) + WT + K;
end


function IOP = NeumannOffScatterer(kappa, eta, obj1, obj2)
% Feld on obj1 (integral ueber obj1) auf obj2 ausgewertet
normdxobj1 = sqrt( obj1.dx(:,1).^2 + obj1.dx(:,2).^2 );
r =     (repmat(obj2.x(:,1),1,obj1.n2)-repmat(obj1.x(:,1)',obj2.n2,1)).^2;
r = r + (repmat(obj2.x(:,2),1,obj1.n2)-repmat(obj1.x(:,2)',obj2.n2,1)).^2;
r = sqrt(r);

J0 = besselj(0,kappa*r);
Y0 = bessely(0,kappa*r);
J1 = besselj(1,kappa*r);
Y1 = bessely(1,kappa*r);

xcrossdxobj2 = obj2.x(:,1).*obj2.dx(:,2) - obj2.x(:,2).*obj2.dx(:,1);
xcrossdxobj2 = -repmat(xcrossdxobj2,1,obj1.n2) + obj2.dx(:,2)*(obj1.x(:,1)')...
               - obj2.dx(:,1)*(obj1.x(:,2)');
Hhelp = kappa * xcrossdxobj2 ./ r;
Hhelp = Hhelp .* repmat(normdxobj1',obj2.n2,1);
H = 0.5 * obj1.pin * Hhelp .* (1i*J1 - Y1);

xcrossdxobj1 = obj1.x(:,1).*obj1.dx(:,2) - obj1.x(:,2).*obj1.dx(:,1);
xcrossdxobj1 = -repmat(xcrossdxobj1',obj2.n2,1) + ...
               obj2.x(:,1)*(obj1.dx(:,2)') - obj2.x(:,2)*(obj1.dx(:,1)');
Thelp1 = -(1i*J0 - Y0) + (1i*J1 - Y1) ./ (kappa * r);
Thelp1 = - kappa^2 * Thelp1 .* xcrossdxobj2 .* xcrossdxobj1 ./ r.^2;
Thelp2 = obj2.dx(:,1)*(obj1.dx(:,1)') + obj2.dx(:,2)*(obj1.dx(:,2)');
Thelp2 = Thelp2 ./ r + xcrossdxobj2 .* xcrossdxobj1 ./ r.^3;
Thelp2 = kappa * (1i*J1 - Y1) .* Thelp2;
T = 0.5 * obj1.pin * (-Thelp1 + Thelp2);

IOP = T - 1i*eta*H;
end


function RHS = EvalRHSDirichlet(kappa, cs, sn, obj)
%
RHS = -2*exp( 1i*kappa* ( obj.x(:,1)*cs + obj.x(:,2)*sn ) );
end


function RHS = EvalRHSNeumann(kappa, cs, sn, obj)
% 
Nudotd = obj.dx(:,2)*cs - obj.dx(:,1)*sn;
RHS = -2*1i*kappa * Nudotd .* exp(1i*kappa*(obj.x(:,1)*cs + obj.x(:,2)*sn));
end


function F = EvalFarField(kappa, eta, cs, sn, sigma, obj)
%
Nudotd = obj.dx(:,2)*cs - obj.dx(:,1)*sn;
normdx = sqrt( obj.dx(:,1).^2 + obj.dx(:,2).^2 );
Normdx = repmat(normdx, 1, obj.n2);
F = exp(1i*kappa*(obj.x(:,1)*cs + obj.x(:,2)*sn));
F = F .* (kappa*Nudotd + eta*Normdx);
F = sigma * obj.pin * F' * obj.PHI;
end


function x = psi(y)
if y == 1
    x = -0.577215664901532;
else
    error('psi undefined')
end
end




