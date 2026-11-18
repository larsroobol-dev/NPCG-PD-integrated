function bankSign = anpcg_bank_direction(x,sigmaMag,previousSign)
% Predict final crossrange for both possible bank directions.
% Keep the current direction unless its predicted crossrange
% is more than 7 times that of the opposite direction.

%% Parameters
P.R     = 3397e3;
P.mu    = 42828e9;
P.omega = 7.088e-5;

P.mass = 58800;
P.S    = 160;
P.CL   = 0.5236;
P.CD   = 0.9696;

P.K  = 1.28;
P.e0 = P.mu/x(1)-0.5*x(4)^2;
P.ef = P.mu/(P.R+2480)-0.5*450^2;

bankSign = previousSign;

if P.e0 >= P.ef
    return;
end

%% Normal to the initial-to-target great-circle plane
rInitial = position_unit(deg2rad(-176.4),deg2rad(-21.3));
rTarget  = position_unit(deg2rad(-175.8),deg2rad(0.276));

normal = cross(rInitial,rTarget);
normal = normal/norm(normal);

%% Predict with current and opposite banking directions
x0 = x(1:6);
x0 = x0(:);

crossCurrent = predict_crossrange( ...
    x0,sigmaMag,previousSign,P,normal);

crossOpposite = predict_crossrange( ...
    x0,sigmaMag,-previousSign,P,normal);

%% Bank-reversal threshold
if crossCurrent > 7*max(crossOpposite,1e-12)
    bankSign = -previousSign;
end

end


function magnitude = predict_crossrange(x0,sigmaMag,bankSign,P,normal)

options = odeset( ...
    'RelTol',1e-7, ...
    'AbsTol',1e-8, ...
    'MaxStep',5, ...
    'Events',@(t,x) lateral_energy_event(t,x,P));

[~,~,te,xe] = ode45( ...
    @(t,x) lateral_dynamics(t,x,sigmaMag,bankSign,P), ...
    [0 1500],x0,options);

if isempty(te)
    error('Lateral predictor did not reach terminal energy.');
end

xf = xe(end,:).';
rFinal = position_unit(xf(2),xf(3));

projection = dot(normal,rFinal);
projection = max(-1,min(1,projection));
magnitude = abs(asin(projection));

end


function dx = lateral_dynamics(~,x,sigmaMag,bankSign,P)

r     = x(1);
phi   = x(3);
V     = x(4);
gamma = x(5);
psi   = x(6);

%% Predicted bank profile
e = P.mu/r-0.5*V^2;
progress = (e-P.e0)/(P.ef-P.e0);

sigma = bankSign*2*sigmaMag/(1+exp(P.K*progress));

%% Atmosphere
h = r-P.R;

T = 1.4e-13*h^3-8.85e-9*h^2 ...
    -1.245e-3*h+205.36;

rho = (559.35/(188.95*T))*exp(-0.000105*h);

%% Accelerations
q = 0.5*rho*V^2;
D = q*P.CD*P.S/P.mass;
L = q*P.CL*P.S/P.mass;
g = P.mu/r^2;

omega = P.omega;

%% Full rotating-planet equations
dx = zeros(6,1);

dx(1) = V*sin(gamma);

dx(2) = V*cos(gamma)*sin(psi)/(r*cos(phi));

dx(3) = V*cos(gamma)*cos(psi)/r;

dx(4) = -D-g*sin(gamma) ...
    + omega^2*r*cos(phi)*(sin(gamma)*cos(phi) ...
    - cos(gamma)*sin(phi)*cos(psi));

dx(5) = L*cos(sigma)/V ...
    + (V/r-g/V)*cos(gamma) ...
    + 2*omega*cos(phi)*sin(psi) ...
    + (omega^2*r/V)*cos(phi)*(cos(gamma)*cos(phi) ...
    + sin(gamma)*sin(phi)*cos(psi));

dx(6) = L*sin(sigma)/(V*cos(gamma)) ...
    + (V/r)*cos(gamma)*sin(psi)*tan(phi) ...
    + 2*omega*(sin(phi)-tan(gamma)*cos(phi)*cos(psi)) ...
    + (omega^2*r/(V*cos(gamma))) ...
    *sin(phi)*cos(phi)*sin(psi);


end


function [value,isterminal,direction] = lateral_energy_event(~,x,P)

value = P.mu/x(1)-0.5*x(4)^2-P.ef;
isterminal = 1;
direction = 1;

end


function u = position_unit(theta,phi)

u = [cos(phi)*cos(theta);
     cos(phi)*sin(theta);
     sin(phi)];

end