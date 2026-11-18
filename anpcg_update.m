function sigma0 = anpcg_update(xCurrent,sigmaGuess)
% Compute the current bank-angle magnitude, in radians.

%% Mission parameters
P.R     = 3397e3;
P.mu    = 42828e9;
P.mass  = 58800;
P.S     = 160;
P.CL    = 0.5236;
P.CD    = 0.9696;

P.K       = 1.28;
P.hTarget = 2480;
P.VTarget = 450;

thetaTarget = deg2rad(-175.8);
phiTarget   = deg2rad(0.276);

%% Current vehicle state
r     = xCurrent(1);
theta = xCurrent(2);
phi   = xCurrent(3);
V     = xCurrent(4);
gamma = xCurrent(5);

%% Remaining range
cosAngle = sin(phi)*sin(phiTarget) ...
    + cos(phi)*cos(phiTarget)*cos(thetaTarget-theta);

cosAngle = max(-1,min(1,cosAngle));
rangeToGo = P.R*acos(cosAngle);

%% Energy and prediction initial conditions
P.e0 = P.mu/r-V^2/2;
P.ef = P.mu/(P.R+P.hTarget)-P.VTarget^2/2;

sigma0 = sigmaGuess;

if P.e0 >= P.ef
    return;
end

y0 = [r; V; gamma; 0];

%% Predictor-corrector
deltaSigma = deg2rad(0.1);

for k = 1:15

    z = predict_range(sigma0,y0,P,rangeToGo);

    if abs(z) < 100
        return;
    end

    zPlus = predict_range( ...
        sigma0+deltaSigma,y0,P,rangeToGo);

    zMinus = predict_range( ...
        sigma0-deltaSigma,y0,P,rangeToGo);

    dz = (zPlus-zMinus)/(2*deltaSigma);

    if ~isfinite(dz) || abs(dz) < 1e-6
        error('ANPCG: invalid range sensitivity.');
    end

    %% Safeguarded Newton update
    newtonStep = -z/dz;
    stepScale = 1;
    accepted = false;

    for backtrack = 1:21

        sigmaTrial = sigma0 + stepScale*newtonStep;

        % Only evaluate physically allowed bank magnitudes
        if isfinite(sigmaTrial) && sigmaTrial >= 0 && sigmaTrial <= pi

            zTrial = predict_range( ...
                sigmaTrial,y0,P,rangeToGo);

            % Accept only if the range error decreases
            if isfinite(zTrial) && abs(zTrial) < abs(z)
                sigma0 = sigmaTrial;
                accepted = true;
                break;
            end
        end

        % Shorten the step and try again
        stepScale = stepScale/2;
    end

    if ~accepted

        % Sample the allowed bank-angle interval
        bankGridDeg = (0:5:180).';
        errorGrid = zeros(size(bankGridDeg));

        for j = 1:numel(bankGridDeg)
            errorGrid(j) = predict_range( ...
                deg2rad(bankGridDeg(j)),y0,P,rangeToGo);
        end

        [bestAbsError,bestIndex] = min(abs(errorGrid));

        fprintf('\nFailed-update diagnostic:\n');
        fprintf('Altitude: %.3f m\n',r-P.R);
        fprintf('Velocity: %.3f m/s\n',V);
        fprintf('Current bank: %.6f deg\n',rad2deg(sigma0));
        fprintf('Current range error: %.6f m\n',z);
        fprintf('Best sampled bank: %.3f deg\n', ...
            bankGridDeg(bestIndex));
        fprintf('Smallest sampled absolute error: %.6f m\n', ...
            bestAbsError);
        fprintf('Sampled residual range: %.6f to %.6f m\n', ...
            min(errorGrid),max(errorGrid));

        save('ANPCG_backtracking_failure.mat', ...
            'xCurrent','sigmaGuess','sigma0','z','dz', ...
            'rangeToGo','P','bankGridDeg','errorGrid');

        error('ANPCG stopped; diagnostic saved for inspection.');
    end

    % Check the accepted correction immediately
    if abs(zTrial) < 100
        return;
    end
end

error('ANPCG failed to converge within 15 iterations.');

end

function [z,yFinal] = predict_range(sigma0,y0,P,rangeToGo)

    options = odeset( ...
        'RelTol',1e-8, ...
        'AbsTol',1e-8, ...
        'MaxStep',2, ...
        'Events',@(t,y) final_energy_event(t,y,P));

    [~,~,te,ye] = ode45( ...
        @(t,y) predictor_dynamics(t,y,sigma0,P), ...
        [0 1500],y0,options);

    if isempty(te)
        error('Predictor did not reach the terminal energy.');
    end

    yFinal = ye(end,:).';
    z = yFinal(4)-rangeToGo;
end

function dy = predictor_dynamics(~,y,sigma0,P)

    r     = y(1);
    V     = y(2);
    gamma = y(3);
    h     = r-P.R;

    %% Atmosphere: h in meters
    T = 1.4e-13*h^3 - 8.85e-9*h^2 ...
        - 1.245e-3*h + 205.36;

    rho = (559.35/(188.95*T))*exp(-0.000105*h);

    %% Accelerations
    q = 0.5*rho*V^2;
    D = q*P.CD*P.S/P.mass;
    L = q*P.CL*P.S/P.mass;
    g = P.mu/r^2;

    %% Logistic bank profile: Equation (10)
    e = P.mu/r - V^2/2;
    progress = (e-P.e0)/(P.ef-P.e0);

    sigma = 2*sigma0/(1+exp(P.K*progress));

    %% Reduced longitudinal prediction model
    dy = zeros(4,1);

    dy(1) = V*sin(gamma);
    dy(2) = -D-g*sin(gamma);
    dy(3) = L*cos(sigma)/V ...
          + (V/r-g/V)*cos(gamma);
    dy(4) = V*cos(gamma);
end

function [value,isterminal,direction] = final_energy_event(~,y,P)

    e = P.mu/y(1)-y(2)^2/2;

    value      = e-P.ef;
    isterminal = 1;
    direction  = 1;
end