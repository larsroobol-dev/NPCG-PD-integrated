%% Mars entry: run, plot, compare and save
runSimulation = true;  % true: simulate; false: reuse existing out
saveResults = true;
modelName = 'MarsEntry2';
stopTime = '500';

%% Run the model
projectFolder = fileparts(mfilename('fullpath'));
if isempty(projectFolder), projectFolder = pwd; end
if runSimulation
    modelFile = fullfile(projectFolder,[modelName '.slx']);
    assert(isfile(modelFile),'Model not found: %s',modelFile);
    addpath(projectFolder,'-begin');
    if bdIsLoaded(modelName)
        assert(strcmpi(get_param(modelName,'FileName'),modelFile), ...
            'A different copy of this model is open. Save and close it first.');
    end
    load_system(modelFile);
    out = sim(modelName,'StopTime',stopTime);
else
    assert(exist('out','var') == 1, ...
        'Load a saved out variable, or set runSimulation = true.');
end

%% Read logged data
[tState,xState] = readHistory(out.StateHistory);
[tBank,sigmaBank] = readHistory(out.BankHistory);
assert(size(xState,2) == 7 && size(sigmaBank,2) == 1, ...
    'Expected seven states and one bank-angle signal.');
assert(numel(tState) > 1 && tState(end) > tState(1), ...
    'Simulation did not advance. Check the stopping condition.');

%% Reporting parameters: must match the model; these do not configure it
P.R = 3397e3; P.mu = 42828e9; P.g0 = 9.8;
P.mass = 58800; P.S = 160; P.CL = 0.5236; P.CD = 0.9696;
P.targetLongitude = deg2rad(-175.8);
P.targetLatitude = deg2rad(0.276);
P.targetAltitude = 2480; P.targetSpeed = 450;
P.loadLimit = 4; P.pressureLimit = 13; P.heatLimit = 500;

%% Trajectory and aerodynamic quantities
r = xState(:,1); longitude = xState(:,2); latitude = xState(:,3);
speed = xState(:,4); altitude = r-P.R;
bankMagnitude = abs(rad2deg(sigmaBank));
temperature = 1.4e-13*altitude.^3 - 8.85e-9*altitude.^2 ...
    - 1.245e-3*altitude + 205.36;
assert(all(temperature > 0),'Atmosphere temperature must be positive.');
density = 559.35./(188.95*temperature).*exp(-0.000105*altitude);
pressure = 0.5*density.*speed.^2;
drag = pressure*P.CD*P.S/P.mass;
lift = pressure*P.CL*P.S/P.mass;
loadG = hypot(lift,drag)/P.g0;
pressureKPa = pressure/1000;
heatFlux = 5.3697e-8*density.^0.5.*speed.^3.15; % kW/m^2

% Great-circle distance to target, in metres.
a = sin((latitude-P.targetLatitude)/2).^2 ...
    + cos(latitude)*cos(P.targetLatitude) ...
    .*sin((longitude-P.targetLongitude)/2).^2;
rangeToTarget = 2*P.R*asin(sqrt(max(0,min(1,a))));
targetMiss = rangeToTarget(end);
energyTarget = P.mu/(P.R+P.targetAltitude)-0.5*P.targetSpeed^2;
energyResidual = P.mu/r(end)-0.5*speed(end)^2-energyTarget;
[peakLoad,iLoad] = max(loadG);
[peakPressure,iPressure] = max(pressureKPa);
[peakHeat,iHeat] = max(heatFlux);

%% Optional digitized reference: time (s), bank magnitude (deg)
referenceBank = [];
referenceFile = fullfile(projectFolder,'ANPCG_blue.csv');
if isfile(referenceFile)
    referenceBank = readmatrix(referenceFile);
    assert(size(referenceBank,2) >= 2,'Reference CSV needs two columns.');
    referenceBank = referenceBank(:,1:2);
    referenceBank = referenceBank(all(isfinite(referenceBank),2),:);
    if any(referenceBank(:,2) < 0)
        warning('Reference CSV has negative magnitudes; check its calibration.');
    end
end

%% Six plots corresponding to Fig. 8
fig = figure('Name','ANPCG validation','Color','w', ...
    'Position',[100 100 1250 740]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
nexttile;
stairs(tBank,bankMagnitude,'b','LineWidth',1.4, ...
    'DisplayName','Simulink ANPCG'); hold on;
if ~isempty(referenceBank)
    plot(referenceBank(:,1),referenceBank(:,2),'ko','MarkerSize',3, ...
        'DisplayName','Digitized paper ANPCG');
    legend('Location','best');
end
xlabel('Time (s)'); ylabel('|\sigma| (deg)');
title('Bank magnitude'); xlim([0 500]); grid on;

nexttile;
plot(speed/1000,altitude/1000,'b','LineWidth',1.4); hold on;
plot(P.targetSpeed/1000,P.targetAltitude/1000,'r*','MarkerSize',8);
xlabel('Velocity (km/s)'); ylabel('Altitude (km)');
title('Altitude versus velocity'); grid on;

nexttile;
plot(rad2deg(longitude),rad2deg(latitude),'b','LineWidth',1.4); hold on;
plot(rad2deg(P.targetLongitude),rad2deg(P.targetLatitude),'r*','MarkerSize',8);
xlabel('Longitude \theta (deg)'); ylabel('Latitude \phi (deg)');
title('Ground track'); grid on;

nexttile;
plot(tState,loadG,'b','LineWidth',1.4); hold on;
yline(P.loadLimit,'r:','LineWidth',1.2);
xlabel('Time (s)'); ylabel('Aerodynamic load (g)');
title('Aerodynamic load'); xlim([0 500]);
ylim([0 max(5,1.1*peakLoad)]); grid on;

nexttile;
plot(tState,pressureKPa,'b','LineWidth',1.4); hold on;
yline(P.pressureLimit,'r:','LineWidth',1.2);
xlabel('Time (s)'); ylabel('Dynamic pressure (kPa)');
title('Dynamic pressure'); xlim([0 500]);
ylim([0 max(15,1.1*peakPressure)]); grid on;

nexttile;
plot(tState,heatFlux,'b','LineWidth',1.4); hold on;
yline(P.heatLimit,'r:','LineWidth',1.2);
xlabel('Time (s)'); ylabel('Heat flux (kW/m^2)');
title('Aerodynamic heating'); xlim([0 500]);
ylim([0 max(550,1.1*peakHeat)]); grid on;
sgtitle('ANPCG - red stars: targets; red dotted lines: limits');

%% Paper initial conditions, targets, results and errors
% Initial range is calculated from the paper's coordinates and Mars radius.
phi0 = deg2rad(-21.3); theta0 = deg2rad(-176.4);
a0 = sin((phi0-P.targetLatitude)/2)^2 + cos(phi0)*cos(P.targetLatitude) ...
    *sin((theta0-P.targetLongitude)/2)^2;
initialRange = 2*P.R*asin(sqrt(max(0,min(1,a0))));
Quantity = ["Altitude";"Longitude";"Latitude";"Speed"; ...
    "Flight-path angle";"Heading angle";"Distance to target"];
Unit = ["m";"deg";"deg";"m/s";"deg";"deg";"m"];
initialValues = [125000;-176.4;-21.3;4700;-10;-2.876;initialRange];
targetValues = [P.targetAltitude;rad2deg(P.targetLongitude); ...
    rad2deg(P.targetLatitude);P.targetSpeed;NaN;NaN;0];
resultValues = [altitude(end);rad2deg(longitude(end)); ...
    rad2deg(latitude(end));speed(end);rad2deg(xState(end,5)); ...
    rad2deg(xState(end,6));targetMiss];
errorValues = resultValues-targetValues;
errorValues(2) = mod(errorValues(2)+180,360)-180;
PaperInitial = string(compose('%.6f',initialValues));
Target = string(compose('%.6f',targetValues));
Result = string(compose('%.6f',resultValues));
Error = string(compose('%+.6f',errorValues));
Target(isnan(targetValues)) = "Not specified";
Error(isnan(errorValues)) = "Not applicable";
summary = table(Quantity,Unit,PaperInitial,Target,Result,Error);
fprintf('\nCONDITIONS AND RESULTS (error = result - target)\n');
disp(summary);

%% Path constraints: evaluate at logged samples, without enforcing them
Constraint = ["Aerodynamic load";"Dynamic pressure";"Heat flux"];
Unit = ["g";"kPa";"kW/m^2"];
Limit = [P.loadLimit;P.pressureLimit;P.heatLimit];
Peak = [peakLoad;peakPressure;peakHeat];
Margin = Limit-Peak;
PeakTime_s = tState([iLoad;iPressure;iHeat]);
Status = repmat("Within limit",3,1);
Status(Margin < 0) = "Limit exceeded";
constraintSummary = table(Constraint,Unit,Limit,Peak,Margin,PeakTime_s,Status);
fprintf('\nPATH CONSTRAINTS\n'); disp(constraintSummary);
fprintf('End time: %.3f s | Energy residual: %.3e J/kg | Last bank: %.3f deg\n', ...
    tState(end),energyResidual,bankMagnitude(end));

%% Save the run and both tables
if saveResults
    resultFolder = fullfile(projectFolder,'results', ...
        ['run_' datestr(now,'yyyymmdd_HHMMSS_FFF')]);
    assert(~isfolder(resultFolder),'Run folder already exists: %s',resultFolder);
    mkdir(resultFolder);
    save(fullfile(resultFolder,'results.mat'),'out','P','summary', ...
        'constraintSummary','tState','xState','tBank','sigmaBank', ...
        'loadG','pressureKPa','heatFlux','referenceBank','modelName');
    writetable(summary,fullfile(resultFolder,'summary.csv'));
    writetable(constraintSummary,fullfile(resultFolder,'path_constraints.csv'));
    savefig(fig,fullfile(resultFolder,'validation.fig'));
    exportgraphics(fig,fullfile(resultFolder,'validation.png'),'Resolution',200);
    fprintf('Saved to: %s\n',resultFolder);
end

%% Read vector or scalar timeseries: rows = time samples
function [t,values] = readHistory(signal)
    assert(isa(signal,'timeseries'),'Use Timeseries for To Workspace blocks.');
    t = signal.Time(:);
    assert(~isempty(t) && all(isfinite(t)) && all(diff(t) >= 0), ...
        'Logged times are empty, nonfinite or out of order.');
    if signal.IsTimeFirst
        values = reshape(signal.Data,numel(t),[]);
    else
        values = reshape(signal.Data,[],numel(t)).';
    end
    assert(all(isfinite(values(:))),'Logged data contains NaN or Inf.');
end
