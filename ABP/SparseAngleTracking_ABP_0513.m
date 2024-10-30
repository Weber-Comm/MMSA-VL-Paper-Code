%% LAST MODIFIED BY Weber, 2023.5.13 %%
% Author: Weibo Wen
% Contact: iweber.wen@gmail.com
% Created: March 24, 2023

% This code performs auxiliary beam pair algorithm for beam tracking.

% D. Zhu, J. Choi, Q. Cheng, W. Xiao, and R. Heath,
% “High-Resolution Angle Tracking for Mobile Wideband Millimeter-Wave Systems
% With Antenna Array Calibration,” IEEE Trans. Wireless Commun.,
% vol. 17, no. 11, pp. 7173–7189, Nov. 2018, doi: 10.1109/TWC.2018.2865759.

% Rest of the code goes here
%% INITIALIZATION %%
clear;
rng(1);
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

fprintf("===== New Run =====\n");
%% PARAMETER SETTING %%

Nt = 16; % number of transmit antennas
Nr = 16; %
antenna_spacing = 0.500; %  d/wavelength

% all varibles related to angle or phase are in radians
% spatial frequencies, not actual angle

% _SPATIAL_FREQUENCY_ = 2*pi*antenna_spacing*sin(_ACTUAL_ANGLE_);

% set parameter "delta" for a shift from boresight angle
delta = pi/16 + eps; % 13

% set parameter "eta_original" for perfect beam training, initially
eta = pi + eps;

% !! THE VARIBLE RELATED TO ANGLE OR PHASE MAY EXCEED BOUNDRIES (IS COMPLEX)
% !! DUE TO ABP TRACKING IN COMPLEX MULTIPATH ENVIRONMENT

%% LOAD DATA FROM SCENE SPARSE
% one CSI struct
name_mat_csi = "Tx2_Rx6_38GHz_sparse_haotian.mat";
s = load(name_mat_csi);
fieldname_mat_csi = "CIR_" + strrep(name_mat_csi,'_38GHz_sparse_haotian.mat',''); %
eval('csi_seq = s.' + fieldname_mat_csi + ';');
max_snapshot = length(csi_seq);

% "csi_seq"

% ULA center position (ground truth)
M = readmatrix('sparse_groundtruth.xlsx','Sheet','AntPos','Range','B3:P1502');
ant2_xy = M(:,1:2);
ant5_xy = M(:,5:6);
ant6_xy = M(:,9:10);
ant8_xy = M(:,13:14);

ant_2to6_xy = ant6_xy - ant2_xy;
angle_2to6_yaxis_original = -atan2(ant_2to6_xy(:,2),ant_2to6_xy(:,1));

ant_2to6_xy(1200:1400,1) = ant_2to6_xy(1200:1400,1) + 20;

% for Tx2 to Rx6, angle from the positive direction of y-axis
angle_2to6_yaxis = -atan2(ant_2to6_xy(:,2),ant_2to6_xy(:,1));

%% PARAMETER CONVERSION DISPLAY %%

delta_real_angle = asin(delta/(2*pi*antenna_spacing));
eta_real_angle = asin(eta/(2*pi*antenna_spacing));
fprintf("Real Angle DELTA = %f in rad, %f in deg\n",...
    delta_real_angle,rad2deg(delta_real_angle));
fprintf("Real Angle ETA = %f in rad, %f in deg\n",...
    eta_real_angle,rad2deg(eta_real_angle));

spatial_freq_zc0 = eta - delta + eps;
zc0_real_angle = asin(spatial_freq_zc0/(2*pi*antenna_spacing));
spatial_freq_zc1 = eta + delta + eps;
zc1_real_angle = asin(spatial_freq_zc1/(2*pi*antenna_spacing));

fprintf("Real Angle ZC0 = %f in rad, %f in deg\n",...
    zc0_real_angle,rad2deg(zc0_real_angle));
fprintf("Real Angle ZC1 = %f in rad, %f in deg\n",...
    zc1_real_angle,rad2deg(zc1_real_angle));

%% CALCULATION ZC %%

% subcarrier index
k = (1:512)';

% ZC sequence index
i0 = 1;
i1 = 2;

% ZC symbols repetetions in time domain
sym_num = 1;

% complex baseband ZC symbols
s_zc0 = exp(-1j*(pi*(k-1).*k*i0/max(k)));
s_zc1 = exp(-1j*(pi*(k-1).*k*i1/max(k)));

% symbols (frequency domain) => time domain samples
x_zc0 = ifft(s_zc0);
x_zc1 = ifft(s_zc1);

% ZC symbols repetetions in time domain
x_zc0 = repmat(x_zc0,1,sym_num);
x_zc1 = repmat(x_zc1,1,sym_num);

%% LOOP TRACKING %%


csi_eng_seq = [];
file_index = 1:max_snapshot;

AoD_rad_estimated_seq = zeros(max_snapshot,1);

for i = file_index(1:1500)

if mod(i-1,1) ~= 0
    AoD_rad_estimated_seq(i) = AoD_rad_estimated_seq(i-1);
    continue;
end


%disp(i);

% VIWI LOAD CSI
% load([file_list(i).folder '\' file_list(i).name]);
% H = data.channel((1:Nt)+32,(1:Nr)+32);

% CROSSING LOAD CSI
H = csi_seq{i}((1:Nt)+0,(1:Nr)+0) + eps;
csi_eng = csi_seq{i}(:,:);
csi_eng = sum(abs(csi_eng).^2,'all');
csi_eng_seq = [csi_eng_seq csi_eng];

spatial_freq_zc0 = eta - delta + eps;
spatial_freq_zc1 = eta + delta + eps;

spatial_freq_r = -eta; % receive beamforming shifted phase

% beamformer
beamformer_zc0 = exp(1j*(0:spatial_freq_zc0:spatial_freq_zc0*(Nt-1))');
beamformer_zc1 = exp(1j*(0:spatial_freq_zc1:spatial_freq_zc1*(Nt-1))');

% spatial_freq_AoD = 2*pi*antenna_spacing*sin(AoD_rad_real); % real
% steering_vector_AoD = exp(1j*(0:spatial_freq_AoD:spatial_freq_AoD*(Nt-1))');
% steering_vector_AoA = exp(1j*(0:spatial_freq_AoA:spatial_freq_AoA*(Nr-1))');

% receive beamforming vector \vartheta
beamformer_r = exp(1j*(0:spatial_freq_r:spatial_freq_r*(Nr-1))');

% y_s_zc0 = beamformer_r'*H'*beamformer_zc0*s_zc0.';
% y_s_zc1 = beamformer_r'*H'*beamformer_zc1*s_zc1.';
y_s = beamformer_r'*H'*[beamformer_zc0 beamformer_zc1]*[s_zc0 s_zc1].';
y_s = y_s.';

% correlating
Lambda_Delta_az = s_zc0'*y_s;
Lambda_Sigma_az = s_zc1'*y_s;

% received signal strength
chi_Delta_az = Lambda_Delta_az'*Lambda_Delta_az;
chi_Sigma_az = Lambda_Sigma_az'*Lambda_Sigma_az;

% ratio metric
zeta_az = (chi_Delta_az-chi_Sigma_az)/(chi_Delta_az+chi_Sigma_az);
psi_esimated = eta - ...
    asin((zeta_az*sin(delta) - zeta_az*sqrt(1-zeta_az^2)*sin(delta)*cos(delta))/...
    (sin(delta)^2 + zeta_az^2*cos(delta)^2));

% disp(zeta_az);

spatial_freq_AoD_estimated = psi_esimated;

% polynomial fitting by curves in codebook
% if UES_CODEBOOK == 1
%     spatial_freq_AoD_estimated = polyval(p,zeta_az) + eta;
% end

AoD_rad_estimated = asin(spatial_freq_AoD_estimated/(2*pi*antenna_spacing));
AoD_rad_estimated_seq(i) = AoD_rad_estimated;

% update AoD
if (-pi<=eta)&&(eta<=pi)
    eta = spatial_freq_AoD_estimated;
elseif eta>pi
    eta = pi;
elseif eta<-pi
    eta = -pi;
    %fprintf("!!!!!ETA EXCEEDING BOUNDARY!!!!!\n");
    %break;
end


end

%% PLOT AND ERROR %%

f1 = figure;

p1 = plot(real(AoD_rad_estimated_seq)-angle_2to6_yaxis_original+angle_2to6_yaxis); %,'LineStyle','--'
p1.LineWidth = 0.5;
p1.LineStyle = "-";
p1.Color = "blue";
p1.Color(4) = 0.8;

hold on;

p2 = plot(angle_2to6_yaxis,"Color", [0.4, 0.4, 0.4, 0.5]);
p2.LineWidth = 1;
p2.LineStyle = "--";
p2.Color = "black";
p2.Color(4) = 0.8;

set(gcf,'Position',[400 400 400 300])

ylim([-pi/2,pi/2]);
legend('Estimated AoD','Real AoD','Location','southwest');
xlabel('Snapshot Index','FontSize',14);
ylabel('AoD, Arrival of Departure (rad)','FontSize',14);
title(' ','FontSize',14);
grid on;

%% ERROR CALCULATION
error_seq = real(AoD_rad_estimated_seq)-angle_2to6_yaxis_original;

mse_error = sum((error_seq.^2))/max_snapshot;
mean_error = sum(abs(error_seq))/max_snapshot;


figure;
plot(error_seq,'Color','red');
hold on;
yline(0,"Color",[0.15 0.15 0.15 0.1],'LineWidth',1);

set(gcf,'Position',[400 400 400 300])

ylim([-0.3,0.3]);
xlabel('Snapshot Index','FontSize',14);
ylabel('Angle Error (rad)','FontSize',14);
grid on;

fprintf("MSE ERROR = %f rad\n",mse_error);
fprintf("MEAN ERROR = %f rad\n, ",mean_error);



