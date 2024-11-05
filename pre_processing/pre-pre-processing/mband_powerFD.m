function fd = mband_powerFD(mov, TR, head)
% This function computes framewise displacement according to the Power
% method. This code is edited for multiband data according to Power et al.
% 2019. It filters out respiratory frequencies, and calculates fd over
% longer delays. Movement parameters can be taken from mcflirts .par output

% This script is adapted from GetFDPower.m (Linden Parkes, Alex Fornito.
% Brain & Mental Health Laboratory, 2016) and GetFDJenk_multiband.m (Linden
% Parkes, Tribikram Thapa, Ben Fulcher & Alex Fornito, 2020).
% According to the filter described in Fair et al., (2020)

% Updated script:
% Kane Pavlovich, & Alex Fornito, 2022.

% INPUTS
% ------
% mov - an N x 6 matrix containing 6 movement parameters and where
% N = length of the time series.
% TR - TR of scan in seconds
% head - head radius (in mm). Default = 50mm
%
% -------
% OUTPUTS
% -------
% fd - N-1 length vector representing the total framewise
% displacement. There will be k zeros at the start of the
% vector.

if nargin < 3
    head=50;
 

% convert degrees to radians to mm
mov(:,4:6) = head*pi/180*mov(:,4:6);

% % Filter out respiratory frequencies based on Fair et al. 2020
stopband=[0.31 0.41];
nyq = (1/TR)/2; % nyquist
Wn = stopband/nyq; % cutoff frequencies
order = 10; % controls steepness of filter roll-off. Set to same value as in Power et al.2019
[B, A] = butter(order,Wn,'stop'); % construct filter

% LPFilter based on Williams et al. 2022
%pb = 0.2./((1/TR)/2); %0.2 Hz LOW PASS FILTER
%[B,A] = butter(2,pb,'low'); %Make butterworth filter

% detrend and filter motion parameters
mov = detrend(mov,'linear'); % HCP already detrended so no need.
mov = filtfilt(B, A, mov);

% differentiate movement parameters
delta_mov = [
            zeros(1,size(mov,2));
		    diff(mov);
			];

% compute total framewise displacement
fd = sum(abs(delta_mov),2);

end
