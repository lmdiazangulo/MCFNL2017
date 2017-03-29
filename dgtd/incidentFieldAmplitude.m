function [res] = incidentFieldAmplitude(As, freq, time, ...
                                        Ag, spread, delay)
% As     := Sine amplitude.
% Freq   := frequency in Hertzs
% Time   := Time in seconds.
% Ag     := Gaussian amplitude.
% spread := Sigma coefficient in gaussian.
% delay  := mu coefficient in gaussian.

expArg = (time - delay) / (spread * sqrt(2));
res = As .* cos(2 .* pi .* freq .* time) .* Ag .* exp( - expArg.^2);


% To obtain a frequency bandwidth fm the spread must be of the order of:
% spread = sqrt(1 ./ (2 * pi^2 * log(sqrt(2)))) ./ fm;

% Coarse calculations:
% To get a pulse with a bandwidth of 1200 THz -> pulse of spread = 30 fs 