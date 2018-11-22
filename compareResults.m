%% compareResults
% This script processes reads the results of an OCT experiment, as given by
% the script processExperiment, and writes a comparison of all statistics
% in excell format.
%
% Consecutive sheets in the excel file contain the results of consecutive
% patients.
%
%% Copyright and License Notices
% This file is part of
% "Assessing the Effect of Contact Lenses on the Image Quality of Retinal
% Spectral Domain Optical Coherence Tomography Using Automated Image
% Analysis" by Yinon Shapira, Talia Aviram, Omer Granak, Igor Viner,
% Erez Ribak, Eitan Z Blumenthal (2018).
%
% Copyright © 2017, 2018 Omer Granek
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% Code
close all; clearvars; clc;
analysisDir = 'ValidationAnalysis/withCut'; % set to desired dir for in/out
main = cd(analysisDir);
mes = dir('**\*.mat');
numOfPatients = 20;

B = 1:numOfPatients;
for b = B
    means=[];
    stds=[];
    SNRs=[];
    CNRs=[];
    ENLs=[];
    DRs=[];
    disp(num2str(b));
    for c = 1:1
        load(mes(b).name, 'stats');
        means = [means, stats.means];
        stds = [stds, stats.stds];
        SNRs = [SNRs, [squeeze(stats.SNRs(1,:,:));...
            squeeze(stats.SNRs(2,:,:))]];
        CNRs = [CNRs, [squeeze(stats.CNRs(1,:,:));...
            squeeze(stats.CNRs(2,:,:))]];
        ENLs = [ENLs, stats.ENLs];
        DRs = [DRs, stats.DRs];
        xlswrite('Mean',means,b);
        xlswrite('Std',stds,b);
        xlswrite('SNR',SNRs,b);
        xlswrite('CNR',CNRs,b);
        xlswrite('ENL',ENLs,b);
        xlswrite('DR',DRs,b);
    end
end