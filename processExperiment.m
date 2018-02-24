%% processExperiment
% This script processes an entire OCT based experiment including multiple
% participants and configurations. It makes use of the functions evaluateB
% and analyzeStats, given seperately.
%
% This program assumes image files exported from Heidelberg Engineering
% software, and requires each image to be given both in segmented and
% unsegmented form.
%
% The script assumes that within each input directory, there are
% numOfParticipants directories labeled as 1,2...,numOfParticipants.
% Within each of the latter there should be two directories, one for each
% configuration labled in '####i' format, where i=1,2. Within each of the
% latter, there should be a constant number of OCT B-Scan images.
%% Copyright and License Notices
% This file is part of
% "Assessing the Effect of Contact Lenses on the Image Quality of Retinal
% Spectral Domain Optical Coherence Tomography Using Automated Image
% Analysis" by Yinon Shapira, Talia Aviram, Omer Granak, Igor Viner,
% Erez Ribak, Eitan Z Blumenthal (2018).
%
% Copyright © 2017, 2018 Omer Granek and Erez Ribak
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

% preliminaries
close all; clearvars; clc;
analysisDir = 'Analysis\With Vessels'; % set to desired dir for output
octDir = 'OCT';
 % set to dir containing images with segmentation
octNoSegDir = [octDir '\without segmentation'];
 % set to dir containing images without sigmentation
octSegDir = [octDir '\with segmentation'];
saveFigures = false;
numOfParticipants = 20;

for b = 1:numOfParticipants
    %read from directory with and without segmentation
    zdSub = num2str(b);
    zd = [octNoSegDir '\' zdSub];
    zdSeg = [octSegDir '\' zdSub];
    if isdir(zd)
        dd = dir(zd);
        ddSeg = dir(zdSeg);
        for c = 1:2 % without / with contact
            pathnameSub = dd(c+2).name;
            pathnameSegSub = ddSeg(c+2).name;
            fsaveSub = [zdSub regexprep(pathnameSub,{'\.',' '},'')];
            pathname = [zd '\' pathnameSub]; disp(pathname)
            pathnameSeg = [zdSeg '\' pathnameSegSub];
            
            % cuts
            dirinfo = dir (pathname); ld = length(dirinfo); ln = 3:ld;
            fsh1 = 0; fsh2 = 0; fsv1 = 0; fsv2 = 0;
            means = zeros(ld - 2, 4);
            stds = zeros(ld - 2, 4);
            DRs = zeros(ld - 2, 1);
            
            zz = cell(1,2);
            for j = 1:2
                fname = [pathname '\' dirinfo(floor((ld-2) * j / 3)).name];
                zz{j} = imread(fname);
            end
            if saveFigures, figs = gobjects(1,ld - 2); end
            
            for l = ln
                %includes face-on and cut
                fname = [pathname '\' dirinfo(l).name];
                %with segmentation
                fnameSeg = [pathnameSeg '\' dirinfo(l).name];
                % skips missing files for standard dev. comp.
                segExists = (exist(fnameSeg, 'file') == 2);
                zz = imread (fname);
                if segExists, zzSeg = imread(fnameSeg); end
                z = double(zz (1:420, 497:1008, 1)); % keep only 2d cut
                
                % horizontal power spectrum
                fth = abs (fft (z)); fth2 = fth .* fth; 
                % add to previous power spectra
                fsh1 = fsh1 + fth; fsh2 = fsh2 + fth2;
                % vertical power spectrum
                ftv = abs (fft (z, 768, 2)); ftv2 = ftv .* ftv;
                % add to previous power spectra
                fsv1 = fsv1 + ftv; fsv2 = fsv2 + ftv2; 
                
                %standard dev. of 4 segments
                if segExists
                    if saveFigures
                        [means(l-2, :), stds(l-2, :), DRs(l-2), ...
                            figs(l-2)] = evaluateB(zz, zzSeg, 420);
                    else
                        [means(l-2, :), stds(l-2, :), DRs(l-2), ~] = ...
                            evaluateB(zz, zzSeg, 420);
                    end
                end
            end
            
            % obtain further stats
            stats = analyzeStats(means, stds, DRs, ld - 2);
            if (b < 10)
                saveNulls = '\00';
            else
                saveNulls = '\0';
            end
            fsave = [analysisDir saveNulls fsaveSub];
            
            % save average power spectra of cuts and standard devs.,
            % change directory name to Fouriers
            if segExists
                save([fsave '.mat'], 'fsh1', 'fsh2', 'fsv1', 'fsv2', 'stats');
                if saveFigures,
                    savefig(figs, [fsave '.fig']);
                end
            else
                save(fsave, 'fsh1', 'fsh2', 'fsv1', 'fsv2');
            end
        end
    end
end