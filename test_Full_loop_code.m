% This scripts tests the Full_loop_func libray with already measured data.

% clear all;

%filename = '/home/ic314/Documents/MATLAB/cOMPRESSED_sENSING_in_se2pc_instrument/dy014926.mat'; % simple decay
%filename = '/home/ic314/Documents/MATLAB/cOMPRESSED_sENSING_in_se2pc_instrument/dy011880.mat'; % phonon measurement & simmetric ibase range
%filename = '/home/shared/SpinEcho/Spinecho2Icon/data/Cu111_4/dy014997.mat'; % phonon measurement & simmetric ibase range
% filename = '/home/shared/SpinEcho/Spinecho2Icon/data/Pd111/dy014034.mat'; % phonon measurement & simmetric ibase range
%filename = 'C:\users\navidor\na364\data\Pd111\dy014036.mat';
%filename = '/home/shared/SpinEcho/Spinecho2Icon/data/Cu111_4/dy015600.mat'; % phonon measurement & simmetric ibase range
filename = 'test_file.mat';

load(filename)

% to test single decay (ibase from 0 A to 10 A) use the following distribution probabilities: 0.1*[10,4,2,2,2,2,2]
  
% to test simmetric ibase range use the following distribution probabilities: 0.1*[2,2,4,10,4,2,2]

% generate ibase_matrix and CS_s structure based on the histogram approach:
[ibase_matrix,CS_s] = Distrib_func.gen_ibase_matrix('Imin',min(meas.ibase),'Imax',max(meas.ibase),'deltaI',min(abs(diff(meas.ibase))),...
     'numloops',meas.numloops,'histogramstr','on','levels',[1/8,1/4,3/8,5/8,3/4,7/8,1],'prob_levels',0.1*[10,10,10,10,10,10,10],'Method','DFT','file_test',filename);

% generate ibase_matrix and CS_s structure based on the continuous probability distribution function approach:
%[ibase_matrix,CS_s] = Distrib_func.gen_ibase_matrix('Imin',min(meas_test.ibase),'Imax',max(meas_test.ibase),'deltaI',(max(meas_test.ibase)-min(meas_test.ibase))/length(meas_test.ibase),...
%    'numloops',meas_test.numloops,'pdf_contstr','on','num_levels',[1 1 1 2 3 4 4 5 6 7],'FWHM',0.5,'fractp_loop',0.1,'Method','DWT','file_test',filename);

% generate ibase_matrix and CS_s structure based on the discrete probability distribution function approach:
%[ibase_matrix,CS_s] = Distrib_func.gen_ibase_matrix('Imin',min(meas_test.ibase),'Imax',max(meas_test.ibase),'deltaI',(max(meas_test.ibase)-min(meas_test.ibase))/length(meas_test.ibase),...
%    'numloops',meas_test.numloops,'pdf_discrstr','on','num_levels',[1 1 1 2 3 4 4 5 6 7],'FWHM',0.5,'fractp_loop',0.35,'Method','DFT','file_test',filename);



% do experiment:
meas = Full_loop_func.CS_scan_current('psu',meas.psu,'ibase_matrix',ibase_matrix, 'numloops',meas.numloops,'savemode',0,'CS_s',CS_s,'file_test',filename);
  