function serial = define_datafile_mfile(prefix, isAlsoFolder)

% Function to redefine and open the data file, in which to store results

global data_path
global data_currentRun_path
global data_filename
global data_file

% Look to see if there is a serial number stored to disk, create if necessary

home_path=CustomFuncs.find_home_path();

number_file = [home_path 'status/' prefix '.num'];

if ~exist([home_path 'status/'],'dir')
    home_path = pwd;
    number_file = [home_path '/' prefix '.num'];
end

fid = fopen (number_file,'r');
if (fid == -1)
  serial = 1;
else
  newnum = fscanf(fid,'%s');
  fclose(fid);
  serial = str2num(newnum);
end
fid2 = fopen(number_file,'w');
fprintf(fid2,'%d',serial+1);
fclose(fid2);

% make the filename for the data file
if strcmp(data_path,'') || isempty(data_path)
    data_path = uigetdir(home_path ,'please define a data_path');
end
    
serstr = sprintf('%.6d',serial);
if exist('isAlsoFolder','var') && isAlsoFolder == 1
    data_currentRun_path = [data_path '/' prefix serstr];
    if ~exist(data_currentRun_path,'dir'), mkdir(data_currentRun_path); end
else
    data_currentRun_path = data_path;
end
data_filename = [prefix serstr '.mat'];
data_file = [data_currentRun_path '/' data_filename];

% open the datafile (append data to the end) and store the time

filetime.creation=fix(clock);
save(data_file,'filetime');
