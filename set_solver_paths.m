if ismac
    addpath(genpath([pwd,'/solver/mac']));
elseif isunix
    addpath(genpath([pwd,'/solver/linux']));
elseif ispc
    addpath(genpath([pwd,'/solver/win64']));
else
    disp('Platform not supported')
end