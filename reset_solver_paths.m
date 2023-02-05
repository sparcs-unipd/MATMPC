if ismac
    rmpath(genpath([pwd,'/solver/mac']));
elseif isunix
    rmpath(genpath([pwd,'/solver/linux']));
elseif ispc
    rmpath(genpath([pwd,'/solver/win64']));
else
    disp('Platform not supported')
end
