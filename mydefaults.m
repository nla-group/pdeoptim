%% Setting up defaults.

close all hidden
clear all

randn('state',0);
rand('state',0);

height = 2/3;
set(0,...
'defaultfigureposition',[180 100 800 800*height],...
'defaultaxeslinewidth',1,...
'defaultaxesfontsize',16,...
'defaultlinelinewidth',2,...
'defaultpatchlinewidth',1,...
'defaultlinemarkersize',10,...
'defaulttextinterpreter','tex');

if exist(fullfile(cd,'chebfun'),'dir') == 7
    addpath(fullfile(cd,'chebfun'))
    savepath
else
    try 
        disp('Please wait: downloading Chebfun... (takes a few seconds but will be done only once!)')
        unzip('https://www.chebfun.org/download/chebfun_v4.2.2889.zip')
        movefile('chebfun_v4.2.2889/chebfun', 'chebfun')
        addpath(fullfile(cd,'chebfun')), savepath
        rmdir('chebfun_v4.2.2889')
        disp('Chebfun has been downloaded and added to your MATLAB path.')
    catch
        disp('Error: Please refer to README and download Chebfun v4.2.2889.')
    end
end
