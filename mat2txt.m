% transfer .mat data file to .txt file for the use of GIMME
clear all,clc,close all

% for i = 1:10
%    filename = strcat('./results/simu1_high2_new_',num2str(i),'.mat');
%    load(filename)
%    folder = strcat('./Datatxt/simu1_ins_',num2str(i));
%    mkdir(folder)
%    for j = 1:100
%        filename_out = strcat(folder,'/simu1_ins_',num2str(i),'_',num2str(j),'.txt');
%        dlmwrite(filename_out,Data{j}');
%    end
% end
   
% for i = 1:10
%    filename = strcat('./results/simu2_high2_new_',num2str(i),'.mat');
%    load(filename)
%    folder = strcat('./Datatxt/simu2_ins_',num2str(i));
%    mkdir(folder)
%    for j = 1:60
%        filename_out = strcat(folder,'/simu2_ins_',num2str(i),'_',num2str(j),'.txt');
%        dlmwrite(filename_out,Data{j}');
%    end
% end


% for i = 1:10
%    filename = strcat('./results/simu3_high2_new_',num2str(i),'.mat');
%    load(filename)
%    folder = strcat('./Datatxt/simu3_ins_',num2str(i));
%    mkdir(folder)
%    for j = 1:80
%        filename_out = strcat(folder,'/simu3_ins_',num2str(i),'_',num2str(j),'.txt');
%        dlmwrite(filename_out,Data{j}');
%    end
% end


% for i = 1:10
%    filename = strcat('./results/simu4_high2_new_',num2str(i),'.mat');
%    load(filename)
%    folder = strcat('./Datatxt/simu4_ins_',num2str(i));
%    mkdir(folder)
%    for j = 1:80
%        filename_out = strcat(folder,'/simu4_ins_',num2str(i),'_',num2str(j),'.txt');
%        dlmwrite(filename_out,Data{j}');
%    end
% end


% for i = 1:10
%    filename = strcat('./results/simu5_high2_new_',num2str(i),'.mat');
%    load(filename)
%    folder = strcat('./Datatxt/simu5_ins_',num2str(i));
%    mkdir(folder)
%    for j = 1:80
%        filename_out = strcat(folder,'/simu5_ins_',num2str(i),'_',num2str(j),'.txt');
%        dlmwrite(filename_out,Data{j}');
%    end
% end


% for i = 1:10
%    filename = strcat('./results/simu13_high2_new_',num2str(i),'.mat');
%    load(filename)
%    folder = strcat('./Datatxt/simu13_ins_',num2str(i));
%    mkdir(folder)
%    for j = 1:100
%        filename_out = strcat(folder,'/simu13_ins_',num2str(i),'_',num2str(j),'.txt');
%        dlmwrite(filename_out,Data{j}');
%    end
% end

% for i = 1:10
%    filename = strcat('./results/simu14_high2_new_',num2str(i),'.mat');
%    load(filename)
%    folder = strcat('./Datatxt/simu14_ins_',num2str(i));
%    mkdir(folder)
%    for j = 1:80
%        filename_out = strcat(folder,'/simu14_ins_',num2str(i),'_',num2str(j),'.txt');
%        dlmwrite(filename_out,Data{j}');
%    end
% end

for i = 1:10
   filename = strcat('./results/simu15_high2_new_',num2str(i),'.mat');
   load(filename)
   folder = strcat('./Datatxt/simu15_ins_',num2str(i));
   mkdir(folder)
   for j = 1:60
       filename_out = strcat(folder,'/simu15_ins_',num2str(i),'_',num2str(j),'.txt');
       dlmwrite(filename_out,Data{j}');
   end
end

% for i = 1:10
%    filename = strcat('./results/simu16_high2_new_',num2str(i),'.mat');
%    load(filename)
%    folder = strcat('./Datatxt/simu16_ins_',num2str(i));
%    mkdir(folder)
%    for j = 1:100
%        filename_out = strcat(folder,'/simu16_ins_',num2str(i),'_',num2str(j),'.txt');
%        dlmwrite(filename_out,Data{j}');
%    end
% end

% for i = 1:10
%    filename = strcat('./results/simu17_high2_new_',num2str(i),'.mat');
%    load(filename)
%    folder = strcat('./Datatxt/simu17_ins_',num2str(i));
%    mkdir(folder)
%    for j = 1:100
%        filename_out = strcat(folder,'/simu17_ins_',num2str(i),'_',num2str(j),'.txt');
%        dlmwrite(filename_out,Data{j}');
%    end
% end