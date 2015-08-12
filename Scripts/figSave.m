function [] = figSave(h,name,desc,e,s,t)
figure(h)

path = [cd, '/StripeFigures/Model/'];
filename = [name, '_E', num2str(e), '_C', num2str(desc),...
    '_S',num2str(s),'_T', num2str(t), '.png'];

String = [path,filename];

export_fig(String)

close all