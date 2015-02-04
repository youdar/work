% Clashcompare compares clash scores produced internaly in
% PDB_interpretation with the scores produced by PROBE
clc
% read data
load('clashscore_compare_reduce_12_6_2013.txt');
raw_data = clashscore_compare_reduce_12_6_2013;
clear clashscore_compare_reduce_12_6_2013;

% file_name,total_nb_clashscore,without_sym_nb_clashscore,clashscore_probe,
% total_nb_clashscore_time,clashscore_probe_time
RM_clashscore = raw_data(:,2); 
RM_clashscore_no_sym = raw_data(:,3); 
clashscore = raw_data(:,4);
time_RM = raw_data(:,5);
time_Probe = raw_data(:,6);
% fitting function

[fittedFun,goodness_of_fit]=fit(RM_clashscore_no_sym,clashscore,'poly1');
% disp(fittedFun)
% prducing y2 for the fitted function
y2 = fittedFun(RM_clashscore_no_sym);
%
[fittedFun,goodness_of_fit]=fit(RM_clashscore,clashscore,'poly1');
% disp(fittedFun)
% prducing y2 for the fitted function
y3 = fittedFun(RM_clashscore);
%
subplot(2,2,1);
plot(RM_clashscore_no_sym,clashscore,'.',RM_clashscore_no_sym,y2);
xlim([0,50]);
ylim([0,50]);
xlabel('clashscore - RM');
ylabel('clashscore - probe');
% title('Compare internal clash-score to PROBE clashscore');
title('Restraints manager vs PROBE - No clahses due to sym. op.');
legend('clashscore','fitted function');
%
subplot(2,2,2);
plot(RM_clashscore,clashscore,'.',RM_clashscore,y3)
xlim([0,50]);
ylim([0,50]);
title('Restraints manager vs PROBE');
legend('clashscore','fitted function');
xlabel('clashscore - RM');
ylabel('clashscore - probe');
% plot run time
x = [0,150];
y = [0,150];
subplot(2,2,3);
plot(time_RM,time_Probe,'.');
xlim([0,100]);
ylim([0,100]);
xlabel('restriant manager - time');
ylabel('PROBE - time');
%
% text(5,95,'The added time by the');
% text(5,90,'new clash score is');
% text(5,85,'Note that the ');
text(130,59,'Note that the oxygen vdw radii are different');
text(130,47,'between Probe and Resitrans manager');
text(130,85,'Hydrogens added using phenix.ready\_set');
