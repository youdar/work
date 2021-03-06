% Clashcompare compares clash scores produced internaly in
% PDB_interpretation with the scores produced by PROBE
clc
% read data
% load('clashscore_compare_and_name.txt');
% x = clashscore_compare_and_name(:,2); 
% y = clashscore_compare_and_name(:,3);
load('clashscore_compare_and_name_no_sym_op_same_Ovdw.txt');
x = clashscore_compare_and_name_no_sym_op_same_Ovdw(:,2); 
y = clashscore_compare_and_name_no_sym_op_same_Ovdw(:,3);
time_RM = clashscore_compare_and_name_no_sym_op_same_Ovdw(:,4);
time_Probe = clashscore_compare_and_name_no_sym_op_same_Ovdw(:,5);
% fitting function

[fittedFun,goodness_of_fit]=fit(x,y,'poly1');
% disp(fittedFun)
% prducing y2 for the fitted function
subplot(1,2,1);
y2 = fittedFun(x);
plot(x,y,'.',x,y2);
xlabel('clashscore - internal');
ylabel('clashscore - probe');
% title('Compare internal clash-score to PROBE clashscore');
title('Compare clashscores : pdb_iterpretation to PROBE ');
xlim([0,50]);
ylim([0,50]);
legend('clashscore - set 1','fitted function');
% plot run time
x = [0,150];
y = [0,150];
subplot(1,2,2);
plot(time_RM,time_Probe,'.',x,y);
xlim([0,100]);
ylim([0,100]);
xlabel('restriant manager - time');
ylabel('PROBE - time');
text(5,95,'The added time by the');
text(5,90,'new clash score is');
text(5,85,'in the range of 0.2 to 2.5sec');

