% using the rotation matrices from the
% c:\Phenix\Dev\Work\work\FAB\FAB_elbow_angle.py
% calculate and draw the axis of each rotation
close
clc

% a = [0.543792658203, 0.404159121759, 0.735489598284; 0.329232559623, -0.908879202461, 0.256016634261; 0.771942657565, 0.10292715693, -0.627304179643];
% b = [-0.944889156632, -0.0645463882279, 0.320964554813; -0.0149904908774, -0.970814509946, -0.239362637982; 0.327047040749, -0.230982577363, 0.916344521505];

% 1dba
a = [-0.944856, 0.028753, 0.326222; 0.022132,-0.988254, 0.151207; 0.326738, 0.150089, 0.933122] % v - ours
b = [-0.918866, 0.053961, 0.390862;-0.000921,-0.990895, 0.134634; 0.394569, 0.123351, 0.910549] % V - site


c = [-0.983924,-0.086648, 0.156161; 0.133151,-0.938654, 0.318119; 0.119017, 0.333798, 0.935101] % C - ours
d = [-0.938461,-0.163679, 0.304139; 0.186099,-0.981452, 0.046043; 0.290962, 0.099810, 0.951514] % C - site

% d = [1,0,0; 0,1,0;0,0, 1] % C - site

% 1bbd
% a = [-0.824573, 0.083149, 0.559612;-0.142942, -0.987668, -0.063871; 0.547400, -0.132658, 0.826290]
% b = [-0.932296  0.074772  0.353882;-0.132476, -0.981001, -0.141729; 0.336562, -0.179014  0.924489]
% 
% 
% c = [ 0.434167, 0.400764, 0.806776; 0.427656, -0.879930, 0.206960; 0.792848, 0.255168, -0.553426]
% d = [ 0.580606  0.097505  0.808325; 0.231973, -0.971465 -0.049439; 0.780439  0.216214, -0.586657]

% get eigen vectors and values
[Va,Da] = eig(a)
[Vb,Db] = eig(b)
[Vc,Dc] = eig(c)
[Vd,Dd] = eig(d)

% Get the eigen vector for the eigen value 1
eigen_v_a = Va(:,1)
eigen_v_b = Vb(:,1)
eigen_v_c = Vc(:,1)
eigen_v_d = Vd(:,1)

% Plot
% o = [0;0;0];
% la = [o,eigen_v_a];
% lb = [o,eigen_v_b];
% figure()
% plot3(la(1,:),la(2,:),la(3,:),'r',lb(1,:),lb(2,:),lb(3,:),'b')

% hold on
% plot3(la(1,:),la(2,:),la(3,:),'r')
% plot3(lb(1,:),lb(2,:),lb(3,:),'b')
% hold off


c = dot(eigen_v_a,eigen_v_b);
disp(['var_ours dot var_site  product is: ',num2str(c)])
disp('   ')
ang = 180*acos(c)/pi;
if ang<90
    ang = 180 - ang;
end;
disp(['The angle between them is: ',num2str(ang)])
disp('   ')
c = dot(eigen_v_c,eigen_v_d);
disp(['const_ours dot const_site   product is: ',num2str(c)])
ang = 180*acos(c)/pi;
if ang<90
    ang = 180 - ang;
end;
disp('   ')
disp(['The angle between them is: ',num2str(ang)])

dot(cross(eigen_v_a,eigen_v_c),cross(eigen_v_b,eigen_v_d))

x = cross(Vb(:,1),Vb(:,2))
sqrt(dot(x,x))
Vb(:,3)
% cross(Vb(:,1),Vb(:,2))
% Vb(:,3)
% cross(Vc(:,1),Vc(:,2))
% Vc(:,3)
% cross(Vd(:,1),Vd(:,2))
% Vd(:,3)


