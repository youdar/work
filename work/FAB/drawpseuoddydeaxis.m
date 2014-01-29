% using the rotation matrices from the
% c:\Phenix\Dev\Work\work\FAB\FAB_elbow_angle.py
% calculate and draw the axis of each rotation
close
clc

a = [0.543792658203, 0.404159121759, 0.735489598284; 0.329232559623, -0.908879202461, 0.256016634261; 0.771942657565, 0.10292715693, -0.627304179643];
b = [-0.944889156632, -0.0645463882279, 0.320964554813; -0.0149904908774, -0.970814509946, -0.239362637982; 0.327047040749, -0.230982577363, 0.916344521505];

% get eigen vectors and values
[Va,Da] = eig(a);
[Vb,Db] = eig(b);

% Get the eigen vector for the eigen value 1
eigen_v_a = Va(:,1)
eigen_v_b = Vb(:,1)

% Plot
o = [0;0;0];
la = [o,eigen_v_a];
lb = [o,eigen_v_b];
figure()
plot3(la(1,:),la(2,:),la(3,:),'r',lb(1,:),lb(2,:),lb(3,:),'b')

% hold on
% plot3(la(1,:),la(2,:),la(3,:),'r')
% plot3(lb(1,:),lb(2,:),lb(3,:),'b')
% hold off


c = dot(eigen_v_a,eigen_v_b);
disp(['The dot product is: ',num2str(c)])
ang = 180*acos(c)/pi;
if ang<90
    ang = 180 - ang;
end;
disp(['The angle between them is: ',num2str(ang)])