function posMatrix_t=spline3dpoints(posMatrix)

YY=posMatrix(:,2:4);
poslist=posMatrix(:,1);

num_cluster=size(YY,1);
F=spline((1:num_cluster),YY');
% Trajectory
step=0.05;
t=[1:step:num_cluster];
YYt=ppval(F,t)';

F=spline((1:num_cluster),poslist);
poslist_t=ppval(F,t);
posMatrix_t=[poslist_t' YYt];
end