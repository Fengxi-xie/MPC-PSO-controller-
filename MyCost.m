function [L, sol]=MyCost(Position, model)
%解码过程
Position = reshape(Position, [], 2);
%%      ------------------参数定义与矩阵声明---------------------------%%
%%  --------模型参数-------%%
psi_ref = 0;%初始转向角 车身横摆角
delta_ref = 0;%初始前轮偏转角
R = 100;%圆半径

%%  --------控制参数-------%%
T = model.T;%采样时间
Np=model.Np;% 预测步长
Nc=model.Nc;% 控制步长

%%  --------矩阵定义-------%%
%状态空间方程x(k+1) = Ak*x(k)+Bk*u;
Ak = [1 0 -model.v_ref*sin(psi_ref)*T; 0 1 model.v_ref*cos(psi_ref)*T;0 0 1];
Bk = [cos(psi_ref)*T 0; sin(psi_ref)*T 0; tan(delta_ref)*T/model.l model.v_ref*T/(model.l*cos(delta_ref))];

% 优化目标参数，加权矩阵
Q=eye(3); R=eye(2);

% 转化为用控制量ut表示的，关于状态量的推导方程的矩阵
At=[]; Bt=[]; temp=[];

% 转换后的加权矩阵
Qt=[]; Rt=[];

%%      ------------------计算---------------------------%%
% 加权矩阵的计算过程，以及推导方程矩阵的叠加过程
for i=1:Nc
    At=[At; Ak^i];
    Bt=[Bt zeros(size(Bt,1), size(Bk,2));
        Ak^(i-1)*Bk temp];
    temp=[Ak^(i-1)*Bk temp];

    Qt=[Qt zeros(size(Qt,1),size(Q,1));
        zeros(size(Q,1),size(Qt,1)) Q];
    Rt=[Rt zeros(size(Rt,1),size(R,1));
        zeros(size(R,1),size(Rt,1)) R];
end

% 控制量ut的上下限
lb = -pi/2*ones(Np,1);
ub = pi/2*ones(Np,1);

% 控制量ut的初始值
u0=zeros(Np,1);

% 误差量的的初始值
x0=model.location;

% 转换后的优化目标函数矩阵，循环优化函数中H后的表达式为优化目标的另一项
H=2*(Bt'*Qt*Bt + Rt);
delta = 0.001;%避免发散的小量
H= [H zeros(length(H),1);
    zeros(1,length(H)) delta];
% 转换后的优化中的不等式约束左边系数矩阵，后面循环中的bi为不等式右边
Bt = [Bt ones(size(Bt,1),1)*delta];
Ai=[Bt; -Bt];
% 声明u来保存每一步采用的控制量
u=[];
x=x0;
xk=x0;

for k=1:Nc  
    % 关于ut的不等式约束，实际上约束的是状态量，常数4就是状态量约束的上下边界
    bi=[5-At*xk; 5+At*xk];
    % 一切准备就绪，进行二次优化
    [ut, fval, exitflag]=quadprog(H,(2*xk'*At'*Qt*Bt)',Ai,bi,[],[],lb,ub,u0);
    % 显示求解结果是否正常
    fprintf('%d\n', exitflag)
    
    % 采用优化得到的控制量的第一个元素作为实际作用的控制量，代入到原系统中得到下一个时刻的状态量
    u1(k) = ut(1)+Position(k,1);
    u2(k) = ut(2)+Position(k,2);
    x(:, k+1) = Ak*x(:, k) + Bk*[u1(k) u2(k)]';
    xk = x(:, k+1);
    % 对优化初始值进行修改，采用预测值的后段作为下一步的初始值
    u0 = [ut(2:model.Np); ut(model.Np)];
end

cost = sum(sqrt((x(1,:)).^2+(x(2,:)).^2));

sol = reshape(Position, [], 2);
L = cost;
