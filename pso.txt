function [ GlobalBest  , BestCost ] =  pso( param   , model , CostFunction )

nVar = model.nVar ;       %  节点个数
VarSize=[1 nVar];   %  决策变量维度

VarMin = 0 ; % 自变量取值   低阶
VarMax= 0.1 ;% 自变量取值   高阶
%% PSO Parameters
MaxIt=  param.MaxIt ;          %    迭代次数
nPop= param.nPop ;           %  种群数目
w= param.w ;                %   权重
wdamp=param.wdamp ;         %    退化率
c1= param.c1 ;             % Personal Learning Coefficient
c2= param.c2 ;             % Global Learning Coefficient

VelMax =  0.01* (VarMax - VarMin );    %   X  最大  速度
VelMin =  -VelMax ;                    %    X 最小速度


%  模拟退火参数
MaxSubIt= param.MaxSubIt ;    
T =  param.T ;  %
alpha1= param.alpha1;    

%%  初始化
rand('seed', sum(clock));

empty_particle.Position=[];
empty_particle.Velocity=[];
empty_particle.Cost=[];
empty_particle.Sol=[];
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];
empty_particle.Best.Sol=[];


GlobalBest.Cost=inf;


particle=repmat(empty_particle,nPop,1);


for i=1:nPop
    
    particle(i).Position = unifrnd(VarMin, VarMax, VarSize);

    particle(i).Velocity = unifrnd(VelMin, VelMax, VarSize);

    [particle(i).Cost, particle(i).Sol]=CostFunction(particle(i).Position);
    

    particle(i).Best.Position=particle(i).Position;
    particle(i).Best.Cost=particle(i).Cost;
    particle(i).Best.Sol=particle(i).Sol;
    
 
    if particle(i).Best.Cost < GlobalBest.Cost
        
        GlobalBest = particle(i).Best;
        
    end
    
end

% 最优解记录
BestCost=zeros(MaxIt,1);

%% PSO Main Loop

%%
    
 realmin=exp(-10*5);
for it=1:MaxIt
    for i=1:nPop
         
         C1 = exp( ( particle(i).Best.Position -  GlobalBest.Position)/ abs( GlobalBest.Position + realmin ) );
         C2 = exp( ( -particle(i).Best.Position + GlobalBest.Position));
         C3 = zeros(1,length(C2));
         for m = 1:1:length(C2)
            C3(m) = levy();
         end
        % Update Velocity
        particle(i).Velocity  = w*particle(i).Velocity ...
            + C1*rand(VarSize).*(particle(i).Best.Position-particle(i).Position) ...
            + C2.*rand(VarSize).*(GlobalBest.Position-particle(i).Position);
        
        % Update Velocity Bounds
        particle(i).Velocity = max(particle(i).Velocity,VelMin);
        particle(i).Velocity = min(particle(i).Velocity,VelMax);
        
        % Update Position
        particle(i).Position = particle(i).Position + particle(i).Velocity + C3.*rand(VarSize).*(GlobalBest.Position-particle(i).Position);
        
        % Update Position Bounds
        particle(i).Position = max(particle(i).Position,VarMin);
        particle(i).Position = min(particle(i).Position,VarMax);
        
        
        % Evaluation
        [particle(i).Cost, particle(i).Sol]=CostFunction(particle(i).Position);
        
        % Update Personal Best
        if particle(i).Cost<particle(i).Best.Cost
            
            particle(i).Best.Position=particle(i).Position;
            particle(i).Best.Cost=particle(i).Cost;
            particle(i).Best.Sol=particle(i).Sol;
            
            % Update Global Best
            if particle(i).Best.Cost < GlobalBest.Cost
                GlobalBest=particle(i).Best;
            end
            
        end
    end
    %%
    
    % 更新最优解
    BestCost(it)=GlobalBest.Cost;
    
    %权重更新
    w=w*wdamp;
    
    
    disp(it);
end


