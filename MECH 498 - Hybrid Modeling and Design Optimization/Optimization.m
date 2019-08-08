%% A SCRIPT TO HOLD THE OPTIMIZATION CODE

fun = @SurrocketOpt;
% Design parameters are as follows:
%   V_tank, C_inj, L, d_port_init, d_th, Area_ratio_nozzle
% lb = [0.004, 10e-6, 0.2, 0.030, 0.02, 3];
% ub = [0.014, 40e-6, 0.6, 0.080, 0.04, 10];
lb = zeros(6,1); % Upper and lower bounds handled internally to improve optimization performance
ub = [6,40,33,1,1,3]; % Weight different parameters according to their effect on objective function UB MUST BE THE SAME HERE AS IN SURROCKETOPT!!!

% load('trialPoints_ga_2019-04-20'); % For evaluation of global patternsearch results

% % Trial Points for Initial Design Selection
% % Trial 1
% x0 = [3.520982997972654,4.991298779553515,11.817555019680903,4.832618428056022,0.501289740415362,0.290456405544833,1.047985791389404];
% % Trial 2
% x0 = [3.498155668156791,4.982939984943386,11.473489106733345,7.117595596711158,0.618911028260646,0.375191020327446,1.628379863306198];
% % Trial 3
% x0 = [3.571995032727452,4.999921021890583,9.941143634521358,3.252181497638357,0.527046261790828,0.080927547193480,1.730744785028823];

% Final Hybrid Optimization
x0 = [3.225862085271187,15.088164637286543,4.251160607304687,0.650685032581768,0.286553121742001,1.279248184119617]; % Result from global patternsearch

fmincon_options = optimoptions('fmincon'); %, 'MaxFunctionEvaluations', 300);
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun,x0,[],[],[],[],lb,ub,[],fmincon_options);

% patternsearch_options = optimoptions('patternsearch'); %, 'MaxTime', 300);
% [x,fval,exitflag,output] = patternsearch(fun,x0,[],[],[],[],lb,ub,[],patternsearch_options);

% ga_options = optimoptions('ga', 'MaxGenerations', 3, 'PopulationSize', 100); % Point Generation for patternsearch
% % ga_options = optimoptions('ga', 'MaxTime', 300, 'PopulationSize', 30, 'PlotFcn', 'gaplotbestf'); % Testing
% % ga_options = optimoptions('ga', 'PopulationSize', 400, 'CrossoverFraction', 0.65, 'MaxStallGenerations', 10, 'FunctionTolerance', 1e-4, 'MaxGenerations', 80, 'PlotFcn', 'gaplotbestf'); % Final
% [x,fval,exitflag,output,population,scores] = ga(fun,length(lb),[],[],[],[],lb,ub,[],[],ga_options);

% simulannealbnd_options = optimoptions('simulannealbnd', 'MaxTime', 300, 'MaxFunctionEvaluations', Inf);
% [x,fval,exitflag,output] = simulannealbnd(fun,x0,lb,ub,simulannealbnd_options);

% % particleswarm_options = optimoptions('particleswarm', 'InitialSwarmSpan', 50, 'MaxTime', 300);
% particleswarm_options = optimoptions('particleswarm', 'InitialSwarmSpan', 800);
% [x,fval,exitflag,output] = particleswarm(fun,length(lb),lb,ub,particleswarm_options);

% patternsearch_options = optimoptions('patternsearch', 'MeshTolerance', 1e-2);
% for k = 1:length(trial) % Run patternsearch repeatedly on trial points
% %     [x(k,:),fval(k),exitflag(k),output{k}] = patternsearch(fun,trial(k,:),[],[],[],[],lb,ub,[],patternsearch_options);
%     [I_sp_g(k), max_alt_g(k), vel_off_rail_g(k), p_g{k}] = SurrocketOpt(x(k,:));
% end



% load handel;
% player = audioplayer(y, Fs); 
% play(player);