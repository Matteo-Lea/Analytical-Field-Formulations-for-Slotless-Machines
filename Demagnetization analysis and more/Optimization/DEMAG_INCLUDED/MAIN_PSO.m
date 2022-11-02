%    [1]Coello, C. A. C., Pulido, G. T., & Lechuga, M. S. (2004). Handling%
%       multiple objectives with particle swarm optimization. IEEE Tran-  %
%       sactions on evolutionary computation, 8(3), 256-279.              %
%                                                                         %
%    [2]Sierra, M. R., & Coello, C. A. C. (2005, March). Improving PSO-   %
%       based multi-objective optimization using crowding, mutation and ?-%
%       dominance. In International Conference on Evolutionary Multi-Crite%
%       rion Optimization (pp. 505-519). Springer Berlin Heidelberg.      %
% ----------------------------------------------------------------------- %
function REP = MAIN_PSO(params,params_mach,Var_max,Var_min,Int_min,Int_max,Const,adm,Constraints)

    % Parameters
    Np      = params.Np;
    el      = Np;
    Nr      = params.Nr;
    maxgen  = params.maxgen;
    W       = params.W;
    C1      = params.C1;
    C2      = params.C2;
    C3      = params.C3;
    C4      = params.C4;
    ngrid   = params.ngrid;
    maxvel  = params.maxvel;
    u_mut   = params.u_mut;
    nVar    = length(Var_max(1,:))+length(Int_max(1,:));
    int     = length(Int_max(1,:));
   	m_PM = params.m_PM; % number of harmonics or series components for the magnetization functions
    r_dis = params.r_dis; % radial discretization of magnets region
    
    % preallocation
    x(el,nVar) = 0;
    integer{int} = [];
    integer_mat{int} = [];
    st_sp_mat{int} = [];
    
    %% Initiation of all the parameters of all the different particles
    x(:,1) = rand(el,1).*(Var_max(:,1)-Var_min(:,1))+Var_min(:,1); % mid-magnet to pole ratio [-]
    x(:,2) = rand(el,1).*(Var_max(:,2)-Var_min(:,2))+Var_min(:,2); % active lenght [m]
    x(:,3) = rand(el,1).*(Var_max(:,3)-Var_min(:,3))+Var_min(:,3); % stator radius (facing winding) [m]
    x(:,4) = rand(el,1).*(Var_max(:,4)-Var_min(:,4))+Var_min(:,4); % winding thickness [m]
    x(:,5) = rand(el,1).*(Var_max(:,5)-Var_min(:,5))+Var_min(:,5); % air-gap thickness [m]
    x(:,6) = rand(el,1).*(Var_max(:,6)-Var_min(:,6))+Var_min(:,6); % magnets thickness [m]
    x(:,7) = rand(el,1).*(Var_max(:,7)-Var_min(:,7))+Var_min(:,7); % phase current density peak value [A]
    % the integer variables are located at the end of the vector
    x(:,8) = sign(rand(el,1).*(Int_max(1,1)-Int_min(1,1))+Int_min(1,1)); % 1 for inrunner -1 for outrunner (this shouldn't vary in the optimizatio n process)
    x(:,9) = randi([Int_min(1,2),Int_max(1,2)],el,1); % pole pairs
    
    % derived parameters
    params_geo.R_sleeve = x(:,3) - x(:,8).*x(:,6)- x(:,8).*x(:,5)- x(:,8).*x(:,4); % magnets array radius (facing magnets support)[m]
    % The chosen initialization allows R_sleeve to assume negative values. Theese
    % values are therefore corrected
    x((params_geo.R_sleeve<0),3) = x((params_geo.R_sleeve<0),3) + 2*abs(params_geo.R_sleeve(params_geo.R_sleeve<0));
    % This definition of R_i shouldn't vary iteratively in the optimiation
    params_geo.R_i =x(:,3) - x(:,8).*x(:,4)- x(:,8).*x(:,6)- x(:,8).*x(:,5); % infinite permeability boundary behind the magnets [m] (magnetic support is default)
    % process to keep the geometries with and without backing iron
%     params_geo.R_i(params_geo.In_Out==1) = params_geo.R_i(params_geo.In_Out==1).*randi([0,1],nnz(params_geo.In_Out==1),1); % non-magnetic support set for some in-runner topologies (randomly)
%     params_geo.R_i(params_geo.In_Out==-1) = params_geo.R_i(params_geo.In_Out==-1).*1./randi([0,1],nnz(params_geo.In_Out==-1),1); % non-magnetic support set for some out-runner topologies (randomly)
%     params_geo.R_sleeve = (x(:,3) - x(:,10).*x(:,7)- x(:,10).*x(:,6)- x(:,10).*x(:,5)- x(:,10).*x(:,8)).*(params_geo.R_i~=0&params_geo.R_i~=Inf); % magnets array radius (facing magnets support)[m]

    var_max = Var_max;
    var_min = Var_min;
    
    p_f = 1; % penalty factor
    toll = 0.005; % target torque constraint margin
    
    
    [M_r_n,M_theta_n] = Magnetization(el,m_PM,x,Const);
    [t_bi,peak,T_avg,t_spec,Magnets_weight,Joule,Tot_weight,Weighted_Perc_demag,Weight_mid,Weight_side] = Field_Solution_Demag(el,m_PM,r_dis,x,params_geo,M_r_n,M_theta_n,params_mach,Const);
    f_lim = [ T_avg, peak,Weighted_Perc_demag];
    f_lim_1 = [ max(abs(f_lim(:,1)-Constraints(1))-toll*Constraints(1),0), max(f_lim(:,2)-Constraints(2),0)];
%     f_lim_2 = p_f*max(f_lim(:,3)-Constraints(3),0);

        % multiple objective optimization example
%         f = [Tot_weight + sum(f_lim_1,2),Weighted_Perc_demag];
        % single objective optimization example
         f = Tot_weight + sum(f_lim_1,2)+Weighted_Perc_demag;
        % functions to be saved [specific torque, active weight, back-iron thickness] 
        real = [ -t_spec,Tot_weight,t_bi];
    % Initialization
    POS = x;
    VEL = zeros(Np,nVar-int); % velocity is solely for the non-integer variables
    POS_fit  = f;
    if size(POS,1) ~= size(POS_fit,1)
        warning(['The objective function is badly programmed. It is not returning' ...
            'a value for each particle, please check it.']);
    end
    
    
    
    PBEST    = POS;
    PBEST_fit= POS_fit;
    DOMINATED= checkDomination(POS_fit);
    REP.pos  = POS(~DOMINATED,:);
    REP.pos_fit = POS_fit(~DOMINATED,:);
    REP.lim = f_lim(~DOMINATED,:);
    REP.real = real(~DOMINATED,:);
    REP      = updateGrid(REP,ngrid);
    maxvel   = (var_max-var_min).*maxvel./100;
    gen      = 1;
    
    % Plotting and verbose
    if(size(POS_fit,2)==2)
        h_fig = figure(1);
        h_par = plot(POS_fit(:,1),POS_fit(:,2),'or'); hold on;
        h_rep = plot(REP.pos_fit(:,1),REP.pos_fit(:,2),'ok'); hold on;
        try
            set(gca,'xtick',REP.hypercube_limits(:,1)','ytick',REP.hypercube_limits(:,2)');
            axis([min(REP.hypercube_limits(:,1)) max(REP.hypercube_limits(:,1)) ...
                  min(REP.hypercube_limits(:,2)) max(REP.hypercube_limits(:,2))]);
            grid on; xlabel('f1'); ylabel('f2');
        end
        drawnow;
    end
    if(size(POS_fit,2)==3)
        h_fig = figure(1);
        h_par = plot3(POS_fit(:,1),POS_fit(:,2),POS_fit(:,3),'or'); hold on;
        h_rep = plot3(REP.pos_fit(:,1),REP.pos_fit(:,2),REP.pos_fit(:,3),'ok'); hold on;
        try
                set(gca,'xtick',REP.hypercube_limits(:,1)','ytick',REP.hypercube_limits(:,2)','ztick',REP.hypercube_limits(:,3)');
                axis([min(REP.hypercube_limits(:,1)) max(REP.hypercube_limits(:,1)) ...
                      min(REP.hypercube_limits(:,2)) max(REP.hypercube_limits(:,2))]);
        end
        grid on; xlabel('f1'); ylabel('f2'); zlabel('f3');
        drawnow;
        axis square;
    end
    display(['Generation #0 - Repository size: ' num2str(size(REP.pos,1))]);
    
    % Main MPSO loop
    stopCondition = false;
    while ~stopCondition
        
        % Select leader
        h = selectLeader(REP);
        
        % Update speeds and positions
        VEL = W.*VEL + C1*rand(Np,nVar-int).*(PBEST(:,1:end-int)-POS(:,1:end-int))...
                     + C2*rand(Np,nVar-int).*(repmat(REP.pos(h,1:end-int),Np,1)-POS(:,1:end-int));
        POS(:,1:end-int) = POS(:,1:end-int) + VEL;
        
        % EVOLUTIONARY STRATEGY FOR THE DISCRETE VARIABLES 
        if gen == 1  
            st_sp = 1./adm.Int; % standard spacing for each possible values of each integer variable   
            for ii = 1:int
                integer{ii} = linspace(Int_min(1,ii),Int_max(1,ii),adm.Int(ii)); % vector holding all the admissible values for each integer variable
                integer_mat{ii} = repmat(integer{ii},el,1); 
                st_sp_mat{ii} = ones(el,adm.Int(ii))*st_sp(ii); % matrix holding the standard spacings between the admissible values for each integer variable
            end
        end 
            for ii = 1:int
                best_ind = (POS(:,end-int+ii)==PBEST(:,end-int+ii)); % vector holding ones for the particles located at their best position so far and zero elsewhere
                Ind_idx = find((POS(:,end-int+ii).*best_ind)==repmat(integer{ii},el,1)); % index of the individuals' best position so far
                st_sp_mat{ii}(Ind_idx) = C4*st_sp_mat{ii}(Ind_idx);
                
                Ind_glo = find(REP.pos(h,end-int+ii)==repmat(integer{ii},el,1)); % index of the global best position so far
                st_sp_mat{ii}(Ind_glo) = C3*st_sp_mat{ii}(Ind_glo);
                
                st_sp_mat{ii} = st_sp_mat{ii}./sum(st_sp_mat{ii},2); % normalization of the spacing matrix
                Intervals = cumsum(st_sp_mat{ii},2); % intervals between the admissible values (from 0 to 1)
                [~, ind_count] = max((rand(el,1)<=Intervals)~=0, [], 2 ,'linear');
                POS(:,end-int+ii) = integer_mat{ii}(ind_count);
%                 rand_int = integer_mat{ii}(find(rand(el,1)<=Intervals));
            end
        
        
        
        % Perform mutation
        POS = mutation(POS,gen,maxgen,Np,Var_max,Var_min,Int_max,Int_min,nVar,u_mut);
        
        % Check boundaries
        [POS,VEL] = checkBoundaries(POS,VEL,maxvel,var_max,var_min,int);  
        
        % adjustment of some geometrical limits (R_s) for the particles
        % coordinates to make sense
        params_geo.R_sleeve = POS(:,3) - POS(:,8).*POS(:,4)- POS(:,8).*POS(:,6)- POS(:,8).*POS(:,5); % magnets array radius (facing magnets support)[m]
        % The chosen initialization allows R_sleeve to assume negative values. Theese
        % values are therefore corrected
        POS((params_geo.R_sleeve<0),3) = POS((params_geo.R_sleeve<0),3) + 2*abs(params_geo.R_sleeve(params_geo.R_sleeve<0));
        % This definition of R_i shouldn't vary iteratively in the optimiation
        params_geo.R_i =POS(:,3) - POS(:,8).*POS(:,4)- POS(:,8).*POS(:,6)- POS(:,8).*POS(:,5); % infinite permeability boundary behind the magnets [m] (magnetic support is default)
        % process to keep the geometries with and without backing iron
    %     params_geo.R_i(params_geo.In_Out==1) = params_geo.R_i(params_geo.In_Out==1).*randi([0,1],nnz(params_geo.In_Out==1),1); % non-magnetic support set for some in-runner topologies (randomly)
    %     params_geo.R_i(params_geo.In_Out==-1) = params_geo.R_i(params_geo.In_Out==-1).*1./randi([0,1],nnz(params_geo.In_Out==-1),1); % non-magnetic support set for some out-runner topologies (randomly)
    
        
        % Evaluate the population
        [M_r_n,M_theta_n] = Magnetization(el,m_PM,POS,Const);
        tic
        [t_bi,peak,T_avg,t_spec,Magnets_weight,Joule,Tot_weight,Weighted_Perc_demag,Weight_mid,Weight_side] =  Field_Solution_Demag(el,m_PM,r_dis,POS,params_geo,M_r_n,M_theta_n,params_mach,Const);
        toc
        f_lim = [ T_avg, peak,Weighted_Perc_demag];
        f_lim_1 = [ max(abs(f_lim(:,1)-Constraints(1))-toll*Constraints(1),0), max(f_lim(:,2)-Constraints(2),0)];
%         f_lim_2 = p_f*max(f_lim(:,3)-Constraints(3),0);

        % multiple objective optimization example
%         f = [Tot_weight + sum(f_lim_1,2),Weighted_Perc_demag];
        % single objective optimization example
         f = Tot_weight + sum(f_lim_1,2)+Weighted_Perc_demag;
        % functions to be saved [specific torque, active weight, back-iron thickness] 
        real = [ -t_spec,Tot_weight,t_bi];

        
        
        POS_fit = f;
        
        % Update the repository
        REP = updateRepository(REP,POS,POS_fit,real,f_lim,ngrid);
        if(size(REP.pos,1)>Nr)
             REP = deleteFromRepository(REP,size(REP.pos,1)-Nr,ngrid);
        end
        
        % Update the best positions found so far for each particle
        pos_best = dominates(POS_fit, PBEST_fit);
        best_pos = ~dominates(PBEST_fit, POS_fit);
        best_pos(rand(Np,1)>=0.5) = 0;
        if(sum(pos_best)>1)
            PBEST_fit(pos_best,:) = POS_fit(pos_best,:);
            PBEST(pos_best,:) = POS(pos_best,:);
        end
        if(sum(best_pos)>1)
            PBEST_fit(best_pos,:) = POS_fit(best_pos,:);
            PBEST(best_pos,:) = POS(best_pos,:);
        end
        
        % Plotting and verbose
        if(size(POS_fit,2)==2)
            figure(h_fig); delete(h_par); delete(h_rep);
            h_par = plot(POS_fit(:,1),POS_fit(:,2),'or'); hold on;
            h_rep = plot(REP.pos_fit(:,1),REP.pos_fit(:,2),'ok'); hold on;
            try
                set(gca,'xtick',REP.hypercube_limits(:,1)','ytick',REP.hypercube_limits(:,2)');
                axis([min(REP.hypercube_limits(:,1)) max(REP.hypercube_limits(:,1)) ...
                      min(REP.hypercube_limits(:,2)) max(REP.hypercube_limits(:,2))]);
            end
%            if(isfield(MultiObj,'truePF'))
%                 try delete(h_pf); end
%                 h_pf = plot(MultiObj.truePF(:,1),MultiObj.truePF(:,2),'.','color',0.8.*ones(1,3)); hold on;
%             end
            grid on; xlabel('f1'); ylabel('f2');
            drawnow;
            axis square;
        end
        if(size(POS_fit,2)==3)
            figure(h_fig); delete(h_par); delete(h_rep); 
            h_par = plot3(POS_fit(:,1),POS_fit(:,2),POS_fit(:,3),'or'); hold on;
            h_rep = plot3(REP.pos_fit(:,1),REP.pos_fit(:,2),REP.pos_fit(:,3),'ok'); hold on;
            try
                set(gca,'xtick',REP.hypercube_limits(:,1)','ytick',REP.hypercube_limits(:,2)','ztick',REP.hypercube_limits(:,3)');
                axis([min(REP.hypercube_limits(:,1)) max(REP.hypercube_limits(:,1)) ...
                      min(REP.hypercube_limits(:,2)) max(REP.hypercube_limits(:,2)) ...
                      min(REP.hypercube_limits(:,3)) max(REP.hypercube_limits(:,3))]);
            end
%             if(isfield(MultiObj,'truePF'))
%                 try delete(h_pf); end
%                 h_pf = plot3(MultiObj.truePF(:,1),MultiObj.truePF(:,2),MultiObj.truePF(:,3),'.','color',0.8.*ones(1,3)); hold on;
%             end
            grid on; xlabel('f1'); ylabel('f2'); zlabel('f3');
            drawnow;
            axis square;
        end
        display(['Generation #' num2str(gen) ' - Repository size: ' num2str(size(REP.pos,1))]);
        
        % Update generation and check for termination
        gen = gen + 1;
        if(gen>maxgen), stopCondition = true; end
    end
    hold off;
end

% Function that updates the repository given a new population and its
% fitness
function REP = updateRepository(REP,POS,POS_fit,real,f_lim,ngrid)
    % Domination between particles
    DOMINATED  = checkDomination(POS_fit);
    REP.pos    = [REP.pos; POS(~DOMINATED,:)];
    REP.pos_fit= [REP.pos_fit; POS_fit(~DOMINATED,:)];
    REP.lim= [REP.lim; f_lim(~DOMINATED,:)];
    REP.real= [REP.real; real(~DOMINATED,:)];
    % Domination between nondominated particles and the last repository
    DOMINATED  = checkDomination(REP.pos_fit);
    
    REP.pos_fit= REP.pos_fit(~DOMINATED,:);
    REP.pos    = REP.pos(~DOMINATED,:);
    REP.lim    = REP.lim(~DOMINATED,:);
    REP.real    = REP.real(~DOMINATED,:);
    % Updating the grid
    REP        = updateGrid(REP,ngrid);
end

% Function that corrects the positions and velocities of the particles that
% exceed the boundaries
function [POS,VEL] = checkBoundaries(POS,VEL,maxvel,var_max,var_min,int)
    % Useful matrices
    Np = size(POS(:,1:end-int),1);
%     MAXLIM   = repmat(var_max(:)',Np,1);
%     MINLIM   = repmat(var_min(:)',Np,1);
%     MAXVEL   = repmat(maxvel(:)',Np,1);
%     MINVEL   = repmat(-maxvel(:)',Np,1);
    MAXLIM   = var_max;
    MINLIM   = var_min;
    MAXVEL   = maxvel;
    MINVEL   = -maxvel;
    
    
    % Correct positions and velocities
    VEL(VEL>MAXVEL) = MAXVEL(VEL>MAXVEL);
    VEL(VEL<MINVEL) = MINVEL(VEL<MINVEL);
    VEL(POS(:,1:end-int)>MAXLIM) = (-1).*VEL(POS(:,1:end-int)>MAXLIM);
    POS(POS(:,1:end-int)>MAXLIM) = MAXLIM(POS(:,1:end-int)>MAXLIM);
    VEL(POS(:,1:end-int)<MINLIM) = (-1).*VEL(POS(:,1:end-int)<MINLIM);
    POS(POS(:,1:end-int)<MINLIM) = MINLIM(POS(:,1:end-int)<MINLIM);
end

% Function for checking the domination between the population. It
% returns a vector that indicates if each particle is dominated (1) or not
function dom_vector = checkDomination(fitness)
    Np = size(fitness,1);
    dom_vector = zeros(Np,1);
    all_perm = nchoosek(1:Np,2);    % Possible permutations
    all_perm = [all_perm; [all_perm(:,2) all_perm(:,1)]];
    
    d = dominates(fitness(all_perm(:,1),:),fitness(all_perm(:,2),:));
    dominated_particles = unique(all_perm(d==1,2));
    dom_vector(dominated_particles) = 1;
end

% Function that returns 1 if x dominates y and 0 otherwise
function d = dominates(x,y)
    d = all(x<=y,2) & any(x<y,2);
end


% Function that updates the hypercube grid, the hypercube where belongs
% each particle and its quality based on the number of particles inside it
function REP = updateGrid(REP,ngrid)
    % Computing the limits of each hypercube
    ndim = size(REP.pos_fit,2);
    REP.hypercube_limits = zeros(ngrid+1,ndim);
    for dim = 1:1:ndim
        REP.hypercube_limits(:,dim) = linspace(min(REP.pos_fit(:,dim)),max(REP.pos_fit(:,dim)),ngrid+1)';
    end
    
    % Computing where belongs each particle
    npar = size(REP.pos_fit,1);
    REP.grid_idx = zeros(npar,1);
    REP.grid_subidx = zeros(npar,ndim);
    for n = 1:1:npar
        idnames = [];
        for d = 1:1:ndim
            REP.grid_subidx(n,d) = find(REP.pos_fit(n,d)<=REP.hypercube_limits(:,d)',1,'first')-1;
            if(REP.grid_subidx(n,d)==0), REP.grid_subidx(n,d) = 1; end
            idnames = [idnames ',' num2str(REP.grid_subidx(n,d))];
        end
        if ndim ==1
            REP.grid_idx(n) = eval(['sub2ind([ngrid,1]' idnames ');']);  
        else
            REP.grid_idx(n) = eval(['sub2ind(ngrid.*ones(1,ndim)' idnames ');']);
        end
    end
    
    % Quality based on the number of particles in each hypercube
    REP.quality = zeros(ngrid,2);
    ids = unique(REP.grid_idx);
    for i = 1:length(ids)
        REP.quality(i,1) = ids(i);  % First, the hypercube's identifier
        REP.quality(i,2) = 10/sum(REP.grid_idx==ids(i)); % Next, its quality
    end
end

% Function that selects the leader performing a roulette wheel selection
% based on the quality of each hypercube
function selected = selectLeader(REP)
    % Roulette wheel
    prob    = cumsum(REP.quality(:,2));     % Cumulated probs
    sel_hyp = REP.quality(find(rand(1,1)*max(prob)<=prob,1,'first'),1); % Selected hypercube
    
    % Select the index leader as a random selection inside that hypercube
    idx      = 1:1:length(REP.grid_idx);
    selected = idx(REP.grid_idx==sel_hyp);
    selected = selected(randi(length(selected)));
end

% Function that deletes an excess of particles inside the repository using
% crowding distances
function REP = deleteFromRepository(REP,n_extra,ngrid)
    % Compute the crowding distances
    crowding = zeros(size(REP.pos,1),1);
    for m = 1:1:size(REP.pos_fit,2)
        [m_fit,idx] = sort(REP.pos_fit(:,m),'ascend');
        m_up     = [m_fit(2:end); Inf];
        m_down   = [Inf; m_fit(1:end-1)];
        distance = (m_up-m_down)./(max(m_fit)-min(m_fit));
        [~,idx]  = sort(idx,'ascend');
        crowding = crowding + distance(idx);
    end
    crowding(isnan(crowding)) = Inf;
    
    % Delete the extra particles with the smallest crowding distances
    [~,del_idx] = sort(crowding,'ascend');
    del_idx = del_idx(1:n_extra);
    REP.pos(del_idx,:) = [];
    REP.pos_fit(del_idx,:) = [];
    REP.lim(del_idx,:) = [];
    REP.real(del_idx,:) = [];
    REP = updateGrid(REP,ngrid); 
end

% Function that performs the mutation of the particles depending on the
% current generation
function POS = mutation(POS,gen,maxgen,Np,Var_max,Var_min,Int_max,Int_min,nVar,u_mut)
    % Sub-divide the swarm in three parts [2]
    fract     = Np/3 - floor(Np/3);
    if(fract<0.5), sub_sizes =[ceil(Np/3) round(Np/3) round(Np/3)];
    else           sub_sizes =[round(Np/3) round(Np/3) floor(Np/3)];
    end
    cum_sizes = cumsum(sub_sizes);
    
    % First part: no mutation
    % Second part: uniform mutation
    nmut = round(u_mut*sub_sizes(2));
    if(nmut>0)
        idx = cum_sizes(1) + randperm(sub_sizes(2),nmut);
%         POS(idx,:) = repmat((var_max-var_min)',nmut,1).*rand(nmut,nVar) + repmat(var_min',nmut,1);
%         POS(idx,:) = (var_max(idx,:)-var_min(idx,:)).*rand(nmut,nVar) + var_min(idx,:);
        POS(idx,1) = rand(nmut,1).*(Var_max(idx,1)-Var_min(idx,1))+Var_min(idx,1); % mid-magnet to pole ratio [-]
        POS(idx,2) = rand(nmut,1).*(Var_max(idx,2)-Var_min(idx,2))+Var_min(idx,2); % active lenght [m]
        POS(idx,3) = rand(nmut,1).*(Var_max(idx,3)-Var_min(idx,3))+Var_min(idx,3); % stator radius (facing winding) [m]
        POS(idx,4) = rand(nmut,1).*(Var_max(idx,4)-Var_min(idx,4))+Var_min(idx,4); % winding thickness [m]
        POS(idx,5) = rand(nmut,1).*(Var_max(idx,5)-Var_min(idx,5))+Var_min(idx,5); % air-gap thickness [m]
        POS(idx,6) = rand(nmut,1).*(Var_max(idx,6)-Var_min(idx,6))+Var_min(idx,6); % magnets thickness [m]
        POS(idx,7) = rand(nmut,1).*(Var_max(idx,7)-Var_min(idx,7))+Var_min(idx,7); % phase current peak value [A] 
        % the integer variables are located at the end of the vector
        POS(idx,8) = sign(rand(nmut,1).*(Int_max(idx,1)-Int_min(idx,1))+Int_min(idx,1)); % 1 for inrunner -1 for outrunner (this shouldn't vary in the optimizatio n process)
        POS(idx,9) = randi([Int_min(1,2),Int_max(1,2)],nmut,1); % pole pairs
    end
    
    % Third part: non-uniform mutation
    per_mut = (1-gen/maxgen)^(5*nVar);     % Percentage of mutation
    nmut    = round(per_mut*sub_sizes(3));
    if(nmut>0)
        idx = cum_sizes(2) + randperm(sub_sizes(3),nmut);
%         POS(idx,:) = repmat((var_max-var_min)',nmut,1).*rand(nmut,nVar) + repmat(var_min',nmut,1);
        POS(idx,1) = rand(nmut,1).*(Var_max(idx,1)-Var_min(idx,1))+Var_min(idx,1); % mid-magnet to pole ratio [-]
        POS(idx,2) = rand(nmut,1).*(Var_max(idx,2)-Var_min(idx,2))+Var_min(idx,2); % active lenght [m]
        POS(idx,3) = rand(nmut,1).*(Var_max(idx,3)-Var_min(idx,3))+Var_min(idx,3); % stator radius (facing winding) [m]
        POS(idx,4) = rand(nmut,1).*(Var_max(idx,4)-Var_min(idx,4))+Var_min(idx,4); % winding thickness [m]
        POS(idx,5) = rand(nmut,1).*(Var_max(idx,5)-Var_min(idx,5))+Var_min(idx,5); % air-gap thickness [m]
        POS(idx,6) = rand(nmut,1).*(Var_max(idx,6)-Var_min(idx,6))+Var_min(idx,6); % magnets thickness [m]
        POS(idx,7) = rand(nmut,1).*(Var_max(idx,7)-Var_min(idx,7))+Var_min(idx,7); % phase current peak value [A]
        % the integer variables are located at the end of the vector
        POS(idx,8) = sign(rand(nmut,1).*(Int_max(idx,1)-Int_min(idx,1))+Int_min(idx,1)); % 1 for inrunner -1 for outrunner (this shouldn't vary in the optimizatio n process)
        POS(idx,9) = randi([Int_min(1,2),Int_max(1,2)],nmut,1); % pole pairs
    end
end