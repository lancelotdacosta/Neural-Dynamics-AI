function MDP = DEM_demo_MDP_X2
% Demo of active inference for trust games
%__________________________________________________________________________
%
% This routine uses a Markov decision process formulation of active
% inference (with variational Bayes) to model foraging for information in a
% three arm maze.  This demo illustrates variational free energy
% minimisation in the context of Markov decision processes, where the agent
% is equipped with prior beliefs that it will minimise expected free energy
% in the future. This free energy is the free energy of future sensory
% states expected under the posterior predictive distribution. It can be
% regarded as a generalisation of the variational formulation of KL control
% in which information gain or epistemic value is formulated explicitly.
%
% In this example, the agent starts at the centre of a three way maze which
% is baited with a reward in one of the two upper arms. However, the
% rewarded arm changes from trial to trial.  Crucially, the agent can
% identify where the reward (US) is located by accessing a cue (CS) in the
% lower arm. This tells the agent whether the reward is on the left or the
% right upper arm.  This means the optimal policy would first involve
% maximising information gain or epistemic value by moving to the lower arm
% and then claiming the reward this signified. Here, there are eight hidden
% states (four locations times right or left reward), four control states
% (that take the agent to the four locations) and four exteroceptive
% outcomes (that depend on the agents locations) plus three interoceptive
% outcomes indicating reward (or not).
%
% This version focuses on factorising the hidden states causing
% (factorised) outcomes. This factorisation is implicit in the tensor
% production used in the companion demo.  Here the factorisation is
% explicit enabling us to model multiple modalities (outcome factors) and
% distinct hidden causes of observation (hidden state factors like what and
% where). The behaviour is formally similar to the vanilla scheme but
% allows a much more intuitive (and possibly flexible) model specification.
%
% see also: DEM_demo_MDP_habits.m and spm_MPD_VB_X.m
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_demo_MDP_X.m 7319 2018-05-29 09:33:01Z karl $
 
% set up and preliminaries
%==========================================================================
  
% outcome probabilities: A
%--------------------------------------------------------------------------
% We start by specifying the probabilistic mapping from hidden states
% to outcomes; where outcome can be exteroceptive or interoceptive: The
% exteroceptive outcomes A{1} provide cues about location and context,
% while interoceptive outcome A{2} denotes different levels of reward
%--------------------------------------------------------------------------
a      = .98; % probability of where the reward is after seeing a cue?
b      = 1 - a;

% exteroceptive outcomes A{1}: location
A{1}(:,:,1) = [... %(:,:,1) means we have right reward
                    %{1} means
    1 0 0 0;    % cue start %don't know where the reward is
    0 1 0 0;    % cue left  
    0 0 1 0;    % cue right
    0 0 0 1     % cue CS right %after accessing cue
    0 0 0 0];   % cue CS left %after accessing cue
A{1}(:,:,2) = [... %(:,:,2) means we have left reward
    1 0 0 0;    % cue start
    0 1 0 0;    % cue left
    0 0 1 0;    % cue right
    0 0 0 0     % cue CS right
    0 0 0 1];   % cue CS left
 
A{2}(:,:,1) = [... %(:,:,1) means we have right reward
    1 0 0 1;    % reward neutral
    0 a b 0;    % reward positive
    0 b a 0];   % reward negative
A{2}(:,:,2) = [... %(:,:,2) means we have left reward
    1 0 0 1;    % reward neutral
    0 b a 0;    % reward positive
    0 a b 0];   % reward negative
 
 
% controlled transitions: B{u}
%--------------------------------------------------------------------------
% Next, we have to specify the probabilistic transitions of hidden states
% for each factor. Here, there are four actions taking the agent directly
% to each of the four locations.
%--------------------------------------------------------------------------
B{1}(:,:,1)  = [1 0 0 1; 0 1 0 0;0 0 1 0;0 0 0 0];
B{1}(:,:,2)  = [0 0 0 0; 1 1 0 1;0 0 1 0;0 0 0 0];
B{1}(:,:,3)  = [0 0 0 0; 0 1 0 0;1 0 1 1;0 0 0 0];
B{1}(:,:,4)  = [0 0 0 0; 0 1 0 0;0 0 1 0;1 0 0 1];
 
% context, which cannot be changed by action
%--------------------------------------------------------------------------
B{2}  = eye(2);
 
% priors: (utility) C
%--------------------------------------------------------------------------
% Finally, we have to specify the prior preferences in terms of log
% probabilities over outcomes. Here, the agent prefers rewards to losses -
% and does not like to be exposed
%--------------------------------------------------------------------------
C{1}  = [-1 -1 -1;
          0  0  0;
          0  0  0;
          0  0  0;
          0  0  0];
c     = 6;
C{2}  = [ 0  0  0;
          c  c  c;
         -c -c -c];
 
% now specify prior beliefs about initial states, in terms of counts. Here
% the hidden states are factorised into location and context:
%--------------------------------------------------------------------------
d{1} = [128 1 1 0]';
d{2} = [2 2]';
 
 
% allowable policies (of depth T).  These are just sequences of actions
% (with an action for each hidden factor)
%--------------------------------------------------------------------------
V(:,:,1) = [1  1  1  2  3 
            1  2  3  2  3 ];
V(:,:,2) = 1;
 
 
% MDP Structure - this will be used to generate arrays for multiple trials
%==========================================================================
mdp.V = V;                       % allowable policies
mdp.A = A;                       % observation model
mdp.B = B;                       % transition probabilities
mdp.C = C;                       % preferred outcomes
mdp.d = d;                       % prior over initial states
mdp.s = [1 1]';                  % true initial state

mdp.Aname = {'exteroceptive','interoceptive'};
mdp.Bname = {'position','context'};
mdp.tau   = 4;
 
% Function to find correct responses
%--------------------------------------------------------------------------
hit = @(MDP) any(MDP.o(3,:) == 2) & ~any(MDP.o(3,:) == 3); % change this to reflect successes, need to understand task better.

% Solve - an example game: a run of reds then an oddball
%==========================================================================
Ns    = 64; %number of subjects
N     = 24; %number of trials per subject
for m = 1:Ns
    clear MDP NatDP
    
j              = [1,3] ;         % change context in a few of trials
[MDP(1:N), NatDP(1:N)]    = deal(mdp);      
[MDP(j).s, NatDP(j).s]     = deal([1 2]');   % puts initial state [1,2] to the trials i;

    % Simulation with natural gradient
    %----------------------------------------------------------------------
    rng(m) %control random number generation
    [NatDP,~]  = spm_MDP_VB_X_nat3(NatDP);
    for i = 1:N
        nat_length(m,i) = cat_inf_length(NatDP(1,i));
    end
   
    % Find correct responses
    %----------------------------------------------------------------------
%     for i = 1:N
%         if hit(NatDP(i)); r(i,m) = 1; else r(i,m) = 0;  end
%     end
    
    % Simulation without natural gradient
    %----------------------------------------------------------------------
    rng(m) %control random number generation
    MDP  = spm_MDP_VB_X(MDP);
    for i = 1:N
        inf_length(m,i) = cat_inf_length(MDP(1,i));
    end
    
    % Find correct responses
    %----------------------------------------------------------------------
% %     for i = 1:N
% %         if hit(MDP(i)); h(i,m) = 1; else h(i,m) = 0;  end
% %     end
%     
%     %Figures
%     
%     %performance
%     %--------------------
% %     spm_figure('GetWin','Figure 11');clf
% %     subplot(4,1,1)
% %   % b = bar(mean(R(:,:),2)); set(b,'EdgeColor','w','FaceColor',[1 1 1]*.0),hold on
% %     b = bar(mean(H(:,:),2)); set(b,'EdgeColor','w','FaceColor',[1 1 1]*.8);%,hold off
% %     xlabel('trial'), ylabel('probability of correct'), axis([1/2 (N + 1/2) 1/3 1]);
% %     title('Average performance','Fontsize',16)
% %     
% %     subplot(4,1,2)
% %     b = bar(mean(R(:,:),2)); set(b,'EdgeColor','w','FaceColor',[1 1 1]*.0);%,hold on
% %     %b = bar(mean(H(:,:),2)); set(b,'EdgeColor','w','FaceColor',[1 1 1]*.8);%,hold off
% %     xlabel('trial'), ylabel('probability of correct'), axis([1/2 (N + 1/2) 1/3 1]);
% %     title('Natural gradient average performance','Fontsize',16)
% %     
% %     subplot(4,1,3)
% %     b = bar(mean(R(:,:),2)); set(b,'EdgeColor','w','FaceColor',[1 1 1]*.0),hold on
% %     b = bar(mean(H(:,:),2)); set(b,'EdgeColor','w','FaceColor',[1 1 1]*.8),hold off
% %     xlabel('trial'), ylabel('probability of correct'), axis([1/2 (N + 1/2) 1/3 1]);
% %     title('Comparison','Fontsize',16)
% %     
% %     subplot(4,1,4)
% %     b = bar(mean(H(:,:),2)); set(b,'EdgeColor','w','FaceColor',[1 1 1]*.8), hold on
% %     b = bar(mean(R(:,:),2)); set(b,'EdgeColor','w','FaceColor',[1 1 1]*.0), hold off
% %     xlabel('trial'), ylabel('probability of correct'), axis([1/2 (N + 1/2) 1/3 1]);
% %     title('Comparison','Fontsize',16)
% %     
% color = linspecer(2);
%     %information length categorical updates plot
%     %-----------------------------------
%     spm_figure('GetWin','Figure 12');clf
%     subplot(2,1,1)
%     b = bar(mean(nat_length(:,:),1)); set(b,'EdgeColor','w','FaceColor',color(1,:));%, hold off
%     %b = bar(mean(R(:,:),2)); set(b,'EdgeColor','w','FaceColor',[1 1 1]*.0), hold off
%     xlabel('trial'), ylabel('Information length');%, axis([1/2 (N + 1/2) 1/3 1]);
%     title('Active inference','Fontsize',16)
%     
%     subplot(2,1,2)
%     b = bar(mean(inf_length(:,:),1)); set(b,'EdgeColor','w','FaceColor',color(2,:));%, hold off
%     %b = bar(mean(R(:,:),2)); set(b,'EdgeColor','w','FaceColor',[1 1 1]*.0), hold off
%     xlabel('trial'), ylabel('Information length');%, axis([1/2 (N + 1/2) 1/3 1]);
%     title('Natural gradient','Fontsize',16)
%     
%     hold off, drawnow
    
%     %differences in information length boxplot
%     %-----------------------------------
%       spm_figure('GetWin','Figure 13');clf
%       subplot(2,1,1)
%       boxplot(nat_length(:)-inf_length(:))
%       ylabel('Active inference - Natural gradient');
%       title('Performance difference','Fontsize',16)
%       
%       subplot(2,1,2)
%       boxplot(abs(nat_length(:)-inf_length(:)))
%       ylabel('|Active inference - Natural gradient|');
%       title('Performance difference','Fontsize',16)
% 
%     hold off, drawnow
      
    simplex2(NatDP, MDP,m)   
% %     spm_figure('GetWin','Figure 14');clf
% %       
% %       subplot(2,1,1)
% 
%       figure(1), clf,
%       
%       trial = randperm(N,1);
%       i = randperm(3,1);
%       j = randperm(3,1);
%       nat_updates = NatDP(trial).xn{1,1}(:,:,i,j);
%       nat_updates(:,4)=[];
%       inf_updates = MDP(trial).xn{1,1}(:,:,i,j);
%       inf_updates(:,4)=[];
%       
%       line_width = 2;
%         
%       %plotting simplex
%       k=2; %2-simplex
%         simplex_vertices = eye(k+1);
% 
%         % for plotting
%         figure(1), clf,
%         simp_vert = [simplex_vertices, simplex_vertices(:,1)];
%         plot3(simp_vert(1,:), simp_vert(2,:), simp_vert(3,:));
%         hold on
%         
%         color = linspecer(2);
%         plot3(nat_updates(:,1), nat_updates(:,2), nat_updates(:,3), 'LineWidth',line_width, 'Color', color(1,:));
%         plot3(inf_updates(:,1), inf_updates(:,2), inf_updates(:,3), 'LineWidth',line_width, 'Color', color(2,:));
% 
% 
%     hold off, drawnow

    m
    disp('agents done.')
end


 
% % illustrate behavioural responses – first trial
% %--------------------------------------------------------------------------
% spm_figure('GetWin','Figure 1'); clf
% spm_MDP_VB_trial(MDP(1));
%  
% % illustrate behavioural responses and neuronal correlates
% %--------------------------------------------------------------------------
% spm_figure('GetWin','Figure 2'); clf
% spm_MDP_VB_game(MDP);
%  
% % illustrate phase-precession and responses to chosen option - 1st trial
% %--------------------------------------------------------------------------
% spm_figure('GetWin','Figure 3'); clf
% spm_MDP_VB_LFP(MDP(1),[2 3;3 3],1);
%  
% % illustrate phase-amplitude (theta-gamma) coupling
% %--------------------------------------------------------------------------
% spm_figure('GetWin','Figure 4'); clf
% spm_MDP_VB_LFP(MDP(1:8));
% 
%  
% % illustrate familiarity (c.f., MMN) and context learning
% %--------------------------------------------------------------------------
% spm_figure('GetWin','Figure 5'); clf
% i = find(ismember(spm_cat({MDP.u}'),[4 2],'rows')); i = (i + 1)/2;
% spm_MDP_VB_LFP(MDP([i(1),i(end)]),[1;1],2)
% subplot(4,1,1), title('Repetition suppression and DA transfer','FontSize',16)
%  
% spm_figure('GetWin','Figure 6'); clf
% n  = size(MDP(1).xn{1},1);
% v  = spm_MDP_VB_LFP(MDP([i(1),i(end)]),[1;1],2);
% t  = ((1:n)*16 + 80)*16/n;
% subplot(2,1,1),plot(t,v{1}{2,1},'b-.',t,v{2}{2,1},'b:',t,v{2}{2,1} - v{1}{2,1})
% xlabel('Time (ms)'),ylabel('LFP'),title('Difference waveform (MMN)','FontSize',16)
% legend({'oddball','standard','MMN'}), grid on, axis square
% 
% w  = [MDP(i(1)).dn MDP(i(end)).dn];
% 
% subplot(2,1,2),bar(w)
% xlabel('Time (bins)'),ylabel(''),title('Phasic DA responses','FontSize',16)
% legend({'oddball','standard'}), grid on
% 

