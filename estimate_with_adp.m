clear;
addpath 'C:\Users\ekhmelni\Dropbox\Seeds and Drugs\Pharma Projects\KatyaJMP\DrugAttrition_Model\Estimation'
%%%%%%%%%%%%%%%%%%%%%%%DESCRIPTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code implements the value funciton approximation approach for ONE
%disease using the TRUE transition probabilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%Checking
    % folder = "C:\Users\ekhmelni\Dropbox\Seeds and Drugs\Pharma Projects\KatyaJMP\DrugAttrition_Model\Estimation\texfiles";
    % %Main data:
    % filename=fullfile(folder,"topandsmall_data_flexible_98_upd.txt");
    % data_orig=tdfread(filename,'tab');
    % vsly=369000;
    % %Scaling the DALYs back:
    % dalys=data_orig.marketsize*1000;
    % marksize=dalys*vsly/(10e+8); %(in billions)
    % tmp_1 = max(marksize(data_orig.indication_7==1,:));
    % tmp_2 = log(tmp_1);
    % 
    % monopoly_profit_1=(1/rho).*(1-exp(-rho)).*(theta_0*10e+8 + tmp_1.*theta_r*10e+8 + log(1+1).*theta_s*10e+8)./10e+8;
    % monopoly_profit_2=(1/rho).*(1-exp(-rho)).*(theta_0*10e+8 + tmp_2.*theta_r*10e+8 + log(1+1).*theta_s*10e+8)./10e+8;

%Loading the list of parameters
folder = "C:\Users\ekhmelni\Dropbox\Seeds and Drugs\Pharma Projects\KatyaJMP\DrugAttrition_Model\Estimation\resultsother";
filename=fullfile(folder,"params_for_counterfactual_all_3.csv"); 
params_all_table=readtable(filename);
params_use=params_all_table(params_all_table.to_use==1 , :);
params=table2struct(params_use,'ToScalar',true);
%Science data:
filename=fullfile(folder,"data_science_generated_4_noop_2_params.csv");
science_data=readtable(filename);

%Some supplemental arrays
state_transition_names={'from0to1','from0to2','from0to3','from0to4','from1to2','from1to0','from2to3','from2to0','from3to4','from3to0','from4to0'};
stage_from={0,0,0,0,1,1,2,2,3,3,4}';
stage_to={1,2,3,4,2,0,3,0,4,0,0}';
ctp_stage_name={'for1','for2','for3','for4'};


%Number of diseases:
n_indications=size(params.indications_num,1);
%I will run the value function approximation for Epilepsy - indication 10%
indication_num_use=10;
%Max firms for each stage:
max_firms_all=[params.max_in_1,params.max_in_2,params.max_in_3,params.max_in_4];
max_firms=max_firms_all(params.indications_num==indication_num_use,:);
%Market size
marketsize = params.marketsize(params.indications_num==indication_num_use,:);

%Creating the state space
%Matlab first iterates the 1-st element in the vector, then second, etc.
%But Python does it the other way round. So I need to re-order.
s_1=0:max_firms(1);
s_2=0:max_firms(2);
s_3=0:max_firms(3);
s_4=0:max_firms(4);
supportS=transpose(combvec(s_4,s_3,s_2,s_1));
supportS=[supportS(:,4) supportS(:,3) supportS(:,2) supportS(:,1)];
clear s_1 s_2 s_3 s_4;


%%%Known parameters (estimated at step 1)%%%

%Stage durations:
lambdas.lambda_1=1/params.stageduration_1_years(params.indications_num==indication_num_use,:);
lambdas.lambda_2=1/params.stageduration_2_years(params.indications_num==indication_num_use,:);
lambdas.lambda_3=1/params.stageduration_3_years(params.indications_num==indication_num_use,:);
lambdas.lambda_4=1/params.stageduration_4_years(params.indications_num==indication_num_use,:);
%Generic entry:
lambdas.lambda_4_e=1/params.years_till_generic(params.indications_num==indication_num_use,:);


%Entry rates:
entry_lambdas.lambda_1 = 1/params.entryduration_1_years(params.indications_num==indication_num_use,:);
entry_lambdas.lambda_2 = 1/params.entryduration_2_years(params.indications_num==indication_num_use,:);
entry_lambdas.lambda_3 = 1/params.entryduration_3_years(params.indications_num==indication_num_use,:);
entry_lambdas.lambda_4 = 0;

%Discount rate:
rho = params.rho(params.indications_num==indication_num_use,:);

%Generic entry:
years_till_generic = params.years_till_generic(params.indications_num==indication_num_use,:);

%Probability of the FDA approval (e.g., transitioning from stage 3 to stage 4)
p_3 = params.p_3(params.indications_num==indication_num_use,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Aggregating the known parameters%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params_known.p_3=p_3;
params_known.rho=rho;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Here are the structural parameters%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The constant in the flow profit:
theta_0 = params.theta_0(params.indications_num==indication_num_use,:);
%The coefficient in front of the DALY (market size) variable
theta_r = params.theta_r(params.indications_num==indication_num_use,:);
%The competition effect:
theta_s = params.theta_s(params.indications_num==indication_num_use,:);

%Stage 3 flow costs:
flow_cost_3 = params.c_3(params.indications_num==indication_num_use,:);
%Stage 2 flow costs:
flow_cost_2 = params.c_2(params.indications_num==indication_num_use,:);

%Probability of receiving a good signal after stage 1:
p_1=params.p_1(params.indications_num==indication_num_use,:);
%Probability of receiving a good signal after stage 2:
p_2=params.p_2(params.indications_num==indication_num_use,:);

%Upload the data that has the states that appeared along the simulated path%
folder = "C:\Users\ekhmelni\Dropbox\Seeds and Drugs\Pharma Projects\KatyaJMP\DrugAttrition_Model\Estimation\resultsother";
filename=fullfile(folder,"data_main_generated_4_noop_2_params_tmp.csv");
data=readtable(filename);

%Separating across decision stages (decisions are made only after stages 1
%and 2)
inds.for1=data.state_index((data.action==11 | data.action==10) & data.indication_num==indication_num_use);
inds.for2=data.state_index((data.action==21 | data.action==20) & data.indication_num==indication_num_use);
actions.for1=data.action((data.action==11 | data.action==10) & data.indication_num==indication_num_use)-10;
actions.for2=data.action((data.action==21 | data.action==20) & data.indication_num==indication_num_use)-20;
%Minus 10 and minus 20 just translates it into ones and zeroes

%FIRST ATTEMPT - WITH THE THEORETICAL CTPS

    %Upload the TRUE THEORETICAL CCPs 
    folder = "C:\Users\ekhmelni\Dropbox\Seeds and Drugs\Pharma Projects\KatyaJMP\DrugAttrition_Model\Estimation\matfiles";
    filename=fullfile(folder,"ccps_orig_allinds_3_tmp");
    ccps = load(filename);
    ccps=ccps.ccps_allinds;
    ccps=ccps.("ind_"+string(indication_num_use));
    
    %Get the corresponding CTPs
    %THESE ARE CCPs FROM
    ctps = ccps;
    ctps(:,1)=params.p_1(params.indications_num==indication_num_use,:)*ccps(:,1);
    ctps(:,2)=params.p_2(params.indications_num==indication_num_use,:)*ccps(:,2);



%Parameters for the VF approximation:
%The number of basis functions:
n_basis = 200;
%The variance for the nurmal distribution from which the basis function
%parameters are drawn:
var_basis = 0.01;



%Getting the state data for Depression:
data_state=[data.n_in_1, data.n_in_2, data.n_in_3, data.n_in_4];
%Here for the states_to_use I will just drop the states where there is zero
%firms in the corresponding stage
states_to_use.for1=data_state((data.action==11 | data.action==10) & data.indication_num==indication_num_use,:);
states_to_use.for2=data_state((data.action==21 | data.action==20) & data.indication_num==indication_num_use,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Creating the adjacent states%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%1.Creating the adjacent states associated with my own trantision (for
%%%transitions from 1 to 2 and transitions from 2 to 3)
states_once_adjacent_i = get_states_once_adjacent_i_adp(states_to_use,max_firms,ctp_stage_name);

%%%2. Creating the adjacent-adjacent states associated with competitor's transitions (that will
%%%follow my own transition from 1 to 2 or from 2 to 3)
states_twice_adjacent_j = get_states_twice_adjacent_j_adp(states_once_adjacent_i,max_firms,state_transition_names,stage_from,stage_to);

%%%3. Creating adjacent-adjacent states associated with my own further
%%%transition from stage 2 to stage 3 (from states_once_adjacent_i.for1)
%%%or from stage 3 to stage 4 (from states_once_adjacent_i.for1)
states_twice_adjacent_i = get_states_twice_adjacent_i_adp(states_once_adjacent_i,max_firms,ctp_stage_name);

%%%4. Creating adjacent-adjacent-adjacent state associated with competitor's transitions (that will
%%%follow my own transition from 3 to 4)
states_thrice_adjacent_j = get_states_thrice_adjacent_j_adp(states_twice_adjacent_i,max_firms,state_transition_names,stage_from,stage_to);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Creating the associated transition probabilities%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%1.Creating the transition probabilities for competitor's transitions (that will
%%%follow my own transition from 1 to 2 or from 2 to 3)
ctps_once_adjacent_i = get_ctps_once_adjacent_i_adp_theory(states_once_adjacent_i,supportS,ctps,ctp_stage_name);

%%%1.Creating the transition probabilities for competitor's transitions (that will
%%%follow my own transition from 3 to 4)
ctps_twice_adjacent_i = get_ctps_twice_adjacent_i_adp_theory(states_twice_adjacent_i,supportS,ctps);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Creating the coefficients for the right-hand side of the Bellmans%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Here I create coefficients in front of the value functions associated 
% with potential competitors' entries/transitions/terminations for the
%right-hand side of the Bellman equations
%Coefficients include:
    %1. Rates lambda - estimated in the first step
    %2. CTPs - estimated in the first step
    %3. State vector - known
%I also create the rates vector that is used to discount the right-hand
%side of the Bellman equation (those are stored in the .sum part of the 
% appropriate structure)

%1. Coefficients following firm i's own transition from 1 to 2 or from 2 to 3
%get_coefficients_base_adp creates coefficients in the same order as%%%
%state_transition_names={'from0to1','from0to2','from0to3','from0to4','from1to2','from1to0','from2to3','from2to0','from3to4','from3to0','from4to0'};
for origin_stage=1:2
    current_stage=origin_stage+1;
    states=states_once_adjacent_i.(ctp_stage_name{origin_stage}).states;
    ctps=ctps_once_adjacent_i.(ctp_stage_name{origin_stage});
    coefficients_once_adjacent_i.(ctp_stage_name{origin_stage}) = get_coefficients_base_adp(states,current_stage,ctps,params_known,entry_lambdas,lambdas,max_firms,state_transition_names);
    %Add the discount factor:
    coefficients_once_adjacent_i.(ctp_stage_name{origin_stage}).sum = coefficients_once_adjacent_i.(ctp_stage_name{origin_stage}).sum + params_known.rho;
end

%2. Coefficients following firm i's own transition from 3 to 4
coefficients_twice_adjacent_i.for2 = get_coefficients_base_adp(states_twice_adjacent_i.for2.states,4,ctps_twice_adjacent_i.for2,params_known,entry_lambdas,lambdas,max_firms,state_transition_names);
%Add the discount factor (augmented due to the probability of generic entry
%as effectively they are equivalent):
coefficients_twice_adjacent_i.for2.sum=coefficients_twice_adjacent_i.for2.sum + lambdas.lambda_4_e + params_known.rho;


%%%%%%%%%%%%%%%%%%%%
%%%Start with V_4%%%
%%%%%%%%%%%%%%%%%%%%

%1. Create the n_basis-by-4 coefficients for the sigmoid function for each
%s_1, s_2, s_3, s_4
%Those will be fixed for V-4
%coefs_basis_4 = normrnd(0,var_basis,[n_basis,4]);
%Initial guess for the weights for V_4:
%weights_4=normrnd(0,1,[1,n_basis]);
%The left-hand side of the Bellman, approximated given the coefficients and
%the weights:
% V_4_lhs=zeros(size(states_to_use.for4,1),1);
% for i=1:size(states_to_use.for4,1)
%     V_4_lhs(i,:)=weights_4*(1./(1+exp(-coefs_basis_4*transpose(states_to_use.for4(i,:)))));
% end

%1. Create the the 4-th order polynomial associated with the state
%variables s_1, s_2, s_3, s_4 for the LEFT-HAND side
%The order in which the resulting variables are stored is in
%create_polyvars(states_twice_adjacent_i.for2.states).varnames;
states_twice_adjacent_i_for2_polyvars=create_polyvars(states_twice_adjacent_i.for2.states).vars;
polynomial_varorder=create_polyvars(states_twice_adjacent_i.for2.states).varnames;

%1. Create the the 4-th order polynomials associated with the state
%variables s_1, s_2, s_3, s_4 for the RIGHT-HAND side
for j = 1:numel(state_transition_names)
    state=states_thrice_adjacent_j.(state_transition_names{j}).states;
    states_thrice_adjacent_j_polyvars.(state_transition_names{j}) = create_polyvars(state).vars;
end


%%%%%%WEIGHTS HERE WILL CHANGE%%%%%%%%%%%%%%%%%%%%

%Initial guess for the weights for V_4:
weights_4=normrnd(0,1,[size(states_twice_adjacent_i_for2_polyvars,2),1]);

%The LEFT-HAND side of the Bellman, approximated given polynomial order
%and the weights:
V_4_lh=states_twice_adjacent_i_for2_polyvars*weights_4;

%The V-s from the RIGHT-HAND side of the Bellman, approximated given polynomial order
%and the weights:
V_4_rh_j_separate=zeros(size(states_twice_adjacent_i_for2_polyvars,1),numel(state_transition_names));
for j = 1:numel(state_transition_names)
    V_4_rh_j_separate(:,j)=states_thrice_adjacent_j_polyvars.(state_transition_names{j})*weights_4;
end

%The full RIGHT-HAND side of the Bellman:

V_4_rh_j_total=sum( V_4_rh_j_separate.*coefficients_twice_adjacent_i.for2.separate,2);
r_4=profits_5(states_twice_adjacent_i.for2.states(:,4),marketsize,theta_0,theta_r,theta_s);
%While in stage 4, the firm receives the flow profit, and at some rate it
%has to exit and get the scrape walue, which is eulergamma in expectation 
V_4_rh_own = r_4 + lambdas.lambda_4.*double(eulergamma);
V_4_rh=(1./(coefficients_twice_adjacent_i.for2.sum)).*(V_4_rh_own+V_4_rh_j_total);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%
%%% V_3 %%%
%%%%%%%%%%%

%1. Create the the 4-th order polynomial for the LEFT-HAND side
states_once_adjacent_i_for2_polyvars=create_polyvars(states_once_adjacent_i.for2.states).vars;

%2. Create the the 4-th order polynomials associated with the state
%variables s_1, s_2, s_3, s_4 for the RIGHT-HAND side
for j = 1:numel(state_transition_names)
    state=states_twice_adjacent_j.for2.(state_transition_names{j}).states;
    states_twice_adjacent_j_for2_polyvars.(state_transition_names{j}) = create_polyvars(state).vars;
end

%%%%%%WEIGHTS HERE WILL CHANGE%%%%%%%%%%%%%%%%%%%%

%Initial guess for the weights for V_3:
weights_3=normrnd(0,1,[size(states_once_adjacent_i_for2_polyvars,2),1]);

%The LEFT-HAND side of the Bellman, approximated given polynomial order
%and the weights:
V_3_lh=states_once_adjacent_i_for2_polyvars*weights_3;

%The V-s from the RIGHT-HAND side of the Bellman, approximated given polynomial order
%and the weights:
V_3_rh_j_separate=zeros(size(states_once_adjacent_i_for2_polyvars,1),numel(state_transition_names));
for j = 1:numel(state_transition_names)
    V_3_rh_j_separate(:,j)=states_twice_adjacent_j_for2_polyvars.(state_transition_names{j})*weights_3;
end

V_3_rh_j_total=sum( V_3_rh_j_separate.*coefficients_once_adjacent_i.for2.separate,2);
%While in stage 3, the firm pays the flow cost, and at some rate it
%has to either transition to the next stage, or to exit and get the scrape walue, which is eulergamma in expectation 
%Notice that if the next stage is full, it has to exit (I increase the
%limits on the number of firms in each stage to ensure that that never
%happens)
V_3_rh_own=flow_cost_3+(lambdas.lambda_3).*(params_known.p_3).*(states_twice_adjacent_i.for2.not_at_constraint).*V_4_lh+(lambdas.lambda_3).*(1-(params_known.p_3).*(states_twice_adjacent_i.for2.not_at_constraint)).*double(eulergamma);
V_3_rh=(1./(coefficients_once_adjacent_i.for2.sum)).*(V_3_rh_own+V_3_rh_j_total);




%%%%%%%%%%%
%%% V_2 %%%
%%%%%%%%%%%

%1. Create the the 4-th order polynomial for the LEFT-HAND side
states_once_adjacent_i_for1_polyvars=create_polyvars(states_once_adjacent_i.for1.states).vars;

%2. Create the the 4-th order polynomials associated with the state
%variables s_1, s_2, s_3, s_4 for the RIGHT-HAND side
for j = 1:numel(state_transition_names)
    state=states_twice_adjacent_j.for1.(state_transition_names{j}).states;
    states_twice_adjacent_j_for1_polyvars.(state_transition_names{j}) = create_polyvars(state).vars;
end

%3. Create the the 4-th order polynomials associated with the transition to
%the next stage, stage 3 (so for V_3) on the RIGHT-HAND side

states_twice_adjacent_i_for1_polyvars=create_polyvars(states_twice_adjacent_i.for1.states).vars;

%%%%%%WEIGHTS HERE WILL CHANGE%%%%%%%%%%%%%%%%%%%%

%Initial guess for the weights for V_3:
weights_2=normrnd(0,1,[size(states_once_adjacent_i_for1_polyvars,2),1]);

%The LEFT-HAND side of the Bellman, approximated given polynomial order
%and the weights:
V_2_lh=states_once_adjacent_i_for1_polyvars*weights_2;

%The V-s from the RIGHT-HAND side of the Bellman, approximated given polynomial order
%and the weights:
V_2_rh_j_separate=zeros(size(states_once_adjacent_i_for1_polyvars,1),numel(state_transition_names));
for j = 1:numel(state_transition_names)
    V_2_rh_j_separate(:,j)=states_twice_adjacent_j_for1_polyvars.(state_transition_names{j})*weights_2;
end

V_2_rh_j_total=sum( V_2_rh_j_separate.*coefficients_once_adjacent_i.for1.separate,2);
%While in stage 2, the firm pays the flow cost, and at some rate it
%has to either transition to the next stage, or to exit and get the scrape walue, which is eulergamma in expectation 
%Notice that if the next stage is full, it has to exit (I increase the
%limits on the number of firms in each stage to ensure that that never
%happens)
V_3_rh_i_stage2=states_twice_adjacent_i_for1_polyvars*weights_3;
expect_max_2=log(1+exp(V_3_rh_i_stage2))+double(eulergamma);
V_2_rh_own=flow_cost_2+(lambdas.lambda_2).*p_2.*(states_twice_adjacent_i.for1.not_at_constraint).*expect_max_2+(lambdas.lambda_2).*(1-p_2.*(states_twice_adjacent_i.for1.not_at_constraint)).*double(eulergamma);
V_2_rh=(1./(coefficients_once_adjacent_i.for1.sum)).*(V_2_rh_own+V_2_rh_j_total);


%%%Next step - write it as a function with respect to the weights%%%



