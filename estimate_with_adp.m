clear;

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

%Number of diseases:
n_indications=size(params.indications_num,1);
%I will run the value function approximation for Depression - indication 9%
indication_num_use=9;
%Max firms for each stage:
max_firms_all=[params.max_in_1,params.max_in_2,params.max_in_3,params.max_in_4];
max_firms=max_firms_all(params.indications_num==indication_num_use,:);
%Market size
marketsize = params.marketsize(params.indications_num==indication_num_use,:);

%%%Known parameters (estimated at step 1)%%%

%Stage durations:
lambdas.lambda_1=1/params.stageduration_1_years(params.indications_num==indication_num_use,:);
lambdas.lambda_2=1/params.stageduration_2_years(params.indications_num==indication_num_use,:);
lambdas.lambda_3=1/params.stageduration_3_years(params.indications_num==indication_num_use,:);
lambdas.lambda_4=1/params.stageduration_4_years(params.indications_num==indication_num_use,:);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Here are the structural parameters%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The constant in the flow profit:
theta_0 = params.theta_0(params.indications_num==indication_num_use,:);
%The coefficient in front of the DALY (market size) variable
theta_r = params.theta_r(params.indications_num==indication_num_use,:);
%The competition effect:
theta_s = params.theta_s(params.indications_num==indication_num_use,:);



%Upload the data that has the states that appeared along the simulated path%
folder = "C:\Users\ekhmelni\Dropbox\Seeds and Drugs\Pharma Projects\KatyaJMP\DrugAttrition_Model\Estimation\resultsother";
filename=fullfile(folder,"data_main_generated_4_noop_2_params.csv");
data=readtable(filename);

%Parameters for the VF approximation:
%The number of basis functions:
n_basis = 200;
%The variance for the nurmal distribution from which the basis function
%parameters are drawn:
var_basis = 0.01;

%Some supplemental arrays
state_transition_names={'from0to1','from0to2','from0to3','from0to4','from1to2','from1to0','from2to3','from2to0','from3to4','from3to0','from4to0'};
stage_from={0,0,0,0,1,1,2,2,3,3,4}';
stage_to={1,2,3,4,2,0,3,0,4,0,0}';
ctp_stage_name={'for1','for2','for3','for4'};

%Getting the state data for Depression:
data_state=[data.n_in_1, data.n_in_2, data.n_in_3, data.n_in_4];
%Here for the states_to_use I will just drop the states where there is zero
%firms in the corresponding stage
states_to_use.for1=data_state(data_state(:,1)>0 & data.indication_num==indication_num_use,:);
states_to_use.for2=data_state(data_state(:,2)>0 & data.indication_num==indication_num_use,:);
states_to_use.for3=data_state(data_state(:,3)>0 & data.indication_num==indication_num_use,:);
states_to_use.for4=data_state(data_state(:,4)>0 & data.indication_num==indication_num_use,:);

%%%1.Creating the adjacent states associated with my own trantision
%Notice that by construction of max_firms, one transition ahead
%should not put me at the boundary
states_adjacent_i = get_states_adjacent_i_adp(states_to_use,ctp_stage_name);

%%%2. Creating states associated with competitor's transitions
states_adjacent_j = get_states_adjacent_j_apd(states_to_use,max_firms,ctp_stage_name,state_transition_names,stage_from,stage_to);

%%%%%%%%%%%%
%%%Checks%%%
%%%%%%%%%%%%

%Total rate of transitioning out of state (taking into account that at the
%max some transitions cannot happen)
not_at_max.for1=states_to_use.for1<max_firms;
not_at_max.for2=states_to_use.for2<max_firms;
not_at_max.for3=states_to_use.for3<max_firms;
not_at_max.for4=states_to_use.for4<max_firms;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Combining the rates (the full discount factor)%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rates_vector.for4=entry_lambdas.lambda_1*not_at_max.for4(:,1)+entry_lambdas.lambda_2*not_at_max.for4(:,2)+entry_lambdas.lambda_3*not_at_max.for4(:,3)+entry_lambdas.lambda_4*not_at_max.for4(:,4)+...
            +lambdas.lambda_1*states_to_use.for4(:,1)+lambdas.lambda_2*states_to_use.for4(:,2)+lambdas.lambda_3*states_to_use.for4(:,3)+lambdas.lambda_4*states_to_use.for4(:,4);

discount_rates.for4 = rates_vector.for4 + rho + 1/years_till_generic;

rates_vector.for3=entry_lambdas.lambda_1*not_at_max.for3(:,1)+entry_lambdas.lambda_2*not_at_max.for3(:,2)+entry_lambdas.lambda_3*not_at_max.for3(:,3)+entry_lambdas.lambda_4*not_at_max.for3(:,4)+...
            +lambdas.lambda_1*states_to_use.for3(:,1)+lambdas.lambda_2*states_to_use.for3(:,2)+lambdas.lambda_3*states_to_use.for3(:,3)+lambdas.lambda_4*states_to_use.for3(:,4);

discount_rates.for3 = rates_vector.for3 + rho;

rates_vector.for2=entry_lambdas.lambda_1*not_at_max.for2(:,1)+entry_lambdas.lambda_2*not_at_max.for2(:,2)+entry_lambdas.lambda_3*not_at_max.for2(:,3)+entry_lambdas.lambda_4*not_at_max.for2(:,4)+...
            +lambdas.lambda_1*states_to_use.for2(:,1)+lambdas.lambda_2*states_to_use.for2(:,2)+lambdas.lambda_3*states_to_use.for2(:,3)+lambdas.lambda_4*states_to_use.for2(:,4);

discount_rates.for2 = rates_vector.for2 + rho;





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
polyvars_4_lh=create_polyvars(states_to_use.for4).vars;

%1. Create the the 4-th order polynomials associated with the state
%variables s_1, s_2, s_3, s_4 for the RIGHT-HAND side

for j = 1:numel(state_transition_names)
    state=states_adjacent_j.for4.(state_transition_names{j});
    states_adjacent_j_to_polyvars.for4.(state_transition_names{j}) = create_polyvars(state).vars;
end




%%%%%%WEIGHTS HERE WILL CHANGE%%%%%%%%%%%%%%%%%%%%

%Initial guess for the weights for V_4:
weights_4=normrnd(0,1,[size(polyvars_4_lh,2),1]);

%The LEFT-HAND side of the Bellman, approximated given polynomial order
%and the weights:
V_4_lh=polyvars_4_lh*weights_4;

%The V-s from the RIGHT-HAND side of the Bellman, approximated given polynomial order
%and the weights:
for j = 1:numel(state_transition_names)
    V_4_rh.(state_transition_names{j})=states_adjacent_j_to_polyvars.for4.(state_transition_names{j})*weights_4;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%













