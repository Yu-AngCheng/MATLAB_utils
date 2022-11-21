%%
% We dot not show the influence of guess rate and starting point.
clear all;
close all;
%%
P_guess = 0;
nUps = 2; % 2-up-1-down
reversal_list = 1:100;
Bias = NaN(max(reversal_list),length(reversal_list));
Variance = NaN(max(reversal_list),length(reversal_list));
Error = NaN(max(reversal_list),length(reversal_list));
Total_trials = NaN(length(reversal_list));
tic
for j = 1:length(reversal_list)
    [Bias(1:reversal_list(j),j), Variance(1:reversal_list(j),j), ...
     Error(1:reversal_list(j),j), Total_trials(j)] = ...
     staircase_simulation_fixed_reversals(reversal_list(j), P_guess, nUps);
end
toc
save simulated_fixed_reversals_optimal_2_down_1_up_ver2.mat
%%
P_guess = 0;
nUps = 3; % 3-up-1-down
reversal_list = 1:100;
Bias = NaN(max(reversal_list),length(reversal_list));
Variance = NaN(max(reversal_list),length(reversal_list));
Error = NaN(max(reversal_list),length(reversal_list));
Total_trials = NaN(length(reversal_list));
tic
for j = 1:length(reversal_list)
    [Bias(1:reversal_list(j),j), Variance(1:reversal_list(j),j), ...
     Error(1:reversal_list(j),j), Total_trials(j)] = ...
     staircase_simulation_fixed_reversals(reversal_list(j), P_guess, nUps);
end
toc
save simulated_fixed_reversals_optimal_3_down_1_up_ver2.mat
%%
for j = 1:length(reversal_list)
    [Error_optimal(j),idx] = min(Error(:,j),[],'omitnan');
    Bias_optimal(j) = Bias(idx,j);
    Variance_optimal(j) = Variance(idx,j);
    temp = 1:1:reversal_list(j);
    Cal_reversals(j) = temp(idx);
end
f = figure();
f.Position = [200,200,1200,600];
subplot(2,3,1);
plot(reversal_list,Bias_optimal,LineWidth=1,Color='k');
ylabel("Bias (log units)",FontSize=12,FontWeight="bold");
xlabel("# of trials",FontSize=12,FontWeight="bold");
box off;
subplot(2,3,2);
plot(reversal_list,Variance_optimal,LineWidth=1,Color='k');
ylabel("Variance (log units)",FontSize=12,FontWeight="bold")
box off;
subplot(2,3,3);
plot(reversal_list,Error_optimal,LineWidth=1,Color='k');
ylabel("Error (log units)",FontSize=12,FontWeight="bold")
box off;
subplot(2,3,4);
plot(reversal_list,Total_trials,LineWidth=1,Color='k');
ylabel("Total Trials",FontSize=12,FontWeight="bold")
box off;
subplot(2,3,5);
plot(reversal_list,reversal_list-Cal_reversals,LineWidth=1,Color='k');
ylabel("Throw Reversals",FontSize=12,FontWeight="bold")
box off
%%
function [Bias, Variance, Error, Total_trials] = staircase_simulation_fixed_reversals(N_reversals, P_guess, nUps)
underlying_threshold = 1;
N_simulation = 100000;
final_threshold_list = [];
total_trials_list = [];
parfor i = 1:N_simulation
    [final_threshold_temp, total_trials_temp] = generate_simulation(underlying_threshold,N_reversals, P_guess,nUps);
    final_threshold_list = [final_threshold_list; final_threshold_temp];
    total_trials_list = [total_trials_list, total_trials_temp];
end
Bias = mean(final_threshold_list, 1); % in log units;
Variance = var(final_threshold_list, [], 1); % in log units;
Error = mean(final_threshold_list.^2, 1); % in log units;
Total_trials = mean(total_trials_list);
end

%%
function [final_threshold, total_trials] = generate_simulation(underlying_threshold, N_reversals, P_guess, nUps)

    starting_point = 1.5 * underlying_threshold;
    convergence_rate = 0.5^(1/nUps);
    SNR = 2 * norminv(convergence_rate);
    internal_noise = underlying_threshold/SNR;
    current_stimuli = starting_point;
    history.testValue = [];
    history.isReversal = 0;
    history.correct = [];
    history.nUp = 0;
    history.UpOrDown = [];

    i_trial = 0;
    i_reversal = 0;
    while i_reversal < N_reversals % Not enough reversals are got
        i_trial = i_trial + 1;
        history.testValue = [history.testValue, current_stimuli];
        response_temp = generate_response(current_stimuli,internal_noise,P_guess);
        history.correct=[history.correct response_temp];
        if response_temp == 1
            if history.nUp(i_trial) + 1 >= nUps
                history.nUp = [history.nUp,0];
                history.UpOrDown=[history.UpOrDown -1];
                current_stimuli = current_stimuli/(10^0.05);
            else
                history.nUp=[history.nUp history.nUp(i_trial)+1];
                current_stimuli = current_stimuli;
                if i_trial == 1
                    history.UpOrDown=[history.UpOrDown,0];
                else
                    history.UpOrDown=[history.UpOrDown,history.UpOrDown(i_trial-1)];
                end
            end

        elseif response_temp == 0
            history.nUp = [history.nUp 0];
            current_stimuli = current_stimuli*(10^0.05);
            history.UpOrDown=[history.UpOrDown,1];
        end

        if  i_trial ~= 1
            if (history.UpOrDown(i_trial-1)+history.UpOrDown(i_trial)==0)&&(history.UpOrDown(i_trial)~=0)
                history.isReversal=[history.isReversal 1];
            else
                history.isReversal=[history.isReversal 0];
            end
        end
        i_reversal = sum(history.isReversal);
    end
    ReversalIndex = find(history.isReversal == 1);
    ReversalValue = history.testValue(ReversalIndex);
    for N_cal = 1:N_reversals
        ReversalCal_temp = ReversalValue(end-N_cal+1:end);
        final_threshold(N_cal) = log10(geomean(ReversalCal_temp));
    end
    total_trials = length(history.correct);

end
%%
function response = generate_response(stimuli, internal_noise, guess_rate)
ACC = normcdf((stimuli/2)/internal_noise)*(1-guess_rate) + 1/2*guess_rate;
response = binornd(1,ACC); % correct or not
end