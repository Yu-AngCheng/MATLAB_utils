clear all;
close all;
%%
% We dot not show the influence of guess rate and starting point.
%%
P_guess = 0;
nUps = 2; % 2-up-1-down
trial_list = 5:5:500;
Bias = NaN(max(trial_list),length(trial_list));
Variance = NaN(max(trial_list),length(trial_list));
Error = NaN(max(trial_list),length(trial_list));
tic
for j = 1:length(trial_list)
    [Bias(1:trial_list(j),j), Variance(1:trial_list(j),j), Error(1:trial_list(j),j)] = staircase_simulation_fixed_trials(trial_list(j), P_guess, nUps);
end
toc
save simulated_fixed_trials_optimal_2_down_1_up_ver2.mat
%%
P_guess = 0;
nUps = 3; % 3-up-1-down
trial_list = 5:5:500;
Bias = NaN(max(trial_list),length(trial_list));
Variance = NaN(max(trial_list),length(trial_list));
Error = NaN(max(trial_list),length(trial_list));
tic
for j = 1:length(trial_list)
    [Bias(1:trial_list(j),j), Variance(1:trial_list(j),j), Error(1:trial_list(j),j)] = staircase_simulation_fixed_trials(trial_list(j), P_guess, nUps);
end
toc
save simulated_fixed_trials_optimal_3_down_1_up_ver2.mat
%%
for j = 1:length(trial_list)
    [Error_optimal(j),idx] = min(Error(:,j),[],'omitnan');
    Bias_optimal(j) = Bias(idx,j);
    Variance_optimal(j) = Variance(idx,j);
    temp = 1:1:trial_list(j);
    Cal_Trials(j) = temp(idx);
end
f = figure();
f.Position = [200,200,1000,400];
subplot(2,2,1);
plot(trial_list,Bias_optimal,LineWidth=1,Color='k');
ylabel("Bias (log units)",FontSize=12,FontWeight="bold");
xlabel("# of trials",FontSize=12,FontWeight="bold");
box off;
subplot(2,2,2);
plot(trial_list,Variance_optimal,LineWidth=1,Color='k');
ylabel("Variance (log units)",FontSize=12,FontWeight="bold")
box off;
subplot(2,2,3);
plot(trial_list,Error_optimal,LineWidth=1,Color='k');
ylabel("Error (log units)",FontSize=12,FontWeight="bold")
box off;
subplot(2,2,4);
plot(trial_list,trial_list-Cal_Trials,LineWidth=1,Color='k');
ylabel("Throw Trials",FontSize=12,FontWeight="bold")
box off;
%%
function [Bias, Variance, Error] = staircase_simulation_fixed_trials(N_trials, P_guess, nUps)
underlying_threshold = 1;
N_simulation = 100000;
final_threshold_list = [];
parfor i = 1:N_simulation
    final_threshold_temp = generate_simulation(underlying_threshold, N_trials, P_guess,nUps);
    final_threshold_list = [final_threshold_list; final_threshold_temp];
end
Bias = mean(final_threshold_list, 1); % in log units;
Variance = var(final_threshold_list, [], 1); % in log units;
Error = mean(final_threshold_list.^2, 1); % in log units;
end

%%
function final_threshold = generate_simulation(underlying_threshold, N_trials, P_guess, nUps)
starting_point = 2 * underlying_threshold;
convergence_rate = 0.5^(1/nUps);
SNR = 2 * norminv(convergence_rate);
internal_noise = underlying_threshold/SNR;
current_stimuli = starting_point;
history.testValue = [];
history.isReversal = 0;
history.correct = [];
history.nUp = 0;
history.UpOrDown = [];
for i_trial = 1:N_trials
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
end
for N_cal = 1:N_trials
    final_threshold(N_cal) = log10(geomean(history.testValue(end-N_cal+1:end)));
end
end
%%
function response = generate_response(stimuli, internal_noise, guess_rate)
ACC = normcdf((stimuli/2)/internal_noise)*(1-guess_rate) + 1/2*guess_rate;
response = binornd(1,ACC); % correct or not
end