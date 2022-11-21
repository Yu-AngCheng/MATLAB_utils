clear all
%% Staircase Parameters
pStaircase.nUps = 2;                     % n-up-1-down, for n: 3-> 79.6%  2->7.70%  1->~50% 
pStaircase.nChanges = 2;                 % Num. of reversals after which the step size changes to 1.
pStaircase.initStepSize = 2;             % Initial step size before nChanges, just to speed up convergence
pStaircase.nPractice = 3;                % Num. of practice trials.
pStaircase.nReversals = 10;              % Num. of reversals to end the staircase
pStaircase.nCal = 6;                     % Num. of reversals used for computing final threshold
pStaircase.initValue = 10;               % Intial stimuli strength, but may be adjusted to the nearest test condition
pStaircase.Scale = 1;                    % This identifies the type of the scale of conditions vector, 0 for linear; 1 for logarithm scale.
pStaircase.testCondition = 10.^(0.5:0.05:1.2);  % The test condition, here is 0.05 log unit from 0.01 to 1
pStaircase.StepSize = 10^0.05;           % For log scale, it should be the ratio; For linear scale, it should be difference

%% Record the results of current data
initValueIndex = find(pStaircase.testCondition >= pStaircase.initValue); % We restrict the testing conditions in pStaircase.testCondition
history.testValue = pStaircase.testCondition(initValueIndex(1));         % Start with the smallest value in test conditions that is greater than initial value 
history.isReversal = 0;       % Whether is a reversal
history.correct = [];         % Correct or not in this trial
history.nUp = 0;              % How many trials is accumulated to be correct for the current stimuli value
history.UpOrDown = [];        % The trend of the psychometrics is up or down, only to calculate the reversals
i_trial = 0;
i_reversal = 0;
while i_reversal < pStaircase.nReversals % Not enough reversals are got
    i_trial = i_trial + 1;
    currentValue = history.testValue(i_trial);
    history.correct = [history.correct, generate_response(currentValue)];
    history = staircaseUpdate(history, pStaircase, i_trial);
    i_reversal = sum(history.isReversal);
end

%% Calculating threshold based on reversals
ReversalIndex = find(history.isReversal == 1);
ReversalValue = history.testValue(ReversalIndex);
ReversalCal = ReversalValue(end-pStaircase.nCal+1:end);
threshold = geomean(ReversalCal);

%%
function response = generate_response(stimuli)
    response = rand>0.5; % correct or not
end