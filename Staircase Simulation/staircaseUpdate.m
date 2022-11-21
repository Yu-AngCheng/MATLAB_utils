function history=staircaseUpdate(history, pStaircase, i_trial)

x = pStaircase.Scale;% x==0,linear; x==1,log;
a = pStaircase.StepSize;

if i_trial > pStaircase.nPractice
    if sum(history.isReversal)>= pStaircase.nChanges
        thisStep = 1;
    else
        thisStep = pStaircase.initStepSize;
    end

    if history.correct(i_trial) == 1 % current trial is correct
        if history.nUp(i_trial) + 1 >= pStaircase.nUps % already statisfy n ups
            history.nUp = [history.nUp 0]; % accumulation goes to 0 since the stimuli value will go lower
            nextTestValue = history.testValue(i_trial) / a^(x*thisStep) - a*(1-x)*thisStep;
            history.UpOrDown=[history.UpOrDown -1];
        else % correct but do not statisfy n ups
            history.nUp=[history.nUp history.nUp(i_trial)+1]; % accumulate one more
            nextTestValue = history.testValue(i_trial);
            if i_trial == pStaircase.nPractice + 1
                history.UpOrDown=[history.UpOrDown 0];
            else
                history.UpOrDown=[history.UpOrDown history.UpOrDown(i_trial-1)];
            end
        end

    else % current trial is not correct
        history.nUp = [history.nUp 0]; % accumulation goes to 0 since the stimuli value will go higher
        nextTestValue=history.testValue(i_trial) * a^(x*thisStep) + a*(1-x)*thisStep;
        history.UpOrDown=[history.UpOrDown 1];
    end

    if  i_trial~=pStaircase.nPractice+1
        if (history.UpOrDown(i_trial-1)+history.UpOrDown(i_trial)==0)&&(history.UpOrDown(i_trial)~=0)
            history.isReversal=[history.isReversal 1];
        else
            history.isReversal=[history.isReversal 0];
        end
    end

else % if it is only practice trial
    history.nUp = [history.nUp 0];
    history.UpOrDown = [history.UpOrDown 0];
    history.isReversal = [history.isReversal 0];
    nextTestValue = history.testValue(i_trial);

end

history.testValue=[history.testValue, nextTestValue];
