function x0 = barebonesSA(x0,T0,NTemp,Nstep,cooling,moveclass,f, report)
    % pick move class and cooling schedule
    switch moveclass
        case 1
            move = @(state) twobondmove(state);
        case 2
            move = @(state) normalmove(state);
    end
    switch cooling
        case 1
            alpha = 0.99;
            Tstep = @(i) exptemp(T0,alpha,i);
    end

    % optimize
    for i = 1:NTemp
        x0  = metropolis(f,T0,x0,Nstep,move); % update state
        T0  = Tstep(i); % update temp

        if exist('report','var')
            if report == 1
                if mod(i,50) == 0 %report progress in terminal
                    fprintf('Iteration no. %i of %i finished \n',i,NTemp)
                end
            end
        end
    end
end
%% moveclass, cooling schedule, metropolis
function newtemp = exptemp(T0,alpha,i)
    newtemp = T0 * alpha^i;
end

function newstate = twobondmove(state)
    i = randi(length(state),1);
    newstate = circshift(state',i-1)';
    j = 2 + randi(length(state)-2);
    newstate = [newstate(1), fliplr(newstate(2:j-1)), newstate(j:end)];
end

function newstate = normalmove(state)
    dist = randn(size(state));
    newstate = state + dist;
end

function state = metropolis(f,T,istate,nsteps,move)
    for i = 1:nsteps
        nstate = move(istate);
        df = f(nstate) - f(istate);
        if df <= 0
            state = nstate;
        else
            R = rand; %uniform between [0,1]
            if R < exp(-df/T)
                state = nstate;
            else
                state = istate;
            end
        end
    end
end