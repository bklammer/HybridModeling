function [out] = thermoNonideal(in, T, p)
%THERMONONIDEAL Calculates the thermodynamic properties of non-saturated N2O in tank
%   'OF' is the temperature in K
%   'p' is the pressure
%   'in' is a string indicating the desired property
%   'out' is a double giving the desired quantity

    p = p/1000000; % Convert Pa to MPa for lookup table
    if T < N2Ononideal.T(1) % Check that query is inside bounds
        error('massFrac:lowOF', 'Outside temperature range: \n T = %f', T)
    elseif T > N2Ononideal.T(end)
        error('massFrac:highOF', 'Outside temperature range: \n T = %f', T)
    elseif p > N2Ononideal.p(1)
        error('massFrac:lowP', 'Outside pressure range: \n p = %f bar', p)
    elseif p < N2Ononideal.p(end)
        error('massFrac:highP', 'Outside pressure range: \n p = %f bar', p)
    end
    for mm = 1:length(N2Ononideal.T)
        if T < N2Ononideal.T(mm) % Find first index that is larger than T
            break
        end
    end
    for nn = 1:length(N2Ononideal.p)
        if p > N2Ononideal.p(nn) % Find first index that is smaller than p
            break
        end
    end
    if ~any(strcmp(in, fieldnames(N2Ononideal))) % If input is not a field name of N2Ononideal
        error('thermoNonideal:input', 'Input ''%s'' is invalid', in)
    end

    % Do some bilinear interpolation (equations from wikipedia)
    T_int = [N2Ononideal.T(mm)-T, T-N2Ononideal.T(mm-1)];
    p_int = [N2Ononideal.p(nn)-p; p-N2Ononideal.p(nn-1)];
    C_int = 1/((N2Ononideal.T(mm)-N2Ononideal.T(mm-1))*(N2Ononideal.p(nn)-N2Ononideal.p(nn-1)));
    out_int = N2Ononideal.(in)(mm-1:mm,nn-1:nn);
    out = C_int.* T_int*out_int*p_int; % property
end