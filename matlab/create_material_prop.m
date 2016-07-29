function aux=create_material_prop(time_dep,values,times,space_dep,space_function_handle)


if strcmp(time_dep,'constant_in_time')
    aux.ft = @(t) (0.*t+values(1));

elseif strcmp(time_dep,'ramp_in_time')
    % ramp that transitions from val1 to val2 when t is between t1 and t2
    if(abs(times(2)-times(1))<1e-10)
        error('must use different time steps in ramp definition');
    end
    if(times(2)<times(1))
        error('we must have t2>t1 in ramp definition')
    end
    aux.ft = @(t) ( values(1) + (values(2)-values(1))*ramp(t,times(1),times(2)) );

elseif strcmp(time_dep,'ramp2_in_time')
    % ramp that transitions from val1 to val2 when t is between t1 and t2
    if(abs(times(2)-times(1))<1e-10)
        error('must use different time steps in ramp definition');
    end
    if(times(2)<times(1))
        error('we must have t2>t1 in ramp definition')
    end
    aux.ft = @(t) ( values(1) + (values(2)-values(1))*ramp(t,times(1),times(2)) + (values(3)-values(2))*ramp(t,times(3),times(4)) );

else
    error('unknown time dependence');
end

if strcmp(space_dep,'constant_in_space')
    aux.fx =@(x) (0.*x+1);
elseif strcmp(space_dep,'spatial_dep')
    aux.fx = space_function_handle;
else
    error('unknown spatial dependence');
end

aux.fT = @(T) 0.*T+1;

return
end
