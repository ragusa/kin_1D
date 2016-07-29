function aux=create_buckled_material_prop(time_dep,values,times,space_dep,space_function_handle,gamma,T0)

aux = create_material_prop(time_dep,values,times,space_dep,space_function_handle);
aux.fT = @(T) (1 + gamma*(sqrt(T)-sqrt(T0)));
