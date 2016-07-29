function values=evaluate_buckled_material_prop(mat_prop,time,x,T)

values = mat_prop.ft(time) * mat_prop.fx(x) * mat_prop.fT(T);

return
end