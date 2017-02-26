function y = IQSfun(y,shape_end,shape_beg,X_beg,tn,TR,NFId_new,NFId_old,IV,Re,eps,C_old)

Rep = residual_IQS(shape_end + eps*y,shape_beg,X_beg,tn,TR,NFId_new,NFId_old,IV,C_old);

y = (Rep - Re) / eps;

return
end


