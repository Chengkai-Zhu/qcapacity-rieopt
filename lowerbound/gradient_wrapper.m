function grad_struct = gradient_wrapper(X_struct, n, d, dR, KADncopy)

    % X_3D = cat(3, X_struct.R1, X_struct.Main, X_struct.R2);

    [gradR1, gradR2, grad] = cohinfo_grad_localU(X_struct, n, d, dR, KADncopy);

    grad_struct.R1   = gradR1;
    grad_struct.Main = grad(:,:,2:end-1);
    grad_struct.R2   = gradR2;
end