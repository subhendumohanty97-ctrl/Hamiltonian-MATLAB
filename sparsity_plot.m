% Compute sparsity.
for k = 1:100
    norm_factor = sqrt(evec_prop(:,k)' * M * evec_prop(:,k));
    evec_prop(:,k) = evec_prop(:,k) / norm_factor;
end

for k = 1:100
    norm_factor = sqrt(evec_diag(:,k)' * M * evec_diag(:,k));
    evec_diag(:,k) = evec_diag(:,k) / norm_factor;
end

Ir = evec_prop' * M * evec_prop;
Ia = evec_diag' * M * evec_diag;

error_thresh = [1e-10, 1e-9, 1e-8, 1e-7, 1e-6, ...
                1e-5, 1e-4, 1e-3, 1e-2, 1e-1];

nnz_vals = zeros(size(error_thresh));
nnz_vals1 = zeros(size(error_thresh));

for i = 1:length(error_thresh)
    threshold = error_thresh(i);
    nnz_vals(i) = nnz(abs(Ir) >= threshold);
    nnz_vals1(i) = nnz(abs(Ia) >= threshold);
end

% Plotting.
figure;
plot(log10(error_thresh), nnz_vals, 'k-s', ...
     log10(error_thresh), nnz_vals1, 'k-o', ...
     'LineWidth', 1.5);

legend('Proposed', 'Diag');
grid on;
axis square;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
```
