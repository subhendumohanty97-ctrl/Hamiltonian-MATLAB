% Assume Ir and Ia are already defined from:
for k = 1:100
    norm_factor = sqrt(evecr(:,k)' * A * evecr(:,k));
    evecr(:,k) = evecr(:,k) / norm_factor;
end
for k = 1:100
    norm_factor = sqrt(eveca(:,k)' * A * eveca(:,k));
    eveca(:,k) = eveca(:,k) / norm_factor;
end

Ir = evecr' * A * evecr;
Ia = eveca' * A * eveca;

error_thresh = [1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1];
nnz_vals = zeros(size(error_thresh));  % To store counts
nnz_vals1 = zeros(size(error_thresh));
for i = 1:length(error_thresh)
    threshold = error_thresh(i);  % DO NOT overwrite variable 'error'
    nnz_vals(i) = nnz(abs(Ir) >= threshold);  % Count how many elements differ by >= threshold
    nnz_vals1(i) = nnz(abs(Ia) >= threshold);
end

% Plotting
figure;
plot(log10(error_thresh), nnz_vals, 'k-s', log10(error_thresh), nnz_vals1, 'k-o','LineWidth',1.5)
legend('Proposed','Diag')
%xlabel('Error threshold')
%ylabel('Number of elements greater than threshold')
%title('Sensitivity of Ir vs Ia to threshold')
grid on;
axis square;
