function  [rhoNext, Pnext, modeLikelihood] = KFaposteriori(rho, P,...
    rhoMeasured, H, rhoJam, percentErrorMeasured)

R = percentErrorMeasured^2*diag(rhoMeasured.^2); % measurement noise
residualCov = H * P * H' + R + 0.001*eye(size(R,1)); % residual cov
gain = P * H' / residualCov; % Kalman gain

residual = rhoMeasured - H * rho; % residual
rhoNext = rho + gain * residual;
rhoNext = max(min(rhoJam,rhoNext),0);

Pnext = (eye(length(rho)) - gain * H) * P;
residualCov = max(residualCov, 0);
residualCov = (residualCov + residualCov')/ 2;
modeLikelihood = mvnpdf(residual, 0, residualCov);