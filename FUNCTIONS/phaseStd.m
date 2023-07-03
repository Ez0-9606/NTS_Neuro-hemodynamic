function phaseStandDev = phaseStd(phaseData)
    phaseStandDev = (-log(mean(sin(phaseData))^2 + mean(cos(phaseData))^2))^0.5;
end