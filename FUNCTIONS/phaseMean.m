function phaseAvg = phaseMean(phaseData)
    phaseAvg = atan2(mean(sin(phaseData)), mean(cos(phaseData)));
end