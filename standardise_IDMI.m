% WRITTEN BY: ELHUM A SHAMSHIRI (Elhum.Shamshiri@unige.ch)
% CODE SHOULD NOT BE DISTRIBUTED AND ANALYSIS SHOULD NOT BE CONDUCTED WITHOUT PRIOR CONSENT FROM  ELHUM A SHAMSHIRI
% If you would like to use this software for publication please contact

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 

function [IDMI_standardised] = standardise_IDMI(patient_IDMI)

load 'D:\Intrisic_Ignition_Model\mean_controls.mat'
load 'D:\Intrisic_Ignition_Model\std_controls.mat'

% z-score the IDMI values for the patient to resemble the controls
for i = 1:length(patient_IDMI)
    z(1,i) = (patient_IDMI(i) - mean_controls(i))/std_controls(i);
end

IDMI_standardised = z;

end