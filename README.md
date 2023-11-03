## ZetaPersei

Welcome to ZetaPersei. This is a preliminary release of two software tools meant to aid in the planned 2024 commisioning of Keck NIRC2 Polarimetry Mode. This repository is divided into two sections - calibration observation planning and science observation planning

# Science Observation Planning

In this folder, there is an example notebook showing the theoretical observation planning of a science target and how digitized versions of Whittet (1992) and Heiles (2000) with added target information (brightnesses in several bands, cross-referencing with the Washington Double Star Catalog etc) can be used to find suitable unpolarized and polarized standards for a specific date and time of observation

# Calibration Observation Planning

In this folder, using the existing exposure time calculator for NIRC2 provided by Keck Observatory (https://github.com/KeckObservatory/exposureTimeCalculator), an example notebook calculates the SNR of two polarized standards. This SNR is then used to calculate the resultant error on two instrument parameters (Q efficiency factor, and U -> Q conversion factor) necessary to develop a Mueller Matrix model for NIRC2 polarimetry mode. By adjusting the exposure times, number of coadds, and exact standard stars used, a suitable polarized standard star pair and observation plan can be found for a max tolerable error on these two instrument terms. 

More details on the derivations of the error terms on the instrument parameters is included in a PDF in this same folder.