# Chemotherapy Analysis Project
The dataset ‘Chemotherapy.csv’ contains the following information:

• Patient: a unique patient id,

• Tumour: level of the tumour marker increase (if positive) or shrinkage (if negative) over a given month of treatment (on the log2 scale), defined as the level of the tumour marker at the end of the month minus its level at the start of the month of interest,

• Line: chemotherapy line. a patient may have several treatments over years, starting with 1st line, then 2nd line a few months later, aso.

• Month: month of treatment for a given line of chemotherapy. Standard chemotherapy treatment last 3 months (per line) but may be shorter/longer depending on the patient health and on the drug combination,

• Sensitivity: each patient had a biopsy of the tumour after diagnosis. The sensitiv- ity of the tumour to each drug combination was assessed in-vitro. The score is high (around 1 and above) for tumour samples which are resistant to the drug combina- tion of interest and low (around 0) for tumour samples which are sensitive to the combination of drug of interest.

It is believed that the in-vitro assessment of the sensitivity of the tumour to a drug combination is a reliable predictor of the treatment effect of the same drug combination on the patient. However, the sensitivity score is believe to best work for initial chemotherapy lines (as the tumour changes afterwards and become resistant to the initial drug combination).
You have been asked by the group leader of a research group at CRUK to check

• if the in-vitro sensitivity score can be used to predict the treatment effect on the patient at the hospital,

• if the in-vitro sensitivity score decreases with time and is stronger for line 1 than for lines 2 and 3+ (you can consider the variable "line" as a 3-level factor with levels "1st line", "2nd line" and "3rd line and more").
