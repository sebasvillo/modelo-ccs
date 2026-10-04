# Data

## `lab2_titration.csv`: thesis bench runs (thesis Table 3.1)

Raw measurements only, copied from the thesis workbook `Resultados Lab(4)_con_calculos.xlsx`
(sha1 721a192f…), sheet `Hoja1`, rows 3–14. Bench: packed column with a 3D-printed TPMS
(Schwarz-D) packing, 15 or 35 cm bed, gas = CO2 + N2 at 25 °C and 1 atm (flowmeters at 25 °C),
solvent NaOH or KOH. Each run: a solvent aliquot (5 mL, diluted with water when noted) titrated
with HCl to the phenolphthalein end point, three times.

| column | meaning |
|---|---|
| run | test number |
| bed_height_cm | packed height |
| solvent, solvent_conc_M | NaOH or KOH and its nominal concentration before the run |
| aliquot_mL, dilution_water_mL | titrated sample volume and added water |
| HCl_conc_M | titrant concentration |
| liquid_flow_mL_min, CO2_flow_mL_min, N2_flow_mL_min | flows |
| titration1..3_mL | HCl volume to the phenolphthalein end point |

With n_HCl = c_HCl·V: n_CO3 = n_OH,0 − n_HCl and n_OH = 2·n_HCl − n_OH,0 (one H⁺ per OH⁻ and per
CO3²⁻ → HCO3⁻ at that end point).

**Open question before using these data for validation:** the workbook converts the aliquot
carbonate into a capture rate as [CO3]·Q_L, i.e. a single pass of liquid. That gives up to 2.8×
the CO2 fed (run 2), so the solvent was most likely recirculated; the run duration and the
reservoir volume are needed to recover the true capture. The capture percentages typed in
rows 4–14 of the workbook are therefore not used here.
