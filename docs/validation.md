# Validation against independent data

Each check is a test in `tests/test_validation.py`. Agreement is stated honestly, including what
cannot be compared yet.

| Quantity | Model | Reference | Agreement | Notes |
|----------|-------|-----------|-----------|-------|
| k2 of CO2 + OH⁻ at 25 °C, I = 0 | 8.05e3 L/(mol·s) | Pohorecki & Moniuk (1988) | by construction | the base case (45 °C) is outside their 291–314 K range and is warned |
| Flooding of 25 mm metal Pall rings, air/water, L/V = 1 | F = 3.0 Pa^0.5 | vendor capacity charts (Kister 1992), Billet & Schultes database | within range | |
| Mellapak 250Y effective area, NaOH absorption (74 mm column, 0.4 M, 10 % CO2) | 214 m⁻¹ (Hanley & Chen) | 170 m⁻¹ fitted by Ellebracht et al. (2023) | +26 % | Ellebracht's packing was 3D-printed resin; k_L not comparable (fast reaction, Ha ≈ 11) |
| Cell at 200 mA/cm² | V 1.422 V, t_Na 0.86, 14.6 m² per t/day (1 e⁻) | Zhang et al. (2024) SI Note 8: 1.42 V, 0.87, 14.58 m² | < 1 % | the thesis j–V and t_Na fits come from Zhang's data |

| Initial CO2 capture, Mellapak 250Y in Ellebracht's column (74 mm, 3 × 150 mm, v_G 0.105 m/s, v_L 45 m/h, 10 % CO2, 0.4 M NaOH, 25 °C; their ESI) | 39.5 % | 67 % measured | **under-predicted** | at this laboratory gas load the Hanley & Chen gas-film correlation gives Sh_V ≈ 0.9 (< 2, below pure diffusion), so the gas film wrongly controls; the model now warns. Industrial loads (F ≈ 1–2 Pa^0.5) are inside their data |

## Not yet validated
- The thesis bench (TPMS Schwarz-D, 150 mm column): `data/lab2_titration.csv` (thesis Table 3.1)
  is not in the repository yet, and there is no correlation for TPMS packings (Ellebracht et al.
  give fitted k_L and a_eff at one condition only).
- Liquid-film coefficient of Hanley & Chen for 250Y at low liquid Reynolds number (0.05 mm/s here
  vs 0.23 mm/s fitted by Ellebracht, weakly identified).
