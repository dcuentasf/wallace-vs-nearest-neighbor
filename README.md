# Melting temperature estimation: Wallace vs Nearest-Neighbor
This project compares two methods for estimating DNA melting temperature (Tm):
- Wallace rule
<img width="1336" height="89" alt="image" src="https://github.com/user-attachments/assets/6618adf8-7749-4c8f-871a-32f44a61aa06" />


- Nearest-neighbor (NN) model + salt correction
  <img width="1369" height="221" alt="image" src="https://github.com/user-attachments/assets/3e4c2c99-6460-4f54-b50f-f2871a7dd44d" />


Simulations were performed in R using randomly generated primers of different lengths.

## Methods
- 3 primer groups: 12–17 bp, 18–24 bp, 25–30 bp
- 1000 simulated primers per group
- Tm calculated using both methods
- Differences (ΔTm and |ΔTm|) analyzed

## Results
- Wallace rule consistently overestimates Tm
- Error increases with primer length
- NN model provides more accurate estimates

## Limitations
- Effect of Mg2+ or dNTPs are not considered
