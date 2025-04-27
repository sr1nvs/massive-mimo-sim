# Massive MIMO System Simulation

This project aims to simulate a massive MIMO system with the following parameters:

- Transmit Antenna: 64-element (8x8) Uniform Rectangular Array
- Receive Antennas: 2 per user
- Number of users: 10
- 128 Subcarriers
- Cell radius: 500m, Min. distance from tx: 35m
- 4-QAM Modulation, and ZF Precoding/Equalization

### Part 1: Power Metrics
- Uses MATLAB's siteviewer tool to visualize the position of BS and users
- Estimates received power at each MS using reflections only (diffraction paths omitted due to time constraints)
- Location: Suburban environment (VIT Chennai)

### Part 2: Performance Metrics
- Model a Rician channel + AWGN (flat fading for each subcarrier)
- Transmitter: Implemented QPSK (4-QAM) modulation and ZF precoding of transmitted data
- Receiver: ZF equalization using channel matrix + demodulation
- Compute error rate (BER/SER) over a range of SNR

## Usage:
- Clone the repository
- Run ``` main.m ```

## Results
- For a cell radius of 500m, the received power at each mobile is -90dBm (avg)
- For all SNR above 0dB, we observe a BER of effectively 0



