To run the optimization code, you need the following files in the same directory:

1. optimize_L.m
2. sim_error.m
3. freq_sim_noisy_arm_3state.m
4. sim_data.m

Run optimize_L.m to perform MLE. optimize_L.m defines the starting point of the optimization process based on Linit. The function sim_error.m produces the error signal to minimize over


sim_error calls freq_sim_noisy_arm_3state to simulate a sum-of-sines tracking trajectory.

make_data_3state.m is used to generate synthetic data.
