#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import matplotlib.pyplot as plt

# Read the Excel file
file1_path = '/Users/mozhganoroujlu/Desktop/SA/score/results/modif/generated/all_opcs/150n_e4_acceptance_rates.xlsx'
df1 = pd.read_excel(file1_path)

# Extract T column (assuming it's the first column)
t_values = df1.iloc[:, 0]

# Calculate average acceptance rate across all runs
avg_acceptance_df1 = df1.iloc[:, 1:].mean(axis=1)

# Create the plot
plt.figure(figsize=(10, 6))

# Plot first 5 individual runs
for col in df1.columns[1:6]:  # Run_1 to Run_5
    plt.plot(t_values, df1[col], linewidth=1, alpha=0.6, label=col)

# Plot the average on top, bolder
plt.plot(t_values, avg_acceptance_df1, 'k-', linewidth=2.5, label='Average across 1068 runs', marker='o')

# Customize
plt.xlabel('T', fontsize=12)
plt.ylabel('Acceptance Rate', fontsize=12)
plt.title('Acceptance Rate per Run', fontsize=14)
plt.legend()
plt.xlim(t_values.iloc[0], t_values.iloc[-1])
plt.tight_layout()

# Save high-quality PNG and vector PDF
plt.savefig('acceptance_rate.png', dpi=400, bbox_inches='tight')
plt.savefig('acceptance_rate.pdf', bbox_inches='tight')

plt.show()