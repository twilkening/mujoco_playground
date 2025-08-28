import matplotlib.pyplot as plt
import pandas as pd

# Read CSV, skip any comment lines
df = pd.read_csv('data.csv', comment='/')
df.columns = df.columns.str.strip()

# Plot 1: Energy (PE, KE, TE) vs Time
plt.figure(figsize=(10, 6))
plt.plot(df['t'], df['PE'], label='Potential Energy (PE)')
plt.plot(df['t'], df['KE'], label='Kinetic Energy (KE)')
plt.plot(df['t'], df['TE'], label='Total Energy (TE)')
plt.xlabel('Time (s)')
plt.ylabel('Energy')
plt.title('Energy vs Time')
plt.legend()
plt.grid(True)

# Plot 2: Joint positions (q1, q2) vs Time
plt.figure(figsize=(10, 6))
plt.plot(df['t'], df['q1'], label='q1')
plt.plot(df['t'], df['q2'], label='q2')
plt.xlabel('Time (s)')
plt.ylabel('Joint Positions')
plt.title('Joint Positions vs Time')
plt.legend()
plt.grid(True)

# Show both plots
plt.show()