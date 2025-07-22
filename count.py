import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Simulate inversion data
np.random.seed(42)
num_genes = 5000
num_inversions = 30

# Generate random inversion ranges
starts = np.random.randint(1, num_genes - 50, size=num_inversions)
ends = starts + np.random.randint(5, 200, size=num_inversions)
ends = np.clip(ends, starts + 1, num_genes)

inversion_data = pd.DataFrame({'inversion_id': range(1, num_inversions+1),
                               'start_idx': starts,
                               'end_idx': ends})
inversion_data['size'] = inversion_data['end_idx'] - inversion_data['start_idx'] + 1

# Histogram of inversion sizes
plt.figure(figsize=(8, 5))
plt.hist(inversion_data['size'], bins=10, edgecolor='black')
plt.title('Inversion Size Distribution')
plt.xlabel('Number of Genes Involved')
plt.ylabel('Frequency')
plt.grid(True)
plt.tight_layout()
plt.show()

# Show top 10 largest inversions
top10_inversions = inversion_data.sort_values(by='size', ascending=False).head(10)

import ace_tools as tools; tools.display_dataframe_to_user(name="Top 10 Largest Inversions", dataframe=top10_inversions)

# Return inversion data for further use
inversion_data.head()
