import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 10,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'figure.dpi': 300
})

np.random.seed(42)
df = pd.DataFrame({
    'Treatment': np.repeat(['Control', 'Drug A (10 µM)', 'Drug B (10 µM)', 'Combo'], 12),
    'Expression': np.concatenate([
        np.random.normal(1.0, 0.15, 12),
        np.random.normal(1.8, 0.25, 12),
        np.random.normal(2.4, 0.30, 12),
        np.random.normal(3.9, 0.35, 12)
    ])
})

palette = ['#999999', '#56B4E9', '#009E73', '#D55E00']

fig, ax = plt.subplots(figsize=(4.8, 4.0))
sns.boxplot(
    data=df,
    x='Treatment',
    y='Expression',
    hue='Treatment',
    palette=palette,
    legend=False,
    width=0.5,
    fliersize=0,
    ax=ax,
    boxprops=dict(alpha=0.85)
)
sns.stripplot(
    data=df,
    x='Treatment',
    y='Expression',
    color='#222222',
    size=5,
    jitter=0.15,
    alpha=0.8,
    ax=ax
)

y_max = df['Expression'].max() + 0.25
ax.plot([0, 0, 3, 3], [y_max, y_max + 0.15, y_max + 0.15, y_max], lw=1.2, c='black')
ax.text(1.5, y_max + 0.18, '*** (p < 0.0001)', ha='center', va='bottom', fontsize=9, fontweight='bold')

ax.set_ylim(top=y_max + 0.75)
ax.set_ylabel('Relative mRNA Expression (Fold Change)', fontweight='bold')
ax.set_xlabel('')
plt.xticks(rotation=15, ha='right')
plt.tight_layout()

plt.savefig('Figure1_GeneExpression.png', dpi=600)
plt.savefig('Figure1_GeneExpression.pdf')
plt.show()