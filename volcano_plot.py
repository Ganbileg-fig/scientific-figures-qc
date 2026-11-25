import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 10,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'figure.dpi': 300
})

np.random.seed(101)
n_genes = 2500
genes = [f"GENE_{i:04d}" for i in range(1, n_genes + 1)]

top_up = {'COX2': (3.4, 1e-14), 'IL6': (2.9, 1e-11), 'TNF': (2.5, 1e-9), 'MMP9': (2.1, 1e-7)}
top_down = {'PPARG': (-3.2, 1e-13), 'SIRT1': (-2.6, 1e-10), 'CAT': (-2.2, 1e-8), 'SOD2': (-1.8, 1e-6)}

log2fc = np.random.normal(0, 0.9, n_genes)
padj = np.random.uniform(1e-4, 1.0, n_genes)

df = pd.DataFrame({
    'Gene': genes,
    'log2FC': log2fc,
    'padj': padj
})

for g, (fc, p) in {**top_up, **top_down}.items():
    df.loc[len(df)] = [g, fc, p]

df['neg_log10_padj'] = -np.log10(df['padj'])

fc_thresh = 1.0
p_thresh = -np.log10(0.05)

conditions = [
    (df['log2FC'] > fc_thresh) & (df['neg_log10_padj'] > p_thresh),
    (df['log2FC'] < -fc_thresh) & (df['neg_log10_padj'] > p_thresh)
]
choices = ['Upregulated', 'Downregulated']
df['Status'] = np.select(conditions, choices, default='Not Significant')

fig, ax = plt.subplots(figsize=(5.5, 4.5))

ns = df[df['Status'] == 'Not Significant']
up = df[df['Status'] == 'Upregulated']
dn = df[df['Status'] == 'Downregulated']

ax.scatter(ns['log2FC'], ns['neg_log10_padj'], c='#CCCCCC', s=10, alpha=0.5, edgecolors='none', label='NS')
ax.scatter(dn['log2FC'], dn['neg_log10_padj'], c='#0072B2', s=14, alpha=0.8, edgecolors='none', label='Downregulated')
ax.scatter(up['log2FC'], up['neg_log10_padj'], c='#D55E00', s=14, alpha=0.8, edgecolors='none', label='Upregulated')

ax.axvline(fc_thresh, color='#666666', linestyle='--', lw=0.9)
ax.axvline(-fc_thresh, color='#666666', linestyle='--', lw=0.9)
ax.axhline(p_thresh, color='#666666', linestyle='--', lw=0.9)

labels_to_show = list(top_up.keys()) + list(top_down.keys())
sub_labels = df[df['Gene'].isin(labels_to_show)]

for _, row in sub_labels.iterrows():
    x_offset = 12 if row['log2FC'] > 0 else -14
    y_offset = 8
    ax.annotate(
        row['Gene'],
        xy=(row['log2FC'], row['neg_log10_padj']),
        xytext=(x_offset, y_offset),
        textcoords='offset points',
        fontsize=8,
        fontweight='bold',
        fontstyle='italic',
        color='#111111',
        arrowprops=dict(arrowstyle='-', color='#444444', lw=0.7, shrinkA=0, shrinkB=2)
    )

ax.set_xlabel('log₂ (Fold Change)', fontweight='bold')
ax.set_ylabel('-log₁₀ (Adjusted p-value)', fontweight='bold')
ax.set_title('Differential Expression Profile', fontweight='bold', pad=10)
ax.legend(frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.14), ncol=3, fontsize=8.5)

plt.tight_layout()
plt.savefig('Figure2_VolcanoPlot.png', dpi=600, bbox_inches='tight')
plt.savefig('Figure2_VolcanoPlot.pdf', bbox_inches='tight')
plt.show()