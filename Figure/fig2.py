# 载入库
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import kruskal

# 载入数据  "D:/OneDrive/OneDrive - International Campus, Zhejiang University/大三下/PoN/ICA/Data/number_ITI_poke.csv"
number_ITI_poke = pd.read_csv("D:/OneDrive/OneDrive - International Campus, Zhejiang University/大三下/PoN/ICA/Data/number_ITI_poke.csv")
experimental_conditions = pd.read_csv("D:/OneDrive/OneDrive - International Campus, Zhejiang University/大三下/PoN/ICA/Data/Experimental_conditions.csv")

# 整理Recall阶段数据
recall_days_iti = number_ITI_poke[["Mouse no.", "Recall day 1", "Recall day 2", "Recall day 3", "Recall day 4", "Recall day 5"]]
iti_recall_long = recall_days_iti.melt(id_vars="Mouse no.", var_name="Recall Day", value_name="iti_poke")
iti_recall_long["Recall Day"] = iti_recall_long["Recall Day"].str.extract(r"(\d+)").astype(int)
iti_recall_long = iti_recall_long.merge(experimental_conditions[["Mouse no.", "NA/NE"]], on="Mouse no.")

# 绘制Boxplot
plt.figure(figsize=(8, 6))
sns.boxplot(data=iti_recall_long, x="NA/NE", y="iti_poke", palette="Set2")
plt.title("Figure 2. Recall Phase ITI Poke Comparison across NA/NE Groups")
plt.xlabel("NA/NE Group")
plt.ylabel("ITI Poke Count")
plt.xticks(ticks=[0, 1, 2], labels=["Clonidine (-1)", "Saline (0)", "Yohimbine (1)"])
plt.grid(True)
plt.tight_layout()
plt.show()

# Kruskal-Wallis总体检验
h_stat, p_value = kruskal(
    iti_recall_long.loc[iti_recall_long["NA/NE"] == -1, "iti_poke"],
    iti_recall_long.loc[iti_recall_long["NA/NE"] == 0, "iti_poke"],
    iti_recall_long.loc[iti_recall_long["NA/NE"] == 1, "iti_poke"]
)
print(f"Kruskal-Wallis H = {h_stat:.2f}, p = {p_value:.4f}")
# Post-hoc Dunn's test（要求scikit_posthocs库）
# 安装: pip install scikit-posthocs

import scikit_posthocs as sp

# 使用Bonferroni校正的Dunn's Test
dunn_results = sp.posthoc_dunn(
    iti_recall_long, 
    val_col="iti_poke", 
    group_col="NA/NE", 
    p_adjust="bonferroni"
)

print(dunn_results)
