import pandas as pd
import scanpy as sc

# 设定路径
input_txt = "/root/Desktop/my_pan/workspace/Data/raw_data/Marrow_raw/BoneMarrow1_dge.txt"
output_h5ad = "/root/Desktop/my_pan/workspace/Data/pipeline_output/BoneMarrow1.h5ad"

# 读取 txt，自动识别分隔符（空格或制表符）
df = pd.read_csv(input_txt, sep=None, engine="python", index_col=0)

adata = sc.AnnData(df.T)

# 添加 obs 和 var 信息
adata.obs["batch"] = "BoneMarrow1"
adata.var_names = df.index.astype(str)  # 基因名
adata.obs_names = df.columns.astype(str)  # 细胞名
adata.var = pd.DataFrame(index=adata.var_names)

# 保存为 .h5ad
adata.write(output_h5ad)
print(f"✅ 转换完成，保存为：{output_h5ad}")

