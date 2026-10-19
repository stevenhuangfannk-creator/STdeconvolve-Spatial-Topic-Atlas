# 数据来源 / Data provenance

核心数据是官方 STdeconvolve 包内 mOB。安装固定提交后，`data("mOB", package="STdeconvolve")` 可加载，无须提供个人数据。下载与保存代码见 `01_Official_Reproduction.R` 和 `environment/install_dependencies.R`。

- 官方源码：[固定提交](https://github.com/JEFworks-Lab/STdeconvolve/tree/13b1953263a5ce34e2471c4defebedbda44c7d08)。
- 原始平台论文：Ståhl et al., Science 2016，[DOI](https://doi.org/10.1126/science.aaf2403)。
- mOB 原始整数矩阵：15,928 genes × 262 spots；坐标 x/y 是原始阵列坐标，不在本项目中宣称微米单位。
- `pos` 与 `counts` 顺序不同，始终以 spot ID 对齐。`annot` 仅覆盖清理后的 260 spots，不能要求原始 262 spots 全有标签，也不能按位置强行配对。
- `annot` 是五个组织层，不是细胞类型比例真值。不能从该数据计算对真实细胞比例的误差。
- 第二切片与功能数据库若下载，来源、下载时刻及校验值由扩展脚本记录；不要把同一切片随机拆分当作独立样本。

`raw/`、二进制模型和本机包库保留在本地，由 `.gitignore` 排除；GitHub 保存小型结果表、参数、图件和下载方法。实际文件校验见 `data_manifest.csv`（生成后可查）。上游数据遵循原始提供方条款，本仓库 GPL-3 许可不重新授权第三方数据。
