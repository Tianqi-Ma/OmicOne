# OmicOne

> 🌐 **语言：** **中文** · [English README](README.md)

**在你自己的电脑上做交互式多组学分析。** 一条命令拉起浏览器界面，选择你手上的数据类型，
然后一步步走完那条流水线——**全部用你自己电脑的 CPU 和内存计算，无需云服务器，数据不上传任何地方**。

每一步都同时照顾新手和老手：

- 💡 通俗的**「这一步是什么？」**解释卡（配一个例子）
- 🔧 **方法可选**（每步都有备选，如降维 UMAP 或 t-SNE）
- 🎚️ **可调阈值**，并给出合理默认值
- 📊 **结果总结** + 大幅**预览图**
- 🧾 **可复现日志**，可导出成 R 脚本或叙述式报告
- 🎞️ 每一步的空白预览区都有一段**示意动画**，演示这一步在做什么
  （全部动画可在 `docs/explain-gallery.html` 预览）
- 🌗 浅色 / 深色主题（跟随系统），English / 中文 —— 顶栏一键切换

> **已改名：** 本项目在 v0.5.0 之前叫 **OMICstudio**。该名字与
> [OmicStudio](https://www.omicstudio.cn)（一个已有的在线分析平台）重名——对一个
> 主打「数据不离开你的机器」的工具来说，这个撞名尤其容易误导。代码不变，只是换名，
> 详见 [NEWS.md](NEWS.md)。

> **当前状态（v0.6.0）：** **单细胞 RNA-seq**（21 步）与 **WES / 体细胞突变**（12 步）
> 两条流程均已完整并通过静态校验。其余三个组学在界面中展示规划路线图，尚未实现。
> 两条流程都**还没有在装好各自引擎的环境里端到端实跑过**，见[注意事项](#注意事项)。

---

## 五条流水线

在起始页点击卡片即可进入对应流程。

| 组学 | 引擎 | 状态 |
|---|---|---|
| **单细胞 RNA-seq** | [scop](https://github.com/mengxu98/scop) + Seurat | ✅ 已完成（21 步） |
| **WES / 体细胞突变** | [maftools](https://bioconductor.org/packages/maftools) | ✅ 已完成（12 步） |
| **Bulk RNA-seq** | TOmicsVis + DESeq2 | 🚧 界面内展示路线图 |
| **空间转录组** | Seurat + SpatialExperiment | 🚧 界面内展示路线图 |
| **多组学整合** | MOFA2 / iClusterPlus / SNFtool | 🚧 界面内展示路线图 |

临床随访数据是**跨流程共享**的：在单细胞的「临床与生存」步骤加载一次队列，WES 的
「突变-预后」步骤会自动接手——同样的曲线、同样的 log-rank 检验、同一套代码。

maftools 是普通的 Bioconductor 二进制包，因此 WES 流程无需源码编译，在禁止编译的机器上也能跑。

---

## 为什么"localhost 优先"？

真实的组学分析（Seurat/Bioconductor）需要原生计算和真实内存，纯浏览器（WASM）跑不动。
所以 OmicOne 采用**本地服务 + 浏览器界面**：界面是网页，但所有计算都在**你的机器**上完成。
这和 `cellxgene launch` 是同一个模式。

---

## 三种运行方式

按你愿意安装多少东西来选。

### A. 你已经装了 R（最轻）
```r
# install.packages("remotes")
remotes::install_github("Tianqi-Ma/OmicOne")
OmicOne::run_app()   # 自动打开浏览器
```
需要 R ≥ 4.1，**shiny ≥ 1.7.4** 与 **bslib ≥ 0.7.0**（推荐 shiny ≥ 1.8.1），
以及你实际要跑的那些步骤所需的重分析包。下载最小，一条命令。

### B. 没有 R、零依赖 → Docker（多数用户推荐）✅
R + Seurat + Bioconductor + scop + 应用**全部打进一个镜像**，你只需装
[Docker](https://www.docker.com/products/docker-desktop/)。
```bash
docker build -t omicone .           # 从本仓库构建一次
docker run --rm -p 3838:3838 -m 16g omicone
# 然后浏览器打开 http://localhost:3838
```
数据通过**浏览器上传**（无需挂载目录）。给 Docker 足够内存（`-m 16g`，大数据还需在
Docker Desktop 里调高内存上限）。

### C. 纯小白双击即用（规划中）
用桌面安装包（Tauri/Electron/electricShine）把 R 和应用打包，用户**双击即可**——无需命令行、
无需 Docker。计划在后续版本推出。

---

## 零门槛试用（无需自己的数据）

在**导入**步骤，保持数据来源为 **演示数据**，点击"加载演示数据"即可。提供三个演示数据集：

| 演示数据 | 内容 | 需要 |
|---|---|---|
| **PBMC 3k**（默认） | 经典的 2,700 细胞 10x PBMC 数据集，**已内置在包里** | 无——即时、离线 |
| **Pancreas** | 小鼠胰腺，含谱系结构与 spliced/unspliced 层，适合轨迹 / 速率 | 需安装 `scop` |
| **微型示例** | 500 × 300 合成矩阵，用于快速检查界面 | 无 |

你也可以粘贴任意 `.rds` / `.h5` 文件的 URL 在线加载。

## 可上传的数据

- **RDS** — 保存的 `Seurat` 或 `SingleCellExperiment` 对象（旧版 Seurat 对象会自动升级）
- **10x `.h5`** — Cell Ranger 的 HDF5
- **count 表格** — CSV/TSV，基因为行、细胞为列
- **临床表格** — CSV/TSV，每位患者一行（ID、随访时间、终点事件）

浏览器上传上限调得很高（默认 5 GB；用 `run_app(max_upload_mb = ...)` 修改）。

---

## 单细胞分析步骤

7 个阶段，21 个步骤。**加粗**为默认方法。

| 阶段 | # | 步骤 | 可选方法 |
|---|---|------|-----------------|
| **数据与质控** | 1 | 导入 | RDS / 10x .h5 / 表格 / URL / 演示数据 |
| | 2 | 质控 QC | **MAD 自适应** / 手动；线粒体、核糖体、血红蛋白、解离基因占比 |
| | 3 | 去双细胞 | **scDblFinder** / DoubletFinder |
| **预处理** | 4 | 归一化 | **LogNormalize** / SCTransform |
| | 5 | 特征选择 / PCA | HVG **vst** / mvp / dispersion |
| | 6 | 批次整合 | **无** / Harmony / CCA / RPCA（Seurat v5），附批次混合度评分 |
| **结构** | 7 | 聚类 | **Leiden** / Louvain |
| | 8 | 降维图 | **UMAP** / t-SNE |
| **身份** | 9 | 标志基因 | **wilcox** / roc / MAST |
| | 10 | 注释 | **手动** / SingleR（按簇，pruned 标签）/ Azimuth |
| | 11 | 富集 / GSEA | ORA + GSEA（经 scop 调用 clusterProfiler） |
| **轨迹与动态** | 12 | 轨迹 | **Slingshot** / Monocle2 / Monocle3 / PAGA / Palantir |
| | 13 | RNA 速率 | scVelo（steady-state / stochastic / dynamical） |
| | 14 | 动态特征 | scop 动态特征 + 热图 |
| **高级** | 15 | 周期与信号 | Seurat 细胞周期打分、UCell / AddModuleScore |
| | 16 | 细胞通讯 | LIANA / CellChat |
| | 17 | 恶性 / CNV | CopyKAT、干性打分 |
| | 18 | **临床与生存** | Kaplan-Meier + log-rank + 单因素 Cox；可按临床列分组，也可按每样本细胞类型占比分组（中位 / 三分位 / 最优切点） |
| **产出** | 19 | 可视化 | UMAP / violin / dotplot / feature / heatmap |
| | 20 | 报告 | 由可复现日志生成叙述式 HTML 报告 |
| | 21 | 导出 | .rds / .h5ad（尽力） / 图 / R 脚本 |

顶栏的 **⤓ 导出**菜单在**每一步**都可用（对象、元数据、表达矩阵、降维坐标），
不必走到最后一步才能把数据取出来。

方法与默认值遵循当前（2023–2025）最佳实践：MAD 质控、scDblFinder、Leiden、Harmony，
全流程使用 scop/SCP 绘图风格。

---

## WES 分析步骤

4 个阶段，12 个步骤，全部基于 maftools。

| 阶段 | # | 步骤 | 内容 |
|---|---|------|-----------------|
| **输入** | 1 | 导入 MAF | MAF（可附临床表），或内置的 TCGA LAML 演示队列 |
| | 2 | 队列概览 | 变异分类、每样本突变负荷、高频突变基因 |
| **突变全景** | 3 | Oncoplot | 瀑布图，可叠加临床注释条 |
| | 4 | TiTv / VAF / rainfall | 突变谱、等位基因频率、kataegis |
| | 5 | TMB | 按你的捕获区域大小计算每 Mb 突变数 |
| | 6 | Lollipop / 结构域 | 单个基因的突变沿蛋白分布，叠加结构域 |
| | 7 | 驱动基因与互作 | oncodrive 位点聚集；共现 / 互斥 |
| **突变特征** | 8 | 突变特征 | 三核苷酸矩阵 → de-novo 特征 → COSMIC 匹配 |
| **临床与预后** | 9 | 临床/通路/药物 | 按临床分组的基因富集、致癌通路、药物-基因相互作用 |
| | 10 | 队列比较 | 两个临床分组间的 Fisher 检验 + 森林图 |
| | 11 | 突变-预后 | 突变型 vs 野生型的 KM + log-rank，走**共享**生存分析层 |
| | 12 | 异质性 | 逐样本 VAF 聚类、MATH 分数 |

第 1 步提供 maftools 自带的 **TCGA LAML** 队列（193 个样本），因此整条流程无需自备数据即可离线体验。

---

## 注意事项

- **哪些已经实跑过、哪些没有（2026-10 审计）。** 单细胞核心链路（导入 → QC → 去双细胞 →
  归一化 → PCA → 整合 → 聚类 → UMAP → 标志基因）和临床/生存统计已在自带数据上用 Seurat 4.3
  实跑，复现日志里的 R 代码重放结果完全一致；WES 流程在 maftools 自带的 TCGA-LAML 示例上实跑过
  （maftools 2.14）。依赖 scop 的步骤（富集、轨迹、速率、动态特征）以及 CellChat / LIANA、copykat、
  SingleR、Azimuth、scDblFinder 已对照各自最新源码核对过参数，但**没有实际运行**；正式使用前请先在
  你自己的环境里跑一遍。
- **`scop` 是单细胞流程的绘图与计算引擎**，它是一个 GitHub 包；如果你的机器禁止源码编译，
  请使用 Docker 镜像——镜像已内置 scop，并预先烤好了 Python 步骤（scVelo、PAGA、Palantir）
  所需的 conda 环境。
- **重依赖在安装时是可选的**（放在 `Suggests`）。缺依赖时 UI 也能起；每个计算步骤会检查所需包，
  缺了会提示你安装。
- **部分步骤首次需联网**（SingleR/Azimuth 下载参考集）。
- **`.h5ad` 导出**在纯 R 下依赖 SeuratDisk（尽力而为）；主导出格式是 `.rds`。
- **生存分析的「最优切点」仅供探索。** 它是为了让分离最大化而挑出来的，因此 p 值偏乐观——
  正式结论请用中位切点，或在独立队列中验证该切点。
- **突变特征需要 BSgenome 包**（`BSgenome.Hsapiens.UCSC.hg19` 或 hg38）以及 `NMF`。
  这是个较大的下载；WES 的其他步骤都不需要它。
- **依赖 VAF 的步骤需要 VAF 列。** 很多变异检出工具不报告该值，此时 VAF、rainfall
  与异质性页签会明确提示，而不是猜一个。
- **内存**随细胞数增长。<10 万细胞在 16–32 GB 上较舒适；更大需更多内存。应用会警告并提供下采样。

---

## 开发

R 包结构（golem 风格）：

```
R/app_ui.R            顶层外壳（顶栏、侧边栏、开屏动画）
R/app_landing.R       组学选择落地页
R/app_server.R        组学路由 + 共享枢纽（rv$obj、rv$clinical、rv$status）
R/steps.R             各组学的步骤注册表——唯一数据源
R/mod_*.R             每步一个模块
R/fct_compute.R       scop/Seurat 计算封装
R/fct_wes.R           maftools 封装 + MAF→生存分析的桥接
R/fct_plots.R         scop 绘图封装 + 共享主题/配色
R/fct_survival.R      与组学无关的生存分析（survival + ggplot2）
R/utils_ui.R          共享的步骤布局与双语 helper
```

开发时本地运行：
```r
pkgload::load_all(); run_app()
```

提交前必须通过的校验：
```r
testthat::test_dir("tests/testthat")   # 其中包含两道防线：每个模块的输出都能求值，
                                       # 以及没有任何函数调用了不存在的名字
```
```sh
tools/check_mirror.sh                  # 共享文件除包名外必须与 scStudio 完全一致
```

单细胞流程同时以 [scStudio](https://github.com/Tianqi-Ma/scStudio) 独立维护；
共通改动会在两个仓库之间同步。

## 许可证

MIT。
