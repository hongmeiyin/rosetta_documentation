#How to prepare structures for use in Rosetta
>本页描述了为 Rosetta 准备蛋白质结构的标准流程。如果你需要准备配体（ligand）用于 Rosetta，请参阅 [[preparing ligands]] 页面以及 教程链接。

*This page describes the standard procedure for preparing protein structures for Rosetta. To prepare ligands for use with Rosetta, see the [[preparing ligands]] page and [tutorial](https://www.rosettacommons.org/demos/latest/tutorials/prepare_ligand/prepare_ligand_tutorial).
>如需准备非肽聚合物结构，请查看 [[preparing PDB files for non-peptide polymers]] 页面。

To prepare non-peptide polymers, see the [[preparing PDB files for non-peptide polymers]] page.
>对于包含蛋白质和 RNA 的 PDB 文件，请参考 [[RNA-protein changes]] 页面。

For PDB files containing both proteins and RNA, see the [[RNA-protein changes]] page.
You can find the introductory tutorial [here](https://www.rosettacommons.org/demos/latest/tutorials/input_and_output/input_and_output#controlling-input_preparing-a-structure-by-refinement). For details, continue reading.

[[_TOC_]]

# Application Purpose
>从 RCSB PDB 下载的结构通常并不直接适用于 Rosetta，常见的问题包括：
>原子之间存在碰撞（clash）
>氨基酸的构象能量过高（如 rotamer 处于不合理状态）
>其他异常错误

Structures derived straight from the PDB are almost never perfectly compatible with Rosetta—it is common for them to have clashes (atom overlaps), amino acid rotamers with terrible Rosetta energy, or other errors. 
It is often beneficial to prepare the structures before doing real work on them to get these errors out of the way beforehand. This provides several benefits:
>结构设计过程中，Rosetta 可能会因为输入结构中原有残基冲突或构象问题而错误替换氨基酸。预处理可保留野生型氨基酸，避免不必要更改。

-   In design, Rosetta often places a certain amino acid solely because the original has a bad clash or rotamer energy, preparing structures first will allow Rosetta to keep the wild-type amino acid and minimize the number of changes made based on incorrect scoring of the native.
>各路径模拟时无需重复修复同一结构问题，提升效率。

-   Less time is spent in each trajectory independently re-relaxing the same errors in the input.
>输出结果更一致，不因路径差异造成误差。

-   The results have less noise caused by errors being handled in different ways in different trajectories.
>最终结构评分更低（更合理）——一个良好折叠的蛋白不应在 score12 或 talaris 能量函数中出现正值。

-   Scores are lower overall—you should never have positive score12 (or talaris) scores for a well folded protein.

# References

Nivón LG, Moretti R, Baker D (2013) A Pareto-Optimal Refinement Method for Protein Design Scaffolds. PLoS ONE 8(4): e59004. [Paper](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0059004)

# How Do I Prepare Structures?
>结构准备方式与你想要进行的模拟类型紧密相关。

Preparing structures is inextricably linked to what you want to do with them. 
>本质上，你的主模拟流程决定了结构准备流程。

In other words, your main protocol dictates your preparation protocol. 
>Rosetta 的 relax 操作并不总是让结构“更接近真实”，而是让结构更符合 Rosetta 的能量函数。你是在“调整结构以让 Rosetta 更喜欢它”。

Remember that all you're really doing here is relaxing into Rosetta's energy function—you're not necessarily making it objectively more correct (although clashes are generally wrong), you're really just making Rosetta like it better.
>以下是一个针对酶设计的 benchmark 流程，也适用于大多数设计场景。其他用途可按需调整。

What follows is a protocol that has been benchmarked for enzyme design and should work for any design situation, and then a series of suggestions from Rosetta developers for other situations.

# Cleaning PDBs for Rosetta
>从 PDB 官网下载的文件常常包含 Rosetta 无法处理的信息，例如水分子（Rosetta 默认不处理），非标准残基等。

Parsing a PDB structure downloaded from the [[PDB|http://www.rcsb.org/pdb/home/home.do]] often results in an error.
This is mostly due to the PDB containing extraneous information (such as waters which are not modeled by Rosetta).
Occasionally, this is due to Rosetta not being able to parse a particular chemical entity.

The simplest way to clean a PDB file is by running the following command-line expression:

```
grep "^ATOM" SOME.pdb > SOME.clean.pdb
```
>这种方式较为粗暴，会删除 PDB 文件中可能有用的信息。

**Warning:** this is rather barbaric and will get rid of PDB information which could be valuable.
>该脚本仅保留 "ATOM" 和 "TER" 行，移除 "HETATM"（除非是 selenomethionine）、0 占据率的残基等。

There is also a script, `tools/protein_tools/scripts/clean_pdb.py`, for cleaning the PDB which will remove most lines that cannot be parsed by Rosetta.
As of June 9th, 2015, the script leaves only lines starting with 'ATOM' or 'TER' and removes 'HETATM' records unless the residue is a selenomethionine.
It also removes residues with 0 occupancy.
See the script header for more information.
To run the script:

```
python clean_pdb.py <pdb> <chain id>
```

## Clean PDBs with Ligand
>如果你的 PDB 文件包含 ligand（HETATM），并希望保留它，可以使用以下脚本：

There are cases where your PDB may contain a HETATM ligand, for obvious reasons you may want to relax the PDB with the ligand.
We often clean structures to replace non-canonical amino acids with their closest counterparts using `clean_pdb_keep_ligand.py`:
```
python clean_pdb_keep_ligand.py your_structure_original.pdb -ignorechain
```
>此脚本可用于保留链和配体，适合单链 PDB。多链结构可分开处理，或使用 -ignorechain 保留所有链。

These protocols are designed for a single-chain PDB. 
For multiple chains we recommend that you split the PDB into one for each chain and run the protocol separately on each. Alternatively, you can try using the `-ignorechain` option of the script, which will keep all chains.

While the previous script (clean_pdb.py) discarded most to all HETATM records, this script keeps HETATMs which are not non-canonical ligands(非规范配体).

# Relax With All-Heavy-Atom Constraints: Introduction

(See also the [[relax documentation|relax]] .)
>我们寻找了一种方法，既能同时最小化Rosetta能量，又能使晶体结构中的所有重原子尽可能保持起始位置。
>正如下文许多帖子——或来之不易的经验——所表明的那样，对结构运行relax操作通常会使主链移动几埃。
>我们目前发现的最佳同步优化方法是：始终开启约束运行relax（通常在relax运行后期循环中约束会逐步减弱），并且不仅约束主链原子，还要约束侧链原子（可通过flags或[[约束文件|constraint-file]]实现，参见[[preparing-structures#generating-constraints-file-for-your-pdb]]）。
>该方案已在51个蛋白质的基准测试集中验证，与原始PDB结构中的设计相比，酶设计的序列恢复率提高了5%。
>在从原始PDB到"带约束relax处理PDB"的过程中，整个蛋白质组的C-alpha原子均方根偏差（RMSD）仅为0.077埃。

We looked for a way to simultaneously minimize Rosetta energy and keep all heavy atoms in a crystal structure as close as possible to their starting positions. 
As many posts below—or hard-won experience—will show, running relax on a structure will often move the backbone a few Angstroms. 
The best way we have found to perform the simultaneous optimization is to run relax with constraints always turned on (typically constraints ramp down in the late cycles of a relax run) and to constrain not just backbone but also side-chain atoms (which can be done with both flags or a [[constraint file|constraint-file]], see [[preparing-structures#generating-constraints-file-for-your-pdb]]).
This protocol has been tested on a benchmark set of 51 proteins and found to increase sequence recovery in enzyme design by 5% as compared with design in raw PDB structures. 
It accomplishes this with only .077 Angstrom RMSD over the set of proteins (C-alpha RMSD) from raw PDB to relaxed-with-csts PDB. 
A more complete description of the data leading to this protocol is below.

# Relax With All-Heavy-Atom Constraints: Protocol
>所需的文件位于：`Rosetta/main/source/src/apps/public/relex_w_allatom_cst`。

The required files are in: `Rosetta/main/source/src/apps/public/relax_w_allatom_cst`. 
>还有一个带有约束文件而不是flags的演示：`Rosetta/dedemos/public/presee_pdb_for_Rosetta_with_relax`。

There is also a demo with a constraint file rather than flags at: `Rosetta/demos/public/prepare_pdb_for_rosetta_with_relax`.

### Short Protocol (recommended, no constraint file)短协议（推荐，无约束文件）
>Relax应用程序本身内置了所有重原子约束的Relax。

Relax with all-heavy-atom constraints is built into the relax application itself. 
>如果这是一个新结构，您可能需要首先使用上述脚本对其进行清理。Relax的流程如下：

If this is a new structure you may want to first clean it up using the above script. Relax proceeds as follows:

```
relax.linuxgccrelease  -database Rosetta/main/database [-extra_res_fa your_ligand.params] -relax:constrain_relax_to_start_coords -relax:coord_constrain_sidechains -relax:ramp_constraints false -s your_structure.pdb
```
>标志文件可以包含您希望使用的任何打包和记分功能标志。我们建议至少：
The flags file can contain whatever packing and scorefunction flags you wish to use. We recommend at least:

```
-ex1
-ex2
-use_input_sc
-flip_HNQ
-no_optH false
```
>使用的约束的强度和类型可以根据"-relax:coord_cst_stdev"和"-relax:coord_cst_width"选项而变化。

The strength and type of constraints uses can be varied with the `-relax:coord_cst_stdev` and `-relax:coord_cst_width` options. 
>默认情况下，relax使用调和约束，其强度由"coord_cst_stdv"调整（更小=更紧）。

By default relax uses a harmonic constraint with the strength adjusted by `coord_cst_stdev` (smaller=tighter). 
>如果指定了`coord_cst_width`，则使用平底线性壁约束，平底井的大小由`coord_cst_widt`控制（更小=更紧），壁的坡度由`coold_cst_stdv`控制（更小=更紧密）。

If `coord_cst_width` is specified, a flat-bottomed, linear-walled constraint is used, with the size of the flat-bottomed well controlled by `coord_cst_width` (smaller=tighter), and the slope of the walls by `coord_cst_stdev` (smaller=tighter). 
More on constraint functions can be found [[here|constraint-file#function-types]].

### Longer Protocol (not recommended, using a constraint file)更长的协议（不建议使用约束文件）
>一般来说，对于大多数应用程序来说，短协议是首选，因为这个版本更复杂，两者给出的结果几乎相同。

In general the short protocol is preferred for most applications, since this version is more complicated and the two give nearly identical results. 
>在此协议中，一个单独的脚本首先从输入PDB生成侧链原子约束，然后使用此预生成的约束文件运行relax协议。

In this protocol an separate script first generates side-chain atom constraints from an input PDB, then the relax protocol is run with this pre-generated constrain file. 
>较短的协议一步完成了这一切，**这个版本在很大程度上被弃用了**。

The shorter protocol does this all in one step, and **this version is largely deprecated**. 
>某些用户可能更喜欢此协议，因为它允许您在放松之前查看所有约束的列表，并可能使用其他脚本/数据修改约束。

Certain users might prefer this protocol because it allows you to see a list of all constraints, and perhaps to modify constraints using other scripts/data, prior to relax.
>使用"sidechain_cst_3.py"在PDB上生成侧链坐标约束：

1. Generate side-chain coordinate constraints on your PDB using `sidechain_cst_3.py`:

    ```
    python sidechain_cst_3.py your_structure.pdb 0.1 0.5
    ``` 
    [output: your\_structure\_sc.cst]
>使用这些约束和自定义的relax脚本运行relax，以强制约束在整个运行过程中保持不变。

2. Run relax using these constraints and with a custom relax script to force constraints to stay on during the entire run.
>对于340个残基的蛋白质，运行时间约为30分钟。

Run time is approximately 30 minutes for a 340 residue protein.

    ```
    relax.linuxgccrelease  -database Rosetta/main/database -relax:script Rosetta/main/source/src/apps/public/relax_w_allatom_cst/always_constrained_relax_script -constrain_relax_to_start_coords -constraints:cst_fa_file your_structure_sc.cst -s your_structure.pdb [-extra_res_fa your_ligand.params]
    ```
>括号内的"extra_res_fa"仅适用于结构中有任何配体的情况。

The bracketed `extra_res_fa` only applies if there are any ligand(s) in the structure. 
The example flags file listed below was used in testing:

>pro_hydroxylated_case1    脯氨酸羟基化（情况1）
>pro_hydroxylated_case2    脯氨酸羟基化（情况2）
>ser_phosphorylated        丝氨酸磷酸化
>thr_phosphorylated        苏氨酸磷酸化
>tyr_phosphorylated        酪氨酸磷酸化
>tyr_sulfated              酪氨酸硫酸化
>lys_dimethylated          赖氨酸二甲基化
>lys_monomethylated        赖氨酸单甲基化
>lys_trimethylated         赖氨酸三甲基化
>lys_acetylated            赖氨酸乙酰化
>glu_carboxylated          谷氨酸羧基化
>cys_acetylated            半胱氨酸乙酰化
>tyr_diiodinated           酪氨酸二碘化
>N_acetylated              N端乙酰化
>C_methylamidated          C端甲基酰胺化
>MethylatedProteinCterm    蛋白质C端甲基化
>这些术语描述了蛋白质翻译后修饰（Post-Translational Modifications, PTMs），即蛋白质在合成后发生的化学修饰，可能影响其功能、定位或稳定性。具体分类如下：
> 1. 氨基酸残基修饰

> 羟基化（Hydroxylation）pro_hydroxylated：脯氨酸（Proline）的羟基化，常见于胶原蛋白，影响结构稳定性。case1/case2 可能指不同位点的修饰（如 3-羟基脯氨酸 vs 4-羟基脯氨酸）。

> 磷酸化（Phosphorylation）ser_phosphorylated、thr_phosphorylated、tyr_phosphorylated：分别在丝氨酸（Serine）、苏氨酸（Threonine）、酪氨酸（Tyrosine）上添加磷酸基团，调控信号传导（如激酶/磷酸酶作用）。

> 硫酸化（Sulfation）tyr_sulfated：酪氨酸的硫酸化，常见于分泌蛋白和细胞表面受体（如趋化因子）。

> 甲基化（Methylation）lys_mono/di/trimethylated：赖氨酸（Lysine）的甲基化（单/二/三），影响组蛋白功能（表观遗传调控）。

> C_methylamidated：C端甲基酰胺化，可能保护蛋白质免受降解。

> 乙酰化（Acetylation）lys_acetylated、cys_acetylated：赖氨酸或半胱氨酸（Cysteine）的乙酰化，调控蛋白质互作或活性（如组蛋白乙酰化激活转录）。

> 羧基化（Carboxylation）glu_carboxylated：谷氨酸（Glutamate）的羧基化，见于凝血因子（如凝血因子Ⅱ、Ⅶ、Ⅸ、Ⅹ依赖维生素K的修饰）。

> 碘化（Iodination）tyr_diiodinated：酪氨酸的二碘化，参与甲状腺激素（T3/T4）合成。

> 3. 蛋白质末端修饰

> N端修饰 N_acetylated：N端乙酰化，常见于真核蛋白，影响稳定性和定位。

> C端修饰 MethylatedProteinCterm：蛋白质C端甲基化，可能调控蛋白质-蛋白质相互作用或半衰期。
```
-ex1
-ex2
-use_input_sc
-correct
-no_his_his_pairE
-score::hbond_params correct_params
-lj_hbond_hdis 1.75
-lj_hbond_OH_donor_dis 2.6
-linmem_ig 10
-nblist_autoupdate true
-no_optH false
-flip_HNQ
-chemical:exclude_patches LowerDNA UpperDNA Cterm_amidation SpecialRotamer
VirtualBB ShoveBB VirtualDNAPhosphate VirtualNTerm CTermConnect sc_orbitals
pro_hydroxylated_case1 pro_hydroxylated_case2 ser_phosphorylated
thr_phosphorylated tyr_phosphorylated tyr_sulfated lys_dimethylated
lys_monomethylated lys_trimethylated lys_acetylated glu_carboxylated
cys_acetylated tyr_diiodinated N_acetylated C_methylamidated
MethylatedProteinCterm
```

Note: Including extra rotamers is important if your goal is to keep all side-chain atoms tightly constrained; if that is not important for your applications, exclude the ex1 and ex2 for speed.
>注意：如果你的目标是保持所有侧链原子的紧密约束，那么包括额外的旋转异构体是很重要的；如果这对您的应用程序不重要，请排除ex1和ex2以提高速度。

# Relax With All-Heavy-Atom Constraints: Data 放松所有重原子约束：数据
> 为了测试松弛输入PDB的协议，我们使用了51个用于酶的支架测试集

To test protocols for relaxation of input PDBs we used the 51 scaffold test set for enzdes. 


We run design over the input structures and calculate the percent of residues which come back with the native identity—sequence recovery, only over the designed residues. 
>我们对输入结构进行设计，并计算出仅在设计的残基上返回的具有天然身份的残基百分比——序列恢复。

This is of course an imperfect metric (the original sequence might not be fully optimal), but it allows us to ask how many residues Rosetta will correctly choose, assuming that the input structure is already at a minimum in sequence space for the ligand in question.
>这当然是一个不完美的度量（原始序列可能不是完全最优的），但它允许我们询问Rosetta将正确选择多少个残基，假设输入结构在所讨论的配体的序列空间中已经处于最小值。

It had already been found that relax alone (with no constraints) will distort most structures, and that those structures will give a much higher sequence recovery in design, but this is a result of the distortion to the input structure.
>已经发现，单独放松（没有约束）会扭曲大多数结构，这些结构在设计中会给出更高的序列恢复，但这是输入结构扭曲的结果。

We decided to test a few protocols to find the one that would best minimize RMSD from the original PDB while maximizing sequence recovery. 
>我们决定测试一些协议，以找到一种在最大限度地减少原始PDB的RMSD的同时最大限度地提高序列恢复的协议。

All calculations are averages over 50 runs of the 51 scaffold set. 
>所有计算均为51个scaffold组50次运行的平均值。

All RMSDs are to native C-alpha. 
>所有RMSD都是native C-alpha。

It is also possible to read in electron density for a PDB and use that as a constraint during relax (described [[here|density-map-scoring]]), but electron density is not uniformly available, and we were looking for a protocol that would work in every case even if the electron density was lacking. 
>也可以读取 PDB 文件的电子密度图（electron density），并在结构松弛（relax）过程中将其作为约束（具体方法参见[[此处|density-map-scoring]]）。但电子密度数据并非普遍可用，因此我们希望开发一种即使在没有电子密度数据时也能适用的方案。

We chose the first protocol as the best because it minimized RMSD (and side-chain motions) while maximizing sequence recovery.  See [[Nivon et al|preparing-structures#References]] for more information.
>我们选择第一种方案作为最佳方案，因为它最大限度地减少了RMSD（和侧链运动），同时最大限度地提高了序列恢复率。

Running relax with side-chain coordinate constraints and backbone coordinate constraints: (ex flags and use native; ex flags on in enzdes)
>使用侧链坐标约束和骨干坐标约束运行relax：（ex flags和使用native；ex flags在enzdes中打开）

-   0.447 sequence recovery (0.077 RMSD) [-557 totalscore]

Running relax with sidechain-sidechain distance constraints at 3.5 distance cutoff. 
>运行松弛时侧链距离限制为3.5距离。

Note that this gets similar RMSD minimization but doesn't maximize sequence recovery or minimize score as well as coordinate constraints. 
>请注意，这得到了类似的RMSD最小化，但并没有最大化序列恢复或最小化分数以及坐标约束。

We tested this protocol with a variety of sidechain-sidechain distance constraint cutoff values and found it to slightly but systematically underperform coordinate constraints:
>我们用各种侧链-侧链距离约束截止值测试了该协议，发现它略微但系统性地低于坐标约束：

-   0.436 sequence recovery (0.0706 RMSD) (-534 totalscore)

Running relax with backbone constraints only: 
>仅对主链（backbone）施加约束运行松弛（relax）

(Note that this does worse in terms of RMSD and that many side-chains are a few angstroms off)
>请注意，就RMSD而言，这更糟糕，许多侧链都偏离了几埃

-   0.488 sequences recovery (0.098 RMSD) (-633 totalscore)

No relax benchmark and Rosetta scoring of native input structures:
>未进行松弛优化的基准测试及天然输入结构的Rosetta评分

-   0.40 (0 RMSD by definition) (-194.7 totalscore)

At this point, the astute reader might ask, what score terms became 438 Rosetta Energy Units (REU) better? 
>此时，敏锐的读者可能会问：究竟是哪些评分项改善了438个Rosetta能量单位（REU）？

We ranked the difference in scores over all structures, comparing the all-atom coordinate constraint protocol and the non-relaxed input structure.
>我们对所有结构的评分差异进行了排序，比较了全原子坐标约束方案与未松弛的输入结构之间的差异。

The biggest difference is fa\_dun (-192.4), followed by fa\_rep (-108.1), pro\_close (-26.4), hbond\_sc (-10.8) and omega (-8.5). 
>能量差异最大的项是 fa_dun（-192.4），其次是 fa_rep（-108.1）、pro_close（-26.4）、hbond_sc（-10.8） 和 omega（-8.5）


|能量项|	说明|	科学意义|
|---|---|---|
|fa_dun|	侧链二面角能量（Dunbrack rotamer能量）|	反映侧链构象优化程度，负值越大说明侧链越接近理想旋转异构体（rotamer）状态
|fa_rep|	原子间范德华排斥力|	负值减少（如-108.1→-200）表示原子冲突降低|
|pro_close|	脯氨酸环闭合能量|	改善脯氨酸的几何构型（特别是环状结构的稳定性）|
|hbond_sc|	侧链参与的氢键能量|	优化极性侧链（如Ser/Thr/Tyr）的氢键网络|
|omega|	肽键平面性（Cα-C-N-Cα二面角）|	维持肽键的平面构象（理想值180°）|

Many input rotamers are close to but not in a "good" Dunbrack rotamer, and the backbone has to be slightly tweaked in order for that residue to get a good dunbrack score. 
>许多输入rotamers接近但不属于“好”的Dunbrack旋转异构体，并且必须对骨架进行轻微调整，才能使残留物获得良好的Dunbrack分数。

Also, many atoms are slightly too close, and they give the fa\_rep contribution.
>此外，许多原子稍微靠得太近，它们会产生fa_rep贡献。

# Developer Discussion: Original Question from Ramesh Jha

Is there a consensus protocol to create the starting PDBs to be used in mini (Rosetta)? 
>在Rosetta的mini（结构最小化）过程中，是否存在用于生成起始PDB文件的共识协议？

It is not unknown that the PDBs right from the Protein Data Bank are composed of artifacts and defects that can give an exceptional jump in energy if happened to be altered during a design protocol. 
>众所周知，直接从蛋白质数据库（Protein Data Bank）获取的PDB文件可能包含人为假象和缺陷，如果在设计流程中这些部分被意外修改，可能会导致能量值异常飙升。

In order to minimize this problem, there are a few things which can be tried and that I am aware of:
>为了尽量减少这个问题，有几件事可以尝试，我知道：

1. Repack (with or without `-use_input_sc`)
2. Repack w/ `sc_min`
3. Relax (fast relax w/ or w/o `use_input_sc`)
4. Idealize

|步骤|	关键参数/操作|	计算目标|	典型应用场景|
|---|---|---|---|
|1.Repack|	-use_input_sc（默认False）|仅优化侧链旋转异构体（rotamer）|-快速评估设计序列的侧链兼容性<br> -对接（docking）前的预处理|                                                                    
|2.Repack+sc_min|-sc_min（侧链局部最小化）|在rotamer采样基础上微调侧链二面角|-修复局部原子冲突<br> -优化氢键网络|                                                                          
|3.Fast Relax|	-relax:fast<br> -use_input_sc |主链+侧链协同优化（受限的能量最小化）|-实验结构的能量优化<br> -设计后的结构验证|         		                        
|4.Idealize|-idealize:skip_chain_closure|强制几何参数标准化（键长/键角/二面角）|-修复AF2预测模型的异常几何<br> -准备分子动力学模拟的初始结构|
                                                                           
                                                                                   
Having tried all of them, I thought the option 3, was the best one, where I used fast relax while using `-use_input_sc` flag. 
>在尝试所有选项后，我认为方案3（启用 -use_input_sc 参数进行快速松弛）是最佳选择。

But recently I observed that though 'relax' is able to substantially decrease the energy of starting PDBs, also result in subtle movements in the backbones and a PDB which could accommodate a ligand could not anymore after being relaxed.
>但最近我发现，尽管"relax"能显著降低初始PDB结构的能量，它也会导致主链的细微移动——以至于一个原本能容纳配体的PDB结构，在松弛后反而无法再结合配体了。
>能量最小化可能牺牲功能相关的非最优构象。建议结合实验数据（如晶体接触分析）或使用约束性协议。

## James's Reply

Try adding the `-constrain_relax_to_start_coords` option to your protocol \#3.

## Ben's Reply

I've been using a protocol that does sc & bb minimization, full packing with `-use_input_sc`, then minimization of bb, rb, and sc. 
It's located in: `rosetta/rosetta\_source/src/apps/pilot/stranges/InterfaceStructMaker.cc`. 
The idea with this is that it keeps things from moving too far from the starting structure. 
There's no backbone sampling so I typically find RMSD to the crystal structure to be \< 1.0. 
Relax actually will do explicit bb sampling thus gives a lower energy structure than my protocol but can also introduce the changes that you observed. 
I'm pasting my typical options file below:

```
-database Rosetta/main/database
-nstruct 20
-ndruns 10
-no_his_his_pairE
-run::min_type dfpmin_armijo_nonmonotone
-ignore_unrecognized_res
-use_input_sc
-ex1
-ex2
-no_optH false
-overwrite
-allow_rbmin true
-min_all_jumps true
-mute protocols.moves.RigidBodyMover protocols.moves.RigidBodyMover core.scoring.etable core.pack.task protocols.docking.DockingInitialPerturbation protocols.TrialMover core.io.database
```

## Sagar's Reply

I just use repack with `-sc_min`, and include the ligand in the process.

## Rocco's Reply

There is a fixed-backbone minimization program that's part of the ligand docking application, ligand\_rpkmin (See section "Preparing the protein receptor for docking" of [[the ligand docking documentation|ligand dock]].

It won't relieve any backbone strain, though, so you may still have issues if the downstream protocol allows for backbone movement.

## Steven C.'s Reply

In general, the Meiler lab does the following:

1. Obtain the protein using a script (see attached python scripts. (scripts were written by I believe James Thompson from the Baker lab and edited by me to work with the Meiler lab configuration) The script cleans, renumbers, and removes multiple conformations of residues from the protein.
2. `relax.linuxgccrelease -ex1 -ex2 -ex1aro -relax:sequence`

This will alleviate clashes in the protein and give a good starting structure for any of the Rosetta applications...unless the protein blows up for some reason.

With an addendum from James Thompson:

Thanks Steven! That's a very reasonable way to gently structures from the PDB, although there are a lot of different ways that you might try this.

Here are two more things that come to mind:

-   If you notice that your protein is moving too much, try adding the `-constrain_relax_to_start_coords` option. This will use coordinate constraints to make your protein stay closer to your input model.
-   Also, Mike Tyka wrote that script for cleaning up PDBs. One of the most useful parts of that script is that it matches non-canonical amino acids (such as selenomethionines) with the appropriate canonical amino acid.

With a second addendum from Andrew Leaver-Fay:

I thought I might point out two things:

1. ex1aro doesn't do anything extra if you already have ex1 on your command line. You can however set the sampling level for ex1aro to be higher than for ex1; e.g. `-ex1 -ex1aro:level 4`. This is stated very explicitly in the option documentation, yet still surprises a lot of people. In Rosetta++, `-ex1aro` behaved as if it were `-ex1aro:level 4`.
2. Mike has observed that extra rotamers will not yield better energy structures out of relax; they will however slow it down.

# Limitations

This is not "normal" relax, that is to say, it will not find the global energy minimum for your structure, it will only find a good energy that is consistent with the input atom positions.

# Post Processing

What post processing steps are typical? 
Are score vs RMSD plots useful? 
Are structures clustered (if so, give a command line)? 
Is it obvious when either the application has succeeded or if it has failed (e.g. if the protocol makes predictions like "This is the docked conformation of proteins A and B"). 
In the case of designs, how should designs be selected?

## See Also

* [[Making Rosetta robust against malformed PDBs|robust]]
* [Controlling Input and Output Tutorial](https://www.rosettacommons.org/demos/latest/tutorials/input_and_output/input_and_output#controlling-input_preparing-a-structure-by-refinement)
* [[Preparing ligands]]: Preparing ligands for use in Rosetta
* [Preparing non-protein residues Tutorial](https://www.rosettacommons.org/demos/wiki/tutorials/prepare_ligand/prepare_ligand_tutorial)
* [[Preparing PDB files for non-peptide polymers]]
* [[Preparing PDB files containing protein and RNA|RNA-protein-changes]]
* [[Running Rosetta with options]]: Instructions for running Rosetta applications on the command line
* [[File types list]]: File types used in Rosetta
