# 翻译规范（必读）

本任务：把 Gholami & Hajiaghayi《Linear Single-Source Shortest Path, Almost》
（arXiv:2611.34911v1）逐节译成中文，并用 Typst 排版。LaTeX 源在
C:\Users\mei\dsh\work\arxiv-2611.34911\src\ ；目标 .typ 文件在本目录。

## 0. 绝对要求

1. 只改你负责的那一个 .typ 文件。不要改 lib.typ、main.typ 或其他人的文件。
2. 保持骨架里已有的标题与 <label> 原样（可以补正文，不要删标题、不要改标签名、不要改顺序）。
   你文件里每个 「#lem[] <lem:x>」 这种占位都要换成真正的完整内容。
3. 忠实、完整：陈述、证明、算法体、图注、枚举列表全部译出，不省略、不概括、不合并。
   公式一律保留（转成 Typst 数学语法）。
4. 交付前整份文档必须能编译，见第 6 节。

## 1. 库函数（lib.typ，文件开头已 #import "lib.typ": *）

| 函数 | 含义 | 用法 |
|---|---|---|
| #thm[...] | 定理（独立编号） | #thm[陈述文本] <thm:main> |
| #lem / #cor / #prop / #defn / #obs / #rem / #claim | 引理/推论/命题/定义/观察/注/断言（共享节内编号） | #lem[文本] <lem:owners>；带小标题：#lem(name: [归属])[文本] |
| #proof[...] | 证明 | #proof[正文]；带标题：#proof(title: [定理 1 的证明])[正文] |
| #alg(...)[...] | 算法框 | #alg(title: [算法名])[#st[...] ...] <alg:sssp> |
| #st[...] | 算法里的一行 | |
| #ind[...] | 算法缩进一层（可嵌套） | |
| #fn("名字", "参数")[...] | 算法里的 Function 头 | #fn("Scan", "X, I, c")[ ... ] |
| #kwin[...] / #kwout[...] | 算法的 Input / Output | |
| #ph[...] | 段内小标题（原文 \paragraph{...}） | 骨架已给出 |
| #sref(<sec:x>) | 引用节 → “第 4.1 节” | 用于原文 \tref{Section}{sec:x} |
| #bref(<par:x>, [文字]) | 引用段内标题 | 用于原文 \hyperref[par:x]{...} |

交叉引用：原文 \tref{Lemma}{lem:owners} → @lem:owners（渲染成“引理 5.3”）；
\tref{Theorem}{thm:main} → @thm:main；\tref{Algorithm}{alg:visit} → @alg:visit；
\tref{Definition}{def:owner} → @def:owner；\tref{Figure}{fig:detour} → @fig:detour。
原文这类引用只写 @标签，不要自己写“定理 1”。

参考文献：键名已改成 Typst 合法形式（去掉 DBLP: 前缀与路径，斜杠换成连字符）：

| 原文 \cite{...} | Typst |
|---|---|
| DBLP:journals/nm/Dijkstra59 | @nm-Dijkstra59 |
| DBLP:journals/jacm/Thorup99 | @jacm-Thorup99 |
| DBLP:journals/jacm/FredmanT87 | @jacm-FredmanT87 |
| DBLP:journals/siamcomp/PettieR05 | @siamcomp-PettieR05 |
| DBLP:conf/stoc/DuanMMSY25 | @stoc-DuanMMSY25 |
| DBLP:conf/focs/DuanMSY23 | @focs-DuanMSY23 |
| DBLP:journals/corr/abs-2602-07868 | @abs-2602-07868 |
| DBLP:journals/corr/abs-2609-15247 | @abs-2609-15247 |
| DBLP:conf/icalp/DuanMSY26 | @icalp-DuanMSY26 |
| DBLP:journals/jcss/FakcharoenpholRT04 | @jcss-FakcharoenpholRT04 |
| DBLP:journals/jacm/Chazelle00 | @jacm-Chazelle00 |
| DBLP:journals/jacm/Tarjan75 | @jacm-Tarjan75 |
| DBLP:journals/tcs/BenderF04 | @tcs-BenderF04 |
| DBLP:journals/siamcomp/HarelT84 | @siamcomp-HarelT84 |
| DBLP:journals/jcss/BlumFPRT73 | @jcss-BlumFPRT73 |
| Yan2026 | @Yan2026 |
| GrishinZabokritskiy2026 | @GrishinZabokritskiy2026 |

原文 \cite[Theorem 37]{key} 这类定位引用，写成 @key 再补文字：“@key 定理 37”。

## 2. LaTeX → Typst 数学语法

- 行内 $x$ → $x$；行间 \[ ... \] 与 align* → 前后各空一行写 $ ... $；
  align* 里的 &= 对齐照抄（Typst 支持 & 对齐点，换行用两个反斜杠）。
- \le → <=，\ge → >=，\ne → !=，\cdot → dot，\times → times，\pm → plus.minus
- \lceil x \rceil → ceil(x)，\lfloor x \rfloor → floor(x)，\lvert x \rvert → abs(x)
- \sum → sum，\prod → prod，\infty → infinity，\varnothing → emptyset
- \gets → arrow.l（或 <-），\to → arrow.r（或 ->），\mapsto → arrow.r.map
- \subset → subset，\subseteq → subset.eq，\in → in，\setminus → without，\cup → union，\cap → sect
- \O(...) → O(...)；\Omega → Omega，\Theta → Theta，\Delta → Delta，\Lambda → Lambda
- 希腊字母（alpha beta gamma delta eta theta kappa lambda mu nu pi rho sigma tau phi psi omega）在 Typst 中同名直接写，大写首字母即大写。
- \hat{x} → hat(x)，\bar{x} → overline(x)，\tilde{x} → tilde(x)
- \operatorname{depth}(v) → op("depth")(v)；\text{abc} → "abc"
- \quad → quad，\qquad → wide（或两个 quad），\; → quad
- \{ \} → { }（Typst 数学里花括号直接写；# 要转义成 \#）
- \ldots → dots，\cdots → dots.c，\Bigl( \Bigr) → lr(( ... ))
- \mathrm{MST}(m,n) → op("MST")(m, n)
- \frac{a}{b} → a/b 或 frac(a, b)
- \begin{cases} → cases(...)
- |x| → abs(x)；\lVert x \rVert → norm(x)
- 下标多字符用 _"max" 或 _(max)：T_{max} → T_max，d_{G}(s,v) → d_G (s, v)
- 反斜杠换行只在行间公式与算法体内少量使用。

## 3. 术语表（必须统一，首次出现可在括号里给英文）

| English | 中文 |
|---|---|
| single-source shortest paths (SSSP) | 单源最短路 |
| comparison–addition model | 比较–加法模型 |
| opaque real register | 不透明实数寄存器 |
| adjacency list / multigraph / loop | 邻接表 / 多重图 / 自环 |
| geodesic | 测地线 |
| tight (arc/path) | 紧（紧弧/紧路径） |
| pivot | 枢纽 |
| owner decomposition / owner cell / owner radius | 归属分解 / 归属单元 / 归属半径 |
| corridor / corridor cut | 走廊 / 走廊割 |
| Steiner branch vertex | 斯坦纳分支顶点 |
| separator size b | 分隔规模 b |
| projected edge | 投影边 |
| request | 请求 |
| certificate / route certificate | 证书 / 路由证书 |
| name (lambda) / dictionary | 名称 / 词典 |
| radial edge / macro edge | 径向边 / 宏观边 |
| auxiliary length / auxiliary graph / auxiliary metric | 辅助长度 / 辅助图 / 辅助度量 |
| owner coordinate J | 归属坐标 J |
| state / state graph / copy | 状态 / 状态图 / 副本 |
| native (inward/outward) copy | 固有（内向/外向）副本 |
| shell / native shell / two-shell window | 壳层 / 固有壳层 / 双壳窗口 |
| near / far transfer | 近邻 / 远距转移 |
| phase switch | 相位切换 |
| reset | 重置 |
| block | 块 |
| ordered hierarchy / seeded ordered hierarchy | 有序层级 / 带种子的有序层级 |
| normalization / unit nu / mass S | 规范化 / 单位 nu / 质量 S |
| light arc / child order / quotient | 轻弧 / 子序 / 商（商图） |
| slot / nominal slot / slot origin | 时间槽 / 名义时间槽 / 槽原点 |
| offer | 候选值 |
| due / deadline / expiry / inclusive | 到期 / 期限 / 过期时刻 / 闭（区间） |
| live / retired / unopened | 活跃 / 已退役 / 未开启 |
| emission / emit | 发出值 / 发出 |
| settle | 确定（settle a vertex：确定某顶点） |
| feasible | 可行 |
| witness / tight witness / protected | 见证 / 紧见证 / 受保护的 |
| run（路径上的） | 段（径向段） |
| run（子节点列表里的连续项） | 连续段 |
| segment | 段 |
| demand | 需求 |
| sample / estimator | 样本 / 估计器 |
| weight band / merge forest | 权重带 / 合并森林 |
| center / partition tree / scale Delta | 中心 / 划分树 / 尺度 Delta |
| table index | 表下标 |
| qualify / qualified | 达标 / 合格的 |
| critical edge | 关键边 |
| halving rule | 折半规则 |
| doubling / dyadic / grid / octave | 倍增 / 二进 / 网格 / 倍频程 |
| atom / cutoff / capped | 原子 / 截断值 / 封顶的 |
| calendar / marker / width class | 日历 / 标记 / 宽度类 |
| family / activation list / cell | 族 / 激活表 / 单元 |
| clip / clipped | 裁剪 / 被裁剪的 |
| strongly connected component (SCC) / condensation | 强连通分量 / 凝聚 |
| join (record) | 连接（记录） |
| union event / merge event | 合并事件 |
| aligned / unaligned | 对齐 / 未对齐 |
| significant weight | 显著权重 |
| replaceable atom / partition hierarchy | 可替换原子 / 剖分层级 |
| native entry | 固有入口 |
| substitution / splice | 替换 / 拼接 |
| spoke（owner spoke） | 辐条（归属辐条） |
| hop bound h | 跳数界 h |
| bottleneck | 瓶颈 |
| growth pass | 生长趟 |
| disjoint set / union–find | 并查集 |
| lowest common ancestor (LCA) / level ancestor | 最近公共祖先 / 层级祖先 |
| split–findmin | 分裂–找最小 |
| Fibonacci heap / decrease-key / extract | 斐波那契堆 / 减键 / 提取 |
| amortized | 摊还 |
| selection（rank-r key） | 选择（第 r 位键） |
| stable counting sort | 稳定计数排序 |
| topological order | 拓扑序 |
| nearly linear | 近似线性 |
| almost surely / expected running time | 几乎必然 / 期望运行时间 |
| reachable / unreachable | 可达 / 不可达 |
| predecessor tree / shortest-path tree | 前驱树 / 最短路树 |
| restore | 恢复 |
| edge index / endpoint | 边下标 / 端点 |
| word（Theta(log n) 位） | 字 |
| rank（并查集秩 / 键名次） | 秩 |

写法约定：中文行文，术语统一，句子通顺，不要逐词硬译；公式变量保持原符号。
\emph → #emph[...]（中文少用斜体，必要时用 #strong[...]）；\textbf → *...*；
\texttt → #raw("...")；\url{x} → #link("x")[x]。

## 4. 算法体写法

原文用 algorithm2e（Function / Input / Output / for each / while / if / else / return）。译法示例：

#alg(title: [把证书分配给枢纽边])[
  #fn("NameRoutes", "G_P")[
    #kwin[如 @def:route-certificate 所述的枢纽路由实例 $G_P$。]
    #kwout[一个已排序的词典、每条输入边对应的词典下标，以及因子 $P$。]
    #st[$u arrow.l ceil(log_2 (t+2))$，$g arrow.l 4u+3$]
    #st[*while* $abs(U) > 32t$：]
    #ind[
      #st[从 $U$ 中均匀随机地取出含 $16t$ 条边的子集 $S$]
      #st[*for each* 边界对应的非空组 $Z$：]
      #ind[
        #st[*if* $a_*$ 存在：]
        #ind[#st[令 $a_e >= a_*$ 的边达标，并取 $lambda_e = a_*$]]
      ]
    ]
    #st[*return* 该词典、下标与 $P$]
  ]
] <alg:name>

注意：算法内的井号一律转义或改写句子避免；内容块里要写方括号时用 #"[" 之类。

## 5. 图

原文两张 TikZ 图已由主控重画为 figs/fig-detour.typ 与 figs/fig-clipped.typ。
你负责的小节若遇到图，在骨架标出的位置写 #include "figs/fig-X.typ" 即可（标签与图注已含在文件里）。

## 6. 编译自检

cd C:\Users\mei\dsh\work\arxiv-2611.34911\typst
typst compile main.typ main.pdf

只看你自己文件（sec-XX-*.typ）里的报错：其他小节由别的代理同时编写，它们的报错可以忽略，
但你自己文件的报错必须清零，且 main.pdf 必须能生成。

常见坑：
- 文本里的井号要写成 \#。
- 行间公式前后必须有空行。
- 标签必须紧跟在元素后：#lem[...] <lem:x>。
- 中文引号直接用“”，英文引号用 \" 转义写法（在内容块里直接写 " 即可）。
