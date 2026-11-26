#import "lib.typ": *

// Typst 的 set/show 规则不会随 #import 传播，因此这里必须重复声明渲染规则。
#show figure.where(kind: KIND-STMT): thmrender
#show figure.where(kind: KIND-THM): thmrender
#show figure.where(kind: KIND-ALG): algrender
#show figure.where(kind: KIND-FIG): figrender

#show ref: it => {
  let el = it.element
  if el != none and el.func() == figure and (el.kind == KIND-STMT or el.kind == KIND-THM) {
    context link(el.location())[#el.supplement #stmt-num(el)]
  } else if el != none and el.func() == figure and el.kind == KIND-ALG {
    context link(el.location())[#el.supplement #flat-num(el)]
  } else if el != none and el.func() == figure and el.kind == KIND-FIG {
    context link(el.location())[#el.supplement #flat-num(el)]
  } else { it }
}

#set document(title: "几乎线性的单源最短路", author: ("Iman Gholami", "MohammadTaghi Hajiaghayi"))
#set page(paper: "a4", margin: (x: 2.4cm, y: 2.5cm), numbering: "1")
#set text(font: ("Times New Roman", "Source Han Serif SC"), size: 10.5pt, lang: "zh", region: "cn")
#set par(justify: true, leading: 0.72em, first-line-indent: (amount: 2em, all: true), spacing: 1.05em)

// 段间距、列表间距、标题间距都略大于行距（leading = 0.72em）
#set list(spacing: 0.6em, indent: 1.8em, body-indent: 0.6em)
#set enum(spacing: 0.6em, indent: 1.8em, body-indent: 0.6em)
#show list: set par(first-line-indent: (amount: 0em, all: true))
#show enum: set par(first-line-indent: (amount: 0em, all: true))
#set heading(numbering: "1.1")
#show heading.where(level: 1): it => {
  // 陈述编号按节重新计数（与原文 \newtheorem{lemma}[corollary] 的 [section] 一致）
  counter(figure.where(kind: KIND-STMT)).update(0)
  v(1.5em)
  set text(size: 15pt, weight: "bold")
  block(above: 0.4em, below: 1.05em)[#it]
}
#show heading.where(level: 2): it => {
  v(1.15em)
  set text(size: 12.5pt, weight: "bold")
  block(above: 0.35em, below: 0.75em)[#it]
}
#show heading.where(level: 3): it => {
  v(0.9em)
  set text(size: 11pt, weight: "bold")
  block(above: 0.3em, below: 0.6em)[#it]
}
#set math.equation(numbering: none)
#show raw: set text(font: ("Consolas", "Source Han Sans SC"), size: 9pt)

// ---------- 标题页 ----------
#align(center)[
  #v(3.2cm)
  #text(size: 20pt, weight: "bold")[几乎线性的单源最短路]
  #v(0.6em)
  #text(size: 13pt)[Linear Single-Source Shortest Path, Almost]
  #v(2.2cm)
  #text(size: 12pt)[Iman Gholami \ MohammadTaghi Hajiaghayi]
  #v(0.4em)
  #text(size: 10pt)[马里兰大学（University of Maryland）]
  #v(2.0cm)
  #text(size: 10pt)[arXiv:2611.34911v1 \ 2026 年 11 月 25 日]
  #v(1.6cm)
  #text(size: 9.5pt, fill: luma(70))[中文译稿（Typst 排版）]
]
#pagebreak()

// ---------- 摘要 ----------
#align(center)[#text(size: 12.5pt, weight: "bold")[摘要]]
#v(0.3em)
#block(width: 100%)[
  #set par(first-line-indent: (amount: 0em, all: true), justify: true)
  #set text(size: 10pt)
  我们给出一个用于无向图单源最短路的比较–加法（comparison-addition）算法，其中每条边的长度是任意的非负实数。对每个含 $n$ 个顶点、$m$ 条边的输入，该算法的期望运行时间为 $O(m + n log log log n)$。算法几乎必然终止。每一次终止的运行都返回精确距离，以及可达顶点集合上的一棵最短路树。该时间界与边长之间的比值无关；并且只要 $m = Omega(n log log log n)$，算法就是线性的。
]
#v(1.2em)
#outline(title: [目录], indent: 1.6em, depth: 2)
#pagebreak()

#include "sec-01-intro.typ"
#include "sec-02-prelim.typ"
#include "sec-03-primitives.typ"
#include "sec-03b-structural.typ"
#include "sec-04-eval.typ"
#include "sec-05-local.typ"
#include "sec-06-routing.typ"
#include "sec-07-aux.typ"
#include "sec-08-outer.typ"

// ---------- AI 使用声明 ----------
#heading(level: 1, outlined: false, numbering: none)[AI 使用声明]
OpenAI 的 GPT Astra 6.2 与 GPT Astra 6 协助了部分证明的推演以及撰写与陈述；Gemini Neon 与 GPT Sol 亦协助了写作。结果正确性与论文最终内容的全部责任由作者承担。本工作的 Lean 形式化文件发布在 #link("https://github.com/igholami/sssp-lean")[https://github.com/igholami/sssp-lean]。

// ---------- 译注 ----------
#heading(level: 1, outlined: false, numbering: none)[译者说明]
本文件是 arXiv:2611.34911v1 的中文译稿，用 Typst 重排。定理、引理、定义、算法与图的编号沿用原文（定理与其后各陈述各自编号，节内按出现顺序编号）。跨语言排版中，公式、算法伪代码与图均按原文结构重画；随译稿还附带一份该算法的 C 语言实现，其与本文各构造步骤的对应关系见实现目录下的 README。

#bibliography("references.bib", style: "ieee")