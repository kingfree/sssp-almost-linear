// lib.typ —— 论文中文版排版库（定理环境、算法框、图表、交叉引用）

#let KIND-STMT = "stmt"
#let KIND-THM = "thm"
#let KIND-ALG = "alg"
#let KIND-FIG = "fig"

// ---------- 编号辅助 ----------
// 取出某个位置之前的节号（"4.1" 形式），标题编号为 1.1
// 陈述编号只用一级节号（"4.3"），与原文 lemma 的 [section] 计数一致
#let secnum-at(loc) = {
  let v = counter(heading).at(loc)
  if v.len() == 0 { none } else { str(v.first()) }
}

// 节引用用完整编号（"4.1" 形式）
#let secnum-full-at(loc) = {
  let v = counter(heading).at(loc)
  if v.len() == 0 { none } else { v.map(x => str(x)).join(".") }
}

#let stmt-num(it) = {
  let sec = secnum-at(it.location())
  let n = it.counter.at(it.location()).first()
  if sec == none { str(n) } else { sec + "." + str(n) }
}

#let flat-num(it) = str(it.counter.at(it.location()).first())

// ---------- 定理类环境 ----------
#let thmfigure(body, kind, supp, name) = figure(
  kind: kind, supplement: supp, numbering: n => n,
  caption: name, placement: none, outlined: false,
  body,
)

#let thmrender(it) = block(
  breakable: true, width: 100%, inset: (x: 0pt, y: 0.5em), above: 0.7em, below: 0.7em,
)[
  #set align(left)
  #set par(first-line-indent: (amount: 0em, all: true), justify: true)
  #context strong[#it.supplement #stmt-num(it).]
  #if it.caption != none { [（#it.caption.body）] }
  #h(0.55em)
  #it.body
]

#let thm(body, name: none) = thmfigure(body, KIND-THM, [定理], name)
#let lem(body, name: none) = thmfigure(body, KIND-STMT, [引理], name)
#let cor(body, name: none) = thmfigure(body, KIND-STMT, [推论], name)
#let prop(body, name: none) = thmfigure(body, KIND-STMT, [命题], name)
#let defn(body, name: none) = thmfigure(body, KIND-STMT, [定义], name)
#let obs(body, name: none) = thmfigure(body, KIND-STMT, [观察], name)
#let rem(body, name: none) = thmfigure(body, KIND-STMT, [注], name)
#let claim(body, name: none) = thmfigure(body, KIND-STMT, [断言], name)

#show figure.where(kind: KIND-THM): thmrender
#show figure.where(kind: KIND-STMT): thmrender

// ---------- 证明 ----------
#let proof(body, title: [证明]) = block(
  breakable: true, above: 0.5em, below: 0.7em,
)[
  #set align(left)
  #set par(first-line-indent: (amount: 2em, all: true), justify: true)
  #emph[#title。]
  #body
  #h(1fr) $square$
]

// ---------- 算法 ----------
#let algrender(it) = block(
  breakable: true, width: 100%, inset: (x: 0.2em, y: 0.5em), above: 0.9em, below: 0.9em,
)[
  #set align(left)
  #set par(first-line-indent: (amount: 0em, all: true), justify: false, leading: 1.58em)
  #set block(spacing: 0.5em)
  #set text(size: 10pt)
  #line(length: 100%, stroke: 0.6pt)
  #block(above: 0.35em, below: 0.35em)[
    #context strong[#it.supplement #flat-num(it)：]
    #if it.caption != none { strong[#it.caption.body] }
  ]
  #line(length: 100%, stroke: 0.4pt)
  #block(above: 0.3em)[#it.body]
  #line(length: 100%, stroke: 0.6pt)
]

#let alg(body, title: none) = figure(
  kind: KIND-ALG, supplement: [算法], numbering: n => n,
  caption: title, placement: none, outlined: false, body,
)

#show figure.where(kind: KIND-ALG): algrender

#let st(body) = block(above: 0.34em, below: 0.34em, width: 100%)[#body]
#let ind(body) = block(inset: (left: 1.15em), above: 0.18em, below: 0.1em)[#body]
#let fn(name, args, body) = block(above: 0.45em, below: 0.15em)[
  #block(above: 0.1em, below: 0.5em)[*函数* #name (#args)：]
  #ind[#body]
]
#let kwin(body) = block(above: 0.45em, below: 0.45em)[*输入*：#body]
#let kwout(body) = block(above: 0.45em, below: 0.45em)[*输出*：#body]

// ---------- 图 ----------
#let figrender(it) = block(
  breakable: false, above: 1em, below: 1em, width: 100%,
)[
  #set align(center)
  #it.body
  #v(0.4em)
  #block(
    width: 100%,
  )[
    #set par(first-line-indent: (amount: 0em, all: true), justify: true)
    #set align(left)
    #set text(size: 9.5pt)
    #context [*#it.supplement #flat-num(it).*] #h(0.4em)
    #if it.caption != none { it.caption.body }
  ]
]

#let fig(body, caption: none) = figure(
  kind: KIND-FIG, supplement: [图], numbering: n => n,
  caption: caption, placement: none, outlined: false, body,
)
#show figure.where(kind: KIND-FIG): figrender

// ---------- 交叉引用 ----------
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

// 节引用：`第 4.1 节`
#let sref(lbl) = context {
  let q = query(lbl)
  if q.len() == 0 { [#lbl] } else {
    let el = q.first()
    let n = secnum-full-at(el.location())
    if n == none { link(el.location())[对应章节] } else { link(el.location())[第 #n 节] }
  }
}

// 段落/其它块级引用
#let bref(lbl, body) = context {
  let q = query(lbl)
  if q.len() == 0 { body } else { link(q.first().location(), body) }
}// ---------- 小节内行内标题 ----------
#let ph(body) = block(above: 1.15em, below: 0.8em, sticky: false)[*#body*]