// 图 1：投影边与经认证的替代路径（原文 fig:detour / fig:projection / fig:shortroute）
#import "../lib.typ": *

#let U = 1cm
#let vdot(x, y) = place(dx: x * U - 1.9pt, dy: y * U - 1.9pt, circle(radius: 1.9pt, fill: white, stroke: 0.5pt))
#let vlab(x, y, body) = place(dx: x * U - 0.28em, dy: y * U - 1.55em, text(size: 9.5pt, body))
#let edge(a, b, thickness: 0.7pt, dashed: false) = place(dx: 0cm, dy: 0cm, line(
  start: (a.at(0) * U, a.at(1) * U), end: (b.at(0) * U, b.at(1) * U),
  stroke: (paint: black, thickness: thickness, dash: if dashed { (2.5pt, 1.7pt) } else { none }),
))
#let lab(x, y, body) = place(dx: x * U - 0.28em, dy: y * U - 0.5em, text(size: 9pt, body))

#let panel(body) = block(width: 4.4cm, height: 2.15cm)[#body]

#let panel-a = panel({
  edge((0, 0.45), (0.8, 1.75), thickness: 1.05pt)
  edge((0.8, 1.75), (2.85, 1.75))
  edge((2.85, 1.75), (3.65, 0.45), thickness: 1.05pt)
  vdot(0, 0.45); vdot(0.8, 1.75); vdot(2.85, 1.75); vdot(3.65, 0.45)
  vlab(0, 0.45, [$p$]); vlab(0.8, 1.75, [$x$]); vlab(2.85, 1.75, [$y$]); vlab(3.65, 0.45, [$q$])
  lab(0.02, 1.1, [$R_x$])
  place(dx: 1.75 * U - 0.2em, dy: 1.75 * U + 0.15em, text(size: 9pt)[$w$])
  lab(3.5, 1.1, [$R_y$])
})

#let panel-b = panel({
  edge((0.8, 1.75), (0, 0.45), thickness: 1.05pt)
  edge((0, 0.45), (3.65, 0.45), thickness: 0.8pt, dashed: true)
  edge((3.65, 0.45), (2.85, 1.75), thickness: 1.05pt)
  edge((0.8, 1.75), (2.85, 1.75))
  vdot(0, 0.45); vdot(0.8, 1.75); vdot(2.85, 1.75); vdot(3.65, 0.45)
  vlab(0, 0.45, [$p$]); vlab(0.8, 1.75, [$x$]); vlab(2.85, 1.75, [$y$]); vlab(3.65, 0.45, [$q$])
  lab(0.02, 1.1, [$J_x$])
  place(dx: 1.75 * U - 0.35em, dy: 0.45 * U - 1.5em, text(size: 9pt)[$d_e$])
  lab(3.5, 1.1, [$J_y$])
  place(dx: 1.75 * U - 0.2em, dy: 1.75 * U + 0.15em, text(size: 9pt)[$w$])
})

#fig(caption: [
  枢纽对为 $p = p(x)$、$q = p(y)$。在 (a) 中，树上的两段路径给出输入权重 $R_x + w + R_y$，
  这是枢纽约的一条边的权重。在 (b) 中，上方那条路径表示若干物理游走，其长度由上标出的量所界定。
  径向判定比较 $2d_e$ 与 $w$；随后的一次删除则比较整个上界 $J_x + d_e + J_y$ 与 $w$。
])[
  #grid(
    columns: (1fr, 1fr), gutter: 0.4cm,
    align(center, block[
      #panel-a
      #v(0.2em)
      #align(center)[#text(size: 9.5pt)[(a) 经 $x$ 与 $y$ 的投影边]]
    ]),
    align(center, block[
      #panel-b
      #v(0.2em)
      #align(center)[#text(size: 9.5pt)[(b) 边 $x y$ 的一条经认证替代路径]]
    ]),
  )
] <fig:detour>