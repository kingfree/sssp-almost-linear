// 图 2：外部宽度单元边界对局部名义时间槽的裁剪（原文 fig:clipped）
#import "../lib.typ": *

#let U = 0.8cm
#let YT = 0.55
#let YB = 1.95
#let hline(y, x0, x1, thickness: 0.7pt) = place(dx: 0cm, dy: 0cm, line(
  start: (x0 * U, y * U), end: (x1 * U, y * U),
  stroke: (paint: black, thickness: thickness),
))
#let tick(x, y) = place(dx: 0cm, dy: 0cm, line(
  start: (x * U, (y - 0.13) * U), end: (x * U, (y + 0.13) * U),
  stroke: (paint: black, thickness: 0.7pt),
))
#let tag(x, y, body) = place(dx: x * U - 0.55em, dy: y * U - 0.42em, text(size: 8.5pt, body))
#let xlabel(x, y, body) = place(dx: x * U - 0.22em, dy: y * U + 0.28em, text(size: 8.5pt, body))

#fig(caption: [
  一个宽度为四的外部边界在时刻 $4$ 切开了局部名义时间槽 $[3, 5)$。
  两次访问依次处理 $[3, 4)$ 与 $[4, 5)$，各自按子节点秩次进行。
  时刻 $4$ 处的一次松弛在第二次访问中即可用。
  在任一次访问中产生的补充候选值都会被填充到至少 $5$。
  每个局部时间槽至多遇到一个这样的切点。
])[
  #block(width: 9.4cm, height: 3.0cm)[
    #hline(YT, 0, 8)
    #for x in (0, 4, 8) { tick(x, YT) }
    #tag(-1.05, YT, [$q = 4$])
    #hline(YB, 0, 8)
    #for x in (0, 1, 2, 3, 4, 5, 6, 7, 8) { tick(x, YB) }
    #tag(-1.05, YB, [$nu = 2$])
    #place(dx: 0cm, dy: 0cm, line(
      start: (3 * U, (YB - 0.22) * U), end: (5 * U, (YB - 0.22) * U),
      stroke: (paint: black, thickness: 1.6pt),
    ))
    #place(dx: 0cm, dy: 0cm, line(
      start: (4 * U, (YT + 0.2) * U), end: (4 * U, (YB - 0.4) * U),
      stroke: (paint: black, thickness: 0.55pt, dash: (2pt, 1.6pt)),
    ))
    #xlabel(0, YB, [$0$]); #xlabel(3, YB, [$3$]); #xlabel(4, YB, [$4$]); #xlabel(5, YB, [$5$]); #xlabel(8, YB, [$8$])
  ]
] <fig:clipped>