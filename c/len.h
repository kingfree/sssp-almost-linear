/* len.h -- 比较–加法模型中的"不透明实数寄存器"。
 *
 * 论文（§2 预备知识）约定：每条输入长度占据一个不透明实数寄存器；
 * 复制、相加、比较各算一次运算；实数寄存器不能提供地址、整数部分、
 * 对数或任何二进制位。乘以 2^j 需要 j 次加倍（doubling）。
 *
 * 本实现把这些约束固化在类型与 API 里：算法代码只能通过 len_* 函数
 * 访问长度；len_ops 统计"寄存器运算次数"，用来观察算法结构是否真的
 * 只用了比较与加法。
 */
#ifndef LEN_H
#define LEN_H

#include <math.h>
#include <stdint.h>
#include <stdio.h>

typedef struct { double v; } len_t;

extern uint64_t len_ops;        /* 比较 + 加法次数 */
extern uint64_t len_doublings;  /* 加倍次数 */

static inline len_t len_make(double x)        { len_t a; a.v = x; return a; }
static inline len_t len_zero(void)            { return len_make(0.0); }
static inline len_t len_inf(void)             { return len_make(INFINITY); }
static inline int   len_is_inf(len_t a)        { return isinf(a.v); }
static inline int   len_is_zero(len_t a)       { return a.v == 0.0; }
static inline int   len_finite(len_t a)        { return !isinf(a.v); }

static inline len_t len_add(len_t a, len_t b)
{
    len_ops++;
    {
        len_t r;
        r.v = a.v + b.v;
        return r;
    }
}

/* 返回 -1 / 0 / 1 */
static inline int len_cmp(len_t a, len_t b)
{
    len_ops++;
    if (a.v < b.v) return -1;
    if (a.v > b.v) return 1;
    return 0;
}

static inline int len_le(len_t a, len_t b) { return len_cmp(a, b) <= 0; }
static inline int len_lt(len_t a, len_t b) { return len_cmp(a, b) < 0; }
static inline int len_ge(len_t a, len_t b) { return len_cmp(a, b) >= 0; }
static inline int len_gt(len_t a, len_t b) { return len_cmp(a, b) > 0; }
static inline int len_eq(len_t a, len_t b) { return len_cmp(a, b) == 0; }
static inline len_t len_min(len_t a, len_t b) { return len_cmp(a, b) <= 0 ? a : b; }
static inline len_t len_max(len_t a, len_t b) { return len_cmp(a, b) >= 0 ? a : b; }

/* 2^j 倍：j >= 0 时做 j 次加倍；j < 0 时按"共同刻度"约定把值缩小
 * （论文允许在固定 2 的幂分母上把参与比较的值都放到公共刻度上）。 */
static inline len_t len_mul_pow2(len_t a, int j)
{
    if (j >= 0) {
        while (j-- > 0) { a.v *= 2.0; len_doublings++; }
    } else {
        while (j++ < 0) { a.v *= 0.5; len_doublings++; }
    }
    return a;
}

/* 带符号量：存成两个非负和（论文 §2 "signed value stored as two
 * nonnegative sums"）。比较 a-b 与 c-d 时比较 a+d 与 c+b。 */
typedef struct { len_t pos, neg; } sgn_t;

static inline sgn_t sgn_from(len_t a)        { sgn_t s; s.pos = a; s.neg = len_zero(); return s; }
static inline sgn_t sgn_neg(sgn_t a)          { sgn_t s; s.pos = a.neg; s.neg = a.pos; return s; }
static inline sgn_t sgn_add(sgn_t a, sgn_t b) { sgn_t s; s.pos = len_add(a.pos, b.pos); s.neg = len_add(a.neg, b.neg); return s; }
static inline sgn_t sgn_sub(sgn_t a, sgn_t b) { sgn_t s; s.pos = len_add(a.pos, b.neg); s.neg = len_add(a.neg, b.pos); return s; }

static inline int sgn_cmp(sgn_t a, sgn_t b)
{
    /* a.pos - a.neg  与  b.pos - b.neg 比较 <=> a.pos + b.neg 与 b.pos + a.neg 比较 */
    return len_cmp(len_add(a.pos, b.neg), len_add(b.pos, a.neg));
}

#endif /* LEN_H */
