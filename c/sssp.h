/* sssp.h -- 《Linear Single-Source Shortest Path, Almost》的 C 实现
 *
 * 论文：Iman Gholami, MohammadTaghi Hajiaghayi, arXiv:2611.34911v1
 *
 * 实现范围（详见 README.md）：
 *   第 3 节  图与数据结构原语：Dijkstra（扫描版/堆版）、有界生成树生长趟（算法 1）、并查集
 *   第 4 节  有序层级的求值：规范化、子节点排序与间距、时间槽扫描（算法 7、8）
 *   第 5 节  局部层级：归属分解、径向/宏观分类、阈值分量与连接（算法 6 的核心）
 *   第 6 节  枢纽路由证书：采样 + 样本距离估计（算法 2、3）
 *   第 7 节  归属坐标 J（用精确 Dijkstra 求同一个数学对象）
 *   第 8 节  外层状态图、阈值连接、主算法（算法 9）
 */
#ifndef SSSP_H
#define SSSP_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "len.h"

#define LEN_BIG 1e300

/* ------------------------------------------------------------------ */
/* 工具：动态整数数组、随机数、并查集、二叉堆                          */
/* ------------------------------------------------------------------ */

typedef struct { int *a; int n, cap; } ivec;
void ivec_init(ivec *v);
void ivec_push(ivec *v, int x);
void ivec_reserve(ivec *v, int cap);
void ivec_free(ivec *v);
void ivec_clear(ivec *v);

typedef struct { len_t *a; int n, cap; } lvec;
void lvec_init(lvec *v);
void lvec_push(lvec *v, len_t x);
void lvec_free(lvec *v);

void *xmalloc(size_t n);
void *xcalloc(size_t n, size_t s);
void *xrealloc(void *p, size_t n);

void rng_seed(uint64_t s);
uint64_t rng_next(void);
int rng_int(int n);          /* [0, n) */
double rng_double(void);     /* [0,1) */

typedef struct { int *p; int *rk; int n; } DSU;
void dsu_init(DSU *d, int n);
int  dsu_find(DSU *d, int x);
int  dsu_unite(DSU *d, int a, int b);   /* 合并则返回 1 */
void dsu_free(DSU *d);

/* 小顶堆，支持 decrease-key（按项的当前位置） */
typedef struct {
    len_t *key;
    int *item;
    int *pos;
    int n, cap;
} Heap;
void heap_init(Heap *h, int maxn);
void heap_free(Heap *h);
void heap_clear(Heap *h);
int  heap_empty(Heap *h);
int  heap_push(Heap *h, len_t k, int item);
void heap_decrease(Heap *h, int item, len_t k);
int  heap_pop(Heap *h, len_t *out);   /* 返回 item，空则 -1 */

/* ------------------------------------------------------------------ */
/* 图                                                                  */
/* ------------------------------------------------------------------ */

typedef struct {
    int n, m;               /* 顶点数、边数 */
    int *eu, *ev;           /* 边的端点 */
    len_t *ew;              /* 边长度 */
    int *eid;               /* 原始边编号 */
    int *head, *nxt, *to, *arc_e;  /* 邻接：弧 -> 下一条弧，弧 -> 边，弧 -> 头 */
    int narc;
} Graph;

void graph_init(Graph *g);
void graph_free(Graph *g);
int  graph_add_edge(Graph *g, int u, int v, len_t w, int eid);
void graph_finalize(Graph *g);

/* 距离数组 */
typedef struct {
    len_t *d;
    int *par_vertex;
    int *par_edge;
} Dist;
void dist_init(Dist *D, int n);
void dist_free(Dist *D);

/* 精确参照：堆版 Dijkstra（算法无关，用于验证） */
void dijkstra_heap(const Graph *g, int s, Dist *D);
/* 论文 §3.1 的"扫描版"Dijkstra：每次 O(k) 找最小值，共 O(k^2+l) */
void dijkstra_scan(const Graph *g, int s, Dist *D);
/* 半径截断版：只确定距离 <= R 的顶点 */
void dijkstra_radius(const Graph *g, int s, len_t R, Dist *D);

/* ------------------------------------------------------------------ */
/* 树查询（论文 §3.4 的 O(1) LCA/层级祖先由 O(log n) 倍增实现）        */
/* ------------------------------------------------------------------ */

typedef struct {
    int n, root, lg;
    int *depth;
    int *up;          /* [lg][n] */
    len_t *dist_root; /* 到根的距离 */
    len_t *up_max;    /* [lg][n] 到 2^k 祖先路径上的最大边长 */
} Tree;

void  tree_build(Tree *T, int n, int root, const int *tu, const int *tv, const len_t *tw, int nt);
void  tree_free(Tree *T);
int   tree_lca(Tree *T, int u, int v);
len_t tree_distance(Tree *T, int u, int v);   /* D_T(u,v) */
len_t tree_bottleneck(Tree *T, int u, int v); /* b_T(u,v) */

/* ------------------------------------------------------------------ */
/* 算法 1：有界生成树生长趟（三趟得到 MST）                            */
/* ------------------------------------------------------------------ */

int  mst_growth_passes(const Graph *g, int *out_edges);  /* 返回边数 */

/* Kruskal 求 MST（用于局部构造与对照） */
int  mst_kruskal(const Graph *g, int *out_edges);

/* 拓扑排序：返回 0 表示有环 */
int  topo_sort(int n, const ivec *adj, int *order);

/* ------------------------------------------------------------------ */
/* 归属分解（§5.1、§7）                                                */
/* ------------------------------------------------------------------ */

typedef struct {
    int n;                 /* 顶点数 */
    int npiv;              /* 枢纽数 */
    int *piv_vertex;       /* 枢纽 -> 顶点 */
    int *owner;            /* 顶点 -> 枢纽下标 */
    len_t *radius;         /* r_v = b_T(v, p(v)) */
    len_t *tdist;          /* R_v = D_T(v, p(v)) */
    int *cell;             /* 顶点 -> 归属单元 */
    int ncell;
    int *cell_pivot;       /* 单元 -> 枢纽下标 */
    int hop;               /* h = b+1 */
    int sep;               /* 分隔规模 b */
    Tree pstree;           /* 枢纽树 T_A（在枢纽上） */
    int *piv_edge_u, *piv_edge_v; len_t *piv_edge_w; int npiv_edge;
} Owner;

void owner_build(Owner *O, const Graph *g, const int *mst_edges, int nmst, int sep);
void owner_free(Owner *O);
len_t owner_sigma(const Owner *O, int p, int q);   /* 枢纽树上的瓶颈 */

/* 树删除后保留的缩减图：非树边若严格重于其树路径则删除 */
int  reduce_tree_dominated(const Graph *g, const int *mst_edges, int nmst, int *keep_edge);

/* ------------------------------------------------------------------ */
/* 枢纽路由证书（算法 2、3）                                           */
/* ------------------------------------------------------------------ */

typedef struct {
    int npiv;
    int nreq;           /* 投影请求（极小数量的无向多重边） */
    int *ru, *rv;       /* 请求端点（枢纽下标） */
    len_t *rw;          /* 投影权重 A_e */
    int *reps_of_edge;  /* 原始边 -> 请求下标 或 -1 */
    len_t *name;        /* 每条请求的证书 lambda_e */
    int ndict;          /* 词典大小 */
    len_t *dict;        /* 已排序去重的证书值 */
    int *dict_of_req;
    len_t P;            /* 证书因子（2 的幂） */
    int trials, rejects;
} Routes;

void routes_build(Routes *R, int npiv, int nreq, const int *ru, const int *rv,
                  const len_t *rw, uint64_t seed);
void routes_free(Routes *R);
void estimate_sample(int t, const int *su, const int *sv, const len_t *sw,
                     int ns, int h, const int *qu, const int *qv,
                     int *out_index, int *out_disc, int *out_ntab, len_t **out_table);

/* ------------------------------------------------------------------ */
/* 外层状态图（§8.1）                                                  */
/* ------------------------------------------------------------------ */

enum { PH_IN = 0, PH_OUT = 1, PH_TIN = 2, PH_TOUT = 3 };

typedef struct {
    int nstates;
    int *proj;        /* 状态 -> 物理顶点 */
    int *phase;       /* PH_* */
    int *shell;       /* 窗口下标（-1 表示枢轴/特殊状态） */
    int *native_in;   /* 顶点 -> 固有内向副本 */
    int *native_out;  /* 顶点 -> 固有外向副本 */
    /* 状态弧 */
    int *arc_head, *arc_next, *arc_to;
    len_t *arc_len;
    int narc;
    int *vhead;       /* 每个状态的第一条弧 */
    int *arc_macro;   /* 1 = 宏观边对应的弧 */
    int cur_macro;    /* 添加弧时的标记 */
} SGraph;

void sgraph_init(SGraph *S);
void sgraph_free(SGraph *S);
int  sgraph_add_state(SGraph *S, int proj, int phase, int shell);
void sgraph_add_arc(SGraph *S, int u, int v, len_t w);
void sgraph_finalize(SGraph *S);

/* ------------------------------------------------------------------ */
/* 有序层级（定义 3.1、4.6）                                           */
/* ------------------------------------------------------------------ */

typedef struct {
    int nnode;
    int *parent;
    int *child_start, *child_cnt, *children;
    int *nchild;
    len_t *h;         /* 高度 */
    len_t *mass;      /* S(C) */
    len_t *nu;        /* 单位 */
    int *is_leaf;
    int *state_of_leaf;  /* 叶子 -> 状态 */
    int *first_leaf, *last_leaf;
    int nleaf;
    int *leaf_state;     /* 叶布局 -> 状态 */
    int *leaf_node;      /* 叶子节点 -> 叶布局下标 */
    int *child_rank;     /* 子节点在其父中的秩次 */
    len_t *gap;          /* 相邻子节点之间的物理间距 g(i-1,i) */
    int root_cnt;
    int *roots;
    int nstate;       /* 状态数 */
    len_t *extra;     /* 每个节点在本类中参与的连接边长度之和 */          /* 根节点（拓扑序） */
} Hierarchy;

void hier_init(Hierarchy *H);
void hier_free(Hierarchy *H);
int  hier_new_node(Hierarchy *H, int is_leaf);

/* 由连接记录构建层级（阈值合并事件） */
typedef struct { int u, v; len_t w; int su, sv; } JoinRec;
int  hier_build_from_joins(Hierarchy *H, const SGraph *S, int nstate,
                           int naux, const int *aux_state,
                           int njoin, JoinRec *joins, int nroot_out, int **roots_out);

/* 外层状态图的顶点映射（outer.c 内部使用，也供主流程查询） */
typedef struct {
    int nstates;
    int *v_in, *v_inlo, *v_out, *v_outhi;
    int *special_in, *special_out;
    int *shell_of;
    double *base_of;
} StateMap;
void build_outer_states(SGraph *S, StateMap *M, const Graph *g, const Owner *O,
                        const len_t *J, const len_t *de, const int *radial,
                        const int *near, const int *mst_edges, int nmst);
void build_requests(const Graph *g, const Owner *O, int nedge,
                    int **out_ru, int **out_rv, len_t **out_rw, int *out_nreq,
                    int **out_rep, len_t *out_de, int *out_is_macro);
void classify_edges(const Graph *g, const Owner *O, const Routes *R, const int *rep,
                    len_t *out_de, int *out_radial);
void compute_J(const Graph *g, const Owner *O, const len_t *de, const int *radial,
               const int *mst_edges, int nmst, len_t *out_J);

/* 阈值分量 + 连接记录 + 层级（见 hier.c / theta.c） */
int outer_build_hierarchy(Hierarchy *H, const SGraph *S, const Owner *O,
                          int nstates, const int *aux_state, int naux);

/* 规范化（引理 3.1） */
void hier_normalize(Hierarchy *H);
/* 子节点排序与间距（引理 3.2、3.3） */
void hier_order_children(Hierarchy *H, const SGraph *S);

/* ------------------------------------------------------------------ */
/* 求值器（算法 7、8）                                                 */
/* ------------------------------------------------------------------ */

typedef struct {
    Hierarchy *H;
    const SGraph *S;
    len_t *tent;      /* 每个状态的试探值 */
    int *settled;
    int *retired;
    int *out_deg;     /* 松弛输出（收集发出的值） */
    len_t *emit_min;  /* 每个状态发出的最小值 */
    /* 段树：叶布局区间最小值 */
    int *seg_arg;     /* 节点 -> 最小叶下标 */
    int seg_n;
    /* 每个节点的扫描状态 */
    int *opened;
    len_t *clock_origin;
    int *slot_k;
    len_t *expiry;
    int *piece_head;  /* 双向链表：piece id */
    int *cursor_piece;
    /* piece */
    int *pc_kind;     /* 0 = 未注册连续段, 1 = 已注册子节点 */
    int *pc_lo, *pc_hi;  /* 叶区间（RUN） */
    int *pc_child;       /* 子节点（CHILD） */
    int *pc_next, *pc_prev;
    int *pc_owner;
    int npiece, cap_piece;
    /* gap 证书（相邻子节点） */
    int *leaf_pos_of_state;
    len_t *seg_val;
} Evaluator;

void eval_init(Evaluator *E, Hierarchy *H, const SGraph *S);
void eval_free(Evaluator *E);
void eval_seed(Evaluator *E, int state, len_t cost);
void eval_run(Evaluator *E);
len_t eval_distance(const Evaluator *E, int vertex, const int *states_of_vertex, int cnt);
void eval_debug_report(void);
/* 在构造出的状态图上做精确搜索（用于验证 §8.1 的构造是否以原长度实现全部测地线） */
void state_graph_dist(const SGraph *S, int seed_state, len_t *out_state_dist);
void hier_validate(const Hierarchy *H);

/* 主流水线（算法 9） */
typedef struct {
    int sep_log;        /* b = max(2, ceil(log2 N)^sep_log)；64 = 论文设定 */
    int seed;
    int verbose;
    int force_hier;     /* 1 = 即使 N <= L_N^64 也走层级路径 */
    int use_slot_eval;  /* 1 = 用算法 7/8 的槽扫描求值；0 = 在状态图上精确搜索 */
} Options;

void sssp_run(const Graph *g, int source, const Options *opt, Dist *out, void *stats);

#endif /* SSSP_H */