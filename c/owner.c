/* owner.c -- 第 5.1 节 归属分解、第 7 节 归属坐标 J
 *
 * 归属分解（§5.1）：
 *   在 MST T 上按后序做 b 规模的标记，得到枢纽集 A0（A0 是最小的、
 *   与每条 b 顶点连通子树相交的集合），补上 A0 生成子树的分支顶点得 A；
 *   对 T-A 的每个"双枢纽邻接"分量，删掉两枢纽走廊上最重的树边；删除后
 *   的连通单元（owner cell）各含一个枢纽。得到 p(v)、R_v、r_v 与枢纽树 T_A。
 *
 * 归属坐标 J（§7）：
 *   把所有枢纽指认成一个公共顶点 a，径向边取辅助长度 w_e + d_e，再保留
 *   归属 MST 辐条；J_v 即 v 到 a 的辅助距离。论文用第 7.3 节的宽度日历 +
 *   局部搜索在近线性时间内算出同一数学对象；这里直接对辅助图跑 Dijkstra，
 *   得到完全相同的 J（差异只在时间复杂度，不影响任何正确性）。
 */
#include "sssp.h"

/* ---------------- 树的辅助结构（根、父、子、后序） ---------------- */

typedef struct {
    int n, root;
    int *parent;
    len_t *wpar;
    int *order;     /* BFS 序 */
    int *pos;       /* 在 order 中的位置 */
    int *head, *nxt, *to;
    int narc;
} RTree;

static void rtree_build(RTree *t, int n, const int *tu, const int *tv, const len_t *tw, int nt, int root)
{
    memset(t, 0, sizeof(*t));
    t->n = n; t->root = root;
    t->parent = (int *)xmalloc((size_t)(n ? n : 1) * sizeof(int));
    t->wpar = (len_t *)xmalloc((size_t)(n ? n : 1) * sizeof(len_t));
    t->order = (int *)xmalloc((size_t)(n ? n : 1) * sizeof(int));
    t->pos = (int *)xmalloc((size_t)(n ? n : 1) * sizeof(int));
    t->head = (int *)xmalloc((size_t)(n ? n : 1) * sizeof(int));
    for (int i = 0; i < n; i++) { t->head[i] = -1; t->parent[i] = -1; }
    t->narc = 0;
    t->nxt = (int *)xmalloc((size_t)(2 * nt + 2) * sizeof(int));
    t->to = (int *)xmalloc((size_t)(2 * nt + 2) * sizeof(int));
    for (int i = 0; i < nt; i++) {
        t->to[t->narc] = tv[i]; t->nxt[t->narc] = t->head[tu[i]]; t->head[tu[i]] = t->narc++;
        t->to[t->narc] = tu[i]; t->nxt[t->narc] = t->head[tv[i]]; t->head[tv[i]] = t->narc++;
    }
    int qh = 0, qt = 0;
    char *vis = (char *)xcalloc((size_t)(n ? n : 1), 1);
    t->order[qt++] = root; vis[root] = 1; t->parent[root] = -1;
    while (qh < qt) {
        int u = t->order[qh++];
        for (int a = t->head[u]; a != -1; a = t->nxt[a]) {
            int v = t->to[a];
            if (vis[v]) continue;
            vis[v] = 1;
            t->parent[v] = u;
            t->order[qt++] = v;
        }
    }
    for (int i = 0; i < qt; i++) t->pos[t->order[i]] = i;
    free(vis);
}
static void rtree_free(RTree *t)
{
    free(t->parent); free(t->wpar); free(t->order); free(t->pos);
    free(t->head); free(t->nxt); free(t->to);
    memset(t, 0, sizeof(*t));
}

/* ---------------- 缩减图：删除被树路径严格统治的非树边 ---------------- */

int reduce_tree_dominated(const Graph *g, const int *mst_edges, int nmst, int *keep_edge)
{
    Tree T;
    int *tu = (int *)xmalloc((size_t)(nmst ? nmst : 1) * sizeof(int));
    int *tv = (int *)xmalloc((size_t)(nmst ? nmst : 1) * sizeof(int));
    len_t *tw = (len_t *)xmalloc((size_t)(nmst ? nmst : 1) * sizeof(len_t));
    for (int i = 0; i < nmst; i++) { tu[i] = g->eu[mst_edges[i]]; tv[i] = g->ev[mst_edges[i]]; tw[i] = g->ew[mst_edges[i]]; }
    tree_build(&T, g->n, 0, tu, tv, tw, nmst);
    char *istree = (char *)xcalloc((size_t)(g->m ? g->m : 1), 1);
    for (int i = 0; i < nmst; i++) istree[mst_edges[i]] = 1;
    for (int e = 0; e < g->m; e++) {
        if (istree[e]) { keep_edge[e] = 1; continue; }
        len_t d = tree_distance(&T, g->eu[e], g->ev[e]);
        keep_edge[e] = !(len_gt(g->ew[e], d));   /* 严格重于树路径则删除 */
    }
    free(istree); free(tu); free(tv); free(tw);
    tree_free(&T);
    return 1;
}

/* ---------------- 归属分解 ---------------- */

static int lca_simple(const Tree *T, int u, int v) { return tree_lca((Tree *)T, u, v); }

void owner_build(Owner *O, const Graph *g, const int *mst_edges, int nmst, int sep)
{
    memset(O, 0, sizeof(*O));
    int n = g->n;
    O->n = n;
    int b = sep < 2 ? 2 : sep;
    O->sep = b;
    O->hop = b + 1;

    int *tu = (int *)xmalloc((size_t)(nmst ? nmst : 1) * sizeof(int));
    int *tv = (int *)xmalloc((size_t)(nmst ? nmst : 1) * sizeof(int));
    len_t *tw = (len_t *)xmalloc((size_t)(nmst ? nmst : 1) * sizeof(len_t));
    for (int i = 0; i < nmst; i++) { tu[i] = g->eu[mst_edges[i]]; tv[i] = g->ev[mst_edges[i]]; tw[i] = g->ew[mst_edges[i]]; }
    Tree T;
    tree_build(&T, n, 0, tu, tv, tw, nmst);
    RTree rt;
    rtree_build(&rt, n, tu, tv, tw, nmst, 0);

    /* 1) A0：后序残量标记 */
    int *resid = (int *)xmalloc((size_t)n * sizeof(int));
    char *inA = (char *)xcalloc((size_t)n, 1);
    ivec pivots; ivec_init(&pivots);
    for (int i = 0; i < n; i++) resid[i] = 1;
    for (int i = n - 1; i >= 0; i--) {
        int v = rt.order[i];
        if (resid[v] >= b) { inA[v] = 1; ivec_push(&pivots, v); resid[v] = 0; }
        if (rt.parent[v] >= 0) resid[rt.parent[v]] += resid[v];
    }
    if (pivots.n == 0) { inA[0] = 1; ivec_push(&pivots, 0); }

    /* 2) A0 生成子树 = 各枢纽到其公共祖先的路径之并 */
    {
        int l = pivots.a[0];
        for (int i = 1; i < pivots.n; i++) l = lca_simple(&T, l, pivots.a[i]);
        char *onpath = (char *)xcalloc((size_t)n, 1);
        for (int i = 0; i < pivots.n; i++) {
            int v = pivots.a[i];
            while (v != l && !onpath[v]) { onpath[v] = 1; v = rt.parent[v]; }
        }
        onpath[l] = 1;
        /* 3) 分支顶点（在生成子树中度 >= 3）加入 A */
        int *deg = (int *)xcalloc((size_t)n, sizeof(int));
        for (int v = 0; v < n; v++) {
            if (!onpath[v]) continue;
            if (rt.parent[v] >= 0 && onpath[rt.parent[v]]) { deg[v]++; deg[rt.parent[v]]++; }
        }
        for (int v = 0; v < n; v++) if (onpath[v] && deg[v] >= 3) inA[v] = 1;
        free(deg); free(onpath);
    }

    /* 4) T - A 的分量，以及每个分量的枢纽邻居 */
    int *comp = (int *)xmalloc((size_t)n * sizeof(int));
    for (int i = 0; i < n; i++) comp[i] = -1;
    int ncomp = 0;
    ivec *pivnb = NULL;
    {
        int *queue = (int *)xmalloc((size_t)n * sizeof(int));
        for (int s = 0; s < n; s++) {
            if (inA[s] || comp[s] >= 0) continue;
            int qh = 0, qt = 0;
            queue[qt++] = s; comp[s] = ncomp;
            while (qh < qt) {
                int u = queue[qh++];
                for (int a = rt.head[u]; a != -1; a = rt.nxt[a]) {
                    int v = rt.to[a];
                    if (inA[v] || comp[v] >= 0) continue;
                    comp[v] = ncomp; queue[qt++] = v;
                }
            }
            ncomp++;
        }
        pivnb = (ivec *)xmalloc((size_t)(ncomp ? ncomp : 1) * sizeof(ivec));
        for (int i = 0; i < ncomp; i++) ivec_init(&pivnb[i]);
        /* 收集相邻的 A 顶点 */
        int *seen = (int *)xcalloc((size_t)n, sizeof(int));
        for (int v = 0; v < n; v++) {
            if (inA[v] || comp[v] < 0) continue;
            for (int a = rt.head[v]; a != -1; a = rt.nxt[a]) {
                int w = rt.to[a];
                if (!inA[w] || seen[w]) continue;
                seen[w] = 1;
                ivec_push(&pivnb[comp[v]], w);
            }
        }
        for (int v = 0; v < n; v++) {
            if (inA[v]) continue;
            for (int a = rt.head[v]; a != -1; a = rt.nxt[a]) seen[rt.to[a]] = 0;
        }
        free(seen); free(queue);
    }

    /* 5) 走廊割：双枢纽邻接分量删去两枢纽走廊上最重的树边 */
    char *deleted = (char *)xcalloc((size_t)n, 1);   /* deleted[v] = 1 表示 v-父 的树边被删 */
    int *del_u = (int *)xmalloc((size_t)(n + 1) * sizeof(int));
    int *del_v = (int *)xmalloc((size_t)(n + 1) * sizeof(int));
    int ndel = 0;
    int *pivmark = (int *)xcalloc((size_t)n, sizeof(int));
    for (int i = 0; i < pivots.n; i++) pivmark[pivots.a[i]] = 1;
    for (int c = 0; c < ncomp; c++) {
        /* 只保留枢纽邻居（A0 成员） */
        int pp[4]; int npp = 0;
        for (int i = 0; i < pivnb[c].n && npp < 4; i++)
            if (pivmark[pivnb[c].a[i]]) pp[npp++] = pivnb[c].a[i];
        if (npp != 2) continue;
        int p = pp[0], q = pp[1];
        int l = lca_simple(&T, p, q);
        /* 沿两端走到 l，找最重的树边 */
        int bestv = -1; len_t bestw = len_zero();
        for (int v = p; v != l; v = rt.parent[v]) {
            len_t w = tree_bottleneck(&T, v, rt.parent[v]);
            if (bestv < 0 || len_gt(w, bestw)) { bestv = v; bestw = w; }
        }
        for (int v = q; v != l; v = rt.parent[v]) {
            len_t w = tree_bottleneck(&T, v, rt.parent[v]);
            if (bestv < 0 || len_gt(w, bestw)) { bestv = v; bestw = w; }
        }
        if (bestv > 0 && !deleted[bestv]) {
            deleted[bestv] = 1;
            del_u[ndel] = rt.parent[bestv]; del_v[ndel] = bestv; ndel++;
        }
    }
    free(pivmark);

    /* 6) 单元 = 删除走廊边后的连通分量 */
    int *cell = (int *)xmalloc((size_t)n * sizeof(int));
    for (int i = 0; i < n; i++) cell[i] = -1;
    int ncell = 0;
    {
        int *queue = (int *)xmalloc((size_t)n * sizeof(int));
        int *inAqueue = (int *)xcalloc((size_t)n, 1);
        for (int s = 0; s < n; s++) {
            if (cell[s] >= 0) continue;
            int qh = 0, qt = 0;
            queue[qt++] = s; cell[s] = ncell;
            while (qh < qt) {
                int u = queue[qh++];
                for (int a = rt.head[u]; a != -1; a = rt.nxt[a]) {
                    int v = rt.to[a];
                    if (cell[v] >= 0) continue;
                    if (deleted[u] && rt.parent[u] == v) continue;
                    if (deleted[v] && rt.parent[v] == u) continue;
                    cell[v] = ncell; queue[qt++] = v;
                }
            }
            ncell++;
        }
        free(queue); free(inAqueue);
    }

    /* 7) 每个单元的枢纽；缺枢纽时退化到最近的枢纽 */
    int *cell_pivot = (int *)xmalloc((size_t)(ncell ? ncell : 1) * sizeof(int));
    for (int i = 0; i < ncell; i++) cell_pivot[i] = -1;
    int badcells = 0;
    for (int i = 0; i < pivots.n; i++) {
        int c = cell[pivots.a[i]];
        if (cell_pivot[c] < 0) cell_pivot[c] = pivots.a[i];
        else if (cell_pivot[c] != pivots.a[i]) badcells++;   /* 同单元多枢纽：不应发生 */
    }
    if (badcells) fprintf(stderr, "[owner] 警告：%d 个单元含多个枢纽（分隔规模 %d）\n", badcells, b);
    /* 无枢纽单元：沿树向上找最近的含枢纽单元 */
    {
        int missing = 0;
        for (int c = 0; c < ncell; c++) if (cell_pivot[c] < 0) missing++;
        if (missing) {
            for (int c = 0; c < ncell; c++) {
                if (cell_pivot[c] >= 0) continue;
                /* 任取该单元一个顶点，沿父链找到有枢纽的单元 */
                int start = -1;
                for (int v = 0; v < n; v++) if (cell[v] == c) { start = v; break; }
                int v = start;
                while (v >= 0 && cell_pivot[cell[v]] < 0) v = rt.parent[v];
                if (v < 0) v = pivots.a[0], cell_pivot[c] = pivots.a[0];
                else cell_pivot[c] = cell_pivot[cell[v]];
            }
            fprintf(stderr, "[owner] 警告：%d 个单元无枢纽，已用最近的枢纽代替\n", missing);
        }
    }

    /* 8) 写回 p(v)、R_v、r_v */
    O->ncell = ncell;
    O->cell = cell;
    O->cell_pivot = cell_pivot;
    O->npiv = pivots.n;
    O->piv_vertex = pivots.a;
    O->owner = (int *)xmalloc((size_t)n * sizeof(int));
    O->radius = (len_t *)xmalloc((size_t)n * sizeof(len_t));
    O->tdist = (len_t *)xmalloc((size_t)n * sizeof(len_t));
    {
        int *pivindex = (int *)xmalloc((size_t)n * sizeof(int));
        for (int i = 0; i < n; i++) pivindex[i] = -1;
        for (int i = 0; i < pivots.n; i++) pivindex[pivots.a[i]] = i;
        for (int v = 0; v < n; v++) {
            int p = cell_pivot[cell[v]];
            O->owner[v] = pivindex[p];
            O->tdist[v] = tree_distance(&T, v, p);
            O->radius[v] = tree_bottleneck(&T, v, p);
        }
        free(pivindex);
    }

    /* 9) 枢纽树 T_A */
    {
        int np = pivots.n;
        int *ptu = (int *)xmalloc((size_t)(np ? np : 1) * sizeof(int));
        int *ptv = (int *)xmalloc((size_t)(np ? np : 1) * sizeof(int));
        len_t *ptw = (len_t *)xmalloc((size_t)(np ? np : 1) * sizeof(len_t));
        int npt = 0;
        int *pivindex = (int *)xmalloc((size_t)n * sizeof(int));
        for (int i = 0; i < n; i++) pivindex[i] = -1;
        for (int i = 0; i < np; i++) pivindex[pivots.a[i]] = i;
        for (int i = 0; i < ndel; i++) {
            int a = del_u[i], c = del_v[i];
            int pa = pivindex[cell_pivot[cell[a]]];
            int pc = pivindex[cell_pivot[cell[c]]];
            if (pa == pc) continue;
            ptu[npt] = pa; ptv[npt] = pc; ptw[npt] = tree_bottleneck(&T, a, c); npt++;
        }
        /* 若枢纽树不连通（退化的分隔规模），用 MST 上的最小瓶颈边补全 */
        if (npt < np - 1) {
            /* 收集所有枢纽对之间的树路径瓶颈，按 Kruskal 补成树 */
            typedef struct { len_t w; int a, b; } PE;
            PE *cands = (PE *)xmalloc((size_t)np * (size_t)np * sizeof(PE));
            int nc = 0;
            for (int i = 0; i < np; i++)
                for (int j = i + 1; j < np; j++) {
                    cands[nc].w = tree_bottleneck(&T, pivots.a[i], pivots.a[j]);
                    cands[nc].a = i; cands[nc].b = j; nc++;
                }
            for (int i = 1; i < nc; i++) {
                PE k = cands[i]; int j = i - 1;
                while (j >= 0 && len_gt(cands[j].w, k.w)) { cands[j + 1] = cands[j]; j--; }
                cands[j + 1] = k;
            }
            DSU d; dsu_init(&d, np);
            for (int i = 0; i < npt; i++) dsu_unite(&d, ptu[i], ptv[i]);
            for (int i = 0; i < nc; i++) {
                if (dsu_unite(&d, cands[i].a, cands[i].b)) {
                    ptu[npt] = cands[i].a; ptv[npt] = cands[i].b; ptw[npt] = cands[i].w; npt++;
                }
            }
            dsu_free(&d);
            free(cands);
        }
        O->piv_edge_u = ptu; O->piv_edge_v = ptv; O->piv_edge_w = ptw; O->npiv_edge = npt;
        tree_build(&O->pstree, np, 0, ptu, ptv, ptw, npt);
        free(pivindex);
    }

    free(del_u); free(del_v); free(deleted);
    for (int i = 0; i < ncomp; i++) ivec_free(&pivnb[i]);
    free(pivnb); free(comp); free(resid); free(inA);
    free(tu); free(tv); free(tw);
    tree_free(&T);
    rtree_free(&rt);
}

void owner_free(Owner *O)
{
    free(O->piv_vertex); free(O->owner); free(O->radius); free(O->tdist);
    free(O->cell); free(O->cell_pivot);
    free(O->piv_edge_u); free(O->piv_edge_v); free(O->piv_edge_w);
    tree_free(&O->pstree);
    memset(O, 0, sizeof(*O));
}

len_t owner_sigma(const Owner *O, int p, int q)
{
    if (p == q) return len_zero();
    return tree_bottleneck((Tree *)&O->pstree, p, q);
}

/* ---------------- 投影请求与归属坐标 J ---------------- */

/* 由投影边构造枢纽请求（每条枢纽对取最轻的投影权重），并写回"边 -> 请求"映射 */
void build_requests(const Graph *g, const Owner *O, int nedge,
                    int **out_ru, int **out_rv, len_t **out_rw, int *out_nreq,
                    int **out_rep, len_t *out_de, int *out_is_macro)
{
    int np = O->npiv;
    int *ru = (int *)xmalloc((size_t)(g->m ? g->m : 1) * sizeof(int));
    int *rv = (int *)xmalloc((size_t)(g->m ? g->m : 1) * sizeof(int));
    len_t *rw = (len_t *)xmalloc((size_t)(g->m ? g->m : 1) * sizeof(len_t));
    int *rep = (int *)xmalloc((size_t)(g->m ? g->m : 1) * sizeof(int));
    for (int i = 0; i < g->m; i++) { ru[i] = rv[i] = -1; rw[i] = len_inf(); rep[i] = -1; }

    /* 分组：按 (p,q) 排序取最小 */
    typedef struct { int a, b; len_t w; int e; } PR;
    int m = g->m;
    PR *pr = (PR *)xmalloc((size_t)(m ? m : 1) * sizeof(PR));
    int npr = 0;
    for (int e = 0; e < m; e++) {
        if (!nedge) break;
        int x = g->eu[e], y = g->ev[e];
        int pi = O->owner[x], qi = O->owner[y];
        if (pi == qi) continue;
        pr[npr].a = pi; pr[npr].b = qi;
        pr[npr].w = len_add(len_add(O->tdist[x], g->ew[e]), O->tdist[y]);
        pr[npr].e = e;
        npr++;
    }
    for (int i = 1; i < npr; i++) {
        PR k = pr[i]; int j = i - 1;
        while (j >= 0) {
            int c = (pr[j].a != k.a) ? (pr[j].a - k.a)
                  : (pr[j].b != k.b) ? (pr[j].b - k.b) : len_cmp(pr[j].w, k.w);
            if (c <= 0) break;
            pr[j + 1] = pr[j]; j--;
        }
        pr[j + 1] = k;
    }
    int nreq = 0;
    {
        int i = 0;
        while (i < npr) {
            int j = i, best = i;
            while (j < npr && pr[j].a == pr[i].a && pr[j].b == pr[i].b) {
                if (len_lt(pr[j].w, pr[best].w)) best = j;
                j++;
            }
            ru[nreq] = pr[best].a; rv[nreq] = pr[best].b; rw[nreq] = pr[best].w;
            for (int k = i; k < j; k++) rep[pr[k].e] = nreq;
            nreq++;
            i = j;
        }
    }
    free(pr);
    *out_ru = ru; *out_rv = rv; *out_rw = rw; *out_rep = rep; *out_nreq = nreq;
    (void)np; (void)out_de; (void)out_is_macro;
}

/* 每条边的证书 d_e 与径向/宏观分类 */
void classify_edges(const Graph *g, const Owner *O, const Routes *R, const int *rep,
                     len_t *out_de, int *out_radial)
{
    for (int e = 0; e < g->m; e++) {
        int x = g->eu[e], y = g->ev[e];
        int pi = O->owner[x], qi = O->owner[y];
        if (pi == qi) { out_de[e] = len_zero(); out_radial[e] = 1; continue; }  /* 单元内：径向，d_e = 0 */
        int r = rep[e];
        len_t lam = (r >= 0 && R) ? R->name[r] : len_inf();
        if (!len_finite(lam)) { out_de[e] = len_inf(); out_radial[e] = 0; continue; }
        len_t d = len_mul_pow2(lam, R->dict ? 0 : 0);
        /* d_e = P * lambda_e，用倍增实现 */
        d = lam;
        {
            double Pv = R->P.v;
            while (Pv > 1.0000001) { d = len_mul_pow2(d, 1); Pv /= 2.0; }
        }
        out_de[e] = d;
        /* 径向 iff 2 d_e <= w */
        out_radial[e] = len_le(len_mul_pow2(d, 1), g->ew[e]) ? 1 : 0;
    }
}

/* 归属坐标：辅助图上从公共根 a 出发的 Dijkstra */
void compute_J(const Graph *g, const Owner *O, const len_t *de, const int *radial,
               const int *mst_edges, int nmst, len_t *out_J)
{
    int n = g->n;
    int na = n + 1;
    Graph aux; graph_init(&aux); aux.n = na;
    int eid = 0;
    for (int e = 0; e < g->m; e++) {
        if (!radial[e]) continue;
        if (!len_finite(de[e])) continue;
        len_t w = len_add(g->ew[e], de[e]);
        graph_add_edge(&aux, g->eu[e], g->ev[e], w, eid++);
    }
    /* 归属 MST 辐条 */
    for (int i = 0; i < nmst; i++) {
        int e = mst_edges[i];
        graph_add_edge(&aux, g->eu[e], g->ev[e], g->ew[e], eid++);
    }
    /* 所有枢纽指认成一个顶点 a */
    for (int i = 0; i < O->npiv; i++) graph_add_edge(&aux, O->piv_vertex[i], n, len_zero(), eid++);
    graph_finalize(&aux);
    Dist D; dist_init(&D, na);
    dijkstra_heap(&aux, n, &D);
    for (int v = 0; v < n; v++) out_J[v] = D.d[v];
    for (int i = 0; i < O->npiv; i++) out_J[O->piv_vertex[i]] = len_zero();
    dist_free(&D);
    graph_free(&aux);
}