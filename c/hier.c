/* hier.c -- 阈值分量、连接记录与有序层级的组装
 *
 * 1) 连接记录：块 MST 边（物理权重）+ 每个状态的 theta 连接到其归属枢纽的
 *    辅助顶点 + 枢纽走廊连接（sigma），见引理 7.3；
 * 2) 按倍增网格的类 [a,2a) 同时合并，得到高为 a 的内部节点（引理 5.2 的
 *    阈值强连通分量刻画）；辅助顶点只识别核，最后被压缩，单子节点被收缩；
 * 3) 质量 S(C)（引理 3.1）、单位 nu(C)（在本构造里等于高度，已在倍增网格上）、
 *    子节点拓扑序、以及相邻子节点代表点的物理间距 g(i,j)（沿连接森林的真实距离）。
 */
#include "sssp.h"

void hier_init(Hierarchy *H) { memset(H, 0, sizeof(*H)); }
void hier_free(Hierarchy *H)
{
    free(H->parent); free(H->child_start); free(H->child_cnt); free(H->children);
    free(H->nchild); free(H->h); free(H->mass); free(H->nu); free(H->is_leaf);
    free(H->state_of_leaf); free(H->first_leaf); free(H->last_leaf);
    free(H->leaf_state); free(H->leaf_node); free(H->child_rank); free(H->gap);
    free(H->roots); free(H->extra);
    memset(H, 0, sizeof(*H));
}
int hier_new_node(Hierarchy *H, int is_leaf)
{
    int id = H->nnode++;
    H->parent = (int *)xrealloc(H->parent, (size_t)H->nnode * sizeof(int));
    H->nchild = (int *)xrealloc(H->nchild, (size_t)H->nnode * sizeof(int));
    H->child_start = (int *)xrealloc(H->child_start, (size_t)H->nnode * sizeof(int));
    H->child_cnt = (int *)xrealloc(H->child_cnt, (size_t)H->nnode * sizeof(int));
    H->h = (len_t *)xrealloc(H->h, (size_t)H->nnode * sizeof(len_t));
    H->mass = (len_t *)xrealloc(H->mass, (size_t)H->nnode * sizeof(len_t));
    H->nu = (len_t *)xrealloc(H->nu, (size_t)H->nnode * sizeof(len_t));
    H->is_leaf = (int *)xrealloc(H->is_leaf, (size_t)H->nnode * sizeof(int));
    H->state_of_leaf = (int *)xrealloc(H->state_of_leaf, (size_t)H->nnode * sizeof(int));
    H->first_leaf = (int *)xrealloc(H->first_leaf, (size_t)H->nnode * sizeof(int));
    H->last_leaf = (int *)xrealloc(H->last_leaf, (size_t)H->nnode * sizeof(int));
    H->child_rank = (int *)xrealloc(H->child_rank, (size_t)H->nnode * sizeof(int));
    H->parent[id] = -1; H->nchild[id] = 0; H->child_start[id] = 0; H->child_cnt[id] = 0;
    H->h[id] = len_zero(); H->mass[id] = len_zero(); H->nu[id] = len_make(1.0);
    H->is_leaf[id] = is_leaf; H->state_of_leaf[id] = -1;
    H->first_leaf[id] = -1; H->last_leaf[id] = -1; H->child_rank[id] = 0;
    return id;
}

static void sort_joins(JoinRec *j, int n)
{
    for (int i = 1; i < n; i++) {
        JoinRec k = j[i]; int p = i - 1;
        while (p >= 0 && len_gt(j[p].w, k.w)) { j[p + 1] = j[p]; p--; }
        j[p + 1] = k;
    }
}

/* 权重在网格上相对 base 的指数 floor(log2(w/base))（只做加倍的比较） */
static int grid_class(len_t w, len_t base)
{
    len_t lim = base;
    int m = 0;
    while (m < 400 && len_ge(w, len_mul_pow2(lim, 1))) { lim = len_mul_pow2(lim, 1); m++; }
    return m;
}

/* ---------------- 连接森林（用于子节点间距） ---------------- */

typedef struct {
    int n;
    int *head, *nxt, *to;
    len_t *w;
    int narc;
} Forest;

static void forest_init(Forest *f, int n, int nedges)
{
    f->n = n;
    f->head = (int *)xmalloc((size_t)(n ? n : 1) * sizeof(int));
    for (int i = 0; i < n; i++) f->head[i] = -1;
    f->nxt = (int *)xmalloc((size_t)(2 * nedges + 2) * sizeof(int));
    f->to = (int *)xmalloc((size_t)(2 * nedges + 2) * sizeof(int));
    f->w = (len_t *)xmalloc((size_t)(2 * nedges + 2) * sizeof(len_t));
    f->narc = 0;
}
static void forest_add(Forest *f, int u, int v, len_t w)
{
    f->to[f->narc] = v; f->w[f->narc] = w; f->nxt[f->narc] = f->head[u]; f->head[u] = f->narc++;
    f->to[f->narc] = u; f->w[f->narc] = w; f->nxt[f->narc] = f->head[v]; f->head[v] = f->narc++;
}
static void forest_free(Forest *f) { free(f->head); free(f->nxt); free(f->to); free(f->w); memset(f, 0, sizeof(*f)); }


/* 叶布局（DFS 序）+ 区间 + 质量。子节点顺序改变后必须重跑。 */
static void hier_layout(Hierarchy *H)
{
    int *stack = (int *)xmalloc((size_t)(H->nnode + 1) * sizeof(int));
    int sp = 0;
    H->nleaf = 0;
    if (!H->leaf_state) H->leaf_state = (int *)xmalloc((size_t)(H->nstate ? H->nstate : 1) * sizeof(int));
    if (!H->leaf_node) H->leaf_node = (int *)xmalloc((size_t)(H->nnode ? H->nnode : 1) * sizeof(int));
    for (int i = 0; i < H->nnode; i++) H->leaf_node[i] = -1;
    for (int i = 0; i < H->nnode; i++) H->first_leaf[i] = H->last_leaf[i] = -1;
    for (int r = H->root_cnt - 1; r >= 0; r--) stack[sp++] = H->roots[r];
    while (sp > 0) {
        int nd = stack[--sp];
        if (H->is_leaf[nd]) {
            H->leaf_node[nd] = H->nleaf;
            H->leaf_state[H->nleaf] = H->state_of_leaf[nd];
            H->nleaf++;
            continue;
        }
        for (int k = H->child_cnt[nd] - 1; k >= 0; k--) {
            int c = H->children[H->child_start[nd] + k];
            if (sp < H->nnode + 1) stack[sp++] = c;
        }
    }
    for (int i = 0; i < H->nnode; i++)
        if (H->is_leaf[i]) H->first_leaf[i] = H->last_leaf[i] = H->leaf_node[i];
    for (int i = 0; i < H->nnode; i++) {
        if (H->is_leaf[i]) { H->mass[i] = len_zero(); continue; }
        int lo = -1, hi2 = -1;
        len_t sm = H->h[i];
        if (H->extra) sm = len_add(sm, H->extra[i]);
        for (int k = 0; k < H->child_cnt[i]; k++) {
            int c = H->children[H->child_start[i] + k];
            if (lo < 0 || H->first_leaf[c] < lo) lo = H->first_leaf[c];
            if (H->last_leaf[c] > hi2) hi2 = H->last_leaf[c];
            sm = len_add(sm, H->mass[c]);
        }
        H->first_leaf[i] = lo; H->last_leaf[i] = hi2; H->mass[i] = sm;
    }
    free(stack);
}

/* ---------------- 主构造 ---------------- */

int hier_build_from_joins(Hierarchy *H, const SGraph *S, int nstate, int naux,
                          const int *aux_state, int njoin, JoinRec *joins,
                          int nroot_out, int **roots_out)
{
    int total = nstate + naux;
    sort_joins(joins, njoin);
    hier_init(H);

    /* ---- 1) 合并树 ---- */
    ivec *mchild = (ivec *)xmalloc((size_t)(njoin + 2) * sizeof(ivec));
    len_t *mh = (len_t *)xmalloc((size_t)(njoin + 2) * sizeof(len_t));
    len_t *msum = (len_t *)xmalloc((size_t)(njoin + 2) * sizeof(len_t));
    for (int i = 0; i < njoin + 2; i++) { ivec_init(&mchild[i]); mh[i] = len_zero(); msum[i] = len_zero(); }
    int nm = 0;
    DSU d; dsu_init(&d, total);
    int *comp_item = (int *)xmalloc((size_t)(total ? total : 1) * sizeof(int));
    for (int i = 0; i < total; i++) comp_item[i] = i;

    int p = 0;
    while (p < njoin) {
        len_t base = joins[p].w;
        int cls = grid_class(joins[p].w, base);
        int q = p;
        while (q < njoin && grid_class(joins[q].w, base) == cls) q++;
        len_t hcls = len_mul_pow2(base, cls);

        /* 本类中出现的当前分量根 */
        int *lmap = (int *)xmalloc((size_t)(total ? total : 1) * sizeof(int));
        for (int i = 0; i < total; i++) lmap[i] = -1;
        ivec roots_local; ivec_init(&roots_local);
        for (int i = p; i < q; i++) {
            int ru = dsu_find(&d, joins[i].u), rv = dsu_find(&d, joins[i].v);
            if (ru == rv) continue;
            if (lmap[ru] < 0) { lmap[ru] = roots_local.n; ivec_push(&roots_local, ru); }
            if (lmap[rv] < 0) { lmap[rv] = roots_local.n; ivec_push(&roots_local, rv); }
        }
        int nl = roots_local.n;
        if (nl > 0) {
            DSU ld; dsu_init(&ld, nl);
            for (int i = p; i < q; i++) {
                int ru = dsu_find(&d, joins[i].u), rv = dsu_find(&d, joins[i].v);
                if (ru == rv) continue;
                dsu_unite(&ld, lmap[ru], lmap[rv]);
            }
            len_t *lsum = (len_t *)xmalloc((size_t)nl * sizeof(len_t));
            for (int i = 0; i < nl; i++) lsum[i] = len_zero();
            for (int i = p; i < q; i++) {
                int ru = dsu_find(&d, joins[i].u), rv = dsu_find(&d, joins[i].v);
                if (ru == rv) continue;
                int r = dsu_find(&ld, lmap[ru]);
                lsum[r] = len_add(lsum[r], joins[i].w);
            }
            ivec *mem = (ivec *)xmalloc((size_t)nl * sizeof(ivec));
            for (int i = 0; i < nl; i++) ivec_init(&mem[i]);
            for (int i = 0; i < nl; i++) {
                int r = dsu_find(&ld, i);
                ivec_push(&mem[r], comp_item[roots_local.a[i]]);
            }
            for (int r = 0; r < nl; r++) {
                if (mem[r].n < 2) continue;
                int nid = nm++;
                mh[nid] = hcls;
                msum[nid] = lsum[r];
                for (int k = 0; k < mem[r].n; k++) ivec_push(&mchild[nid], mem[r].a[k]);
                /* 全局合并成一项，代表用 total + nid 表示 */
                int rep = -1;
                for (int i = 0; i < nl; i++) if (dsu_find(&ld, i) == r) { rep = i; break; }
                int repitem = roots_local.a[rep];
                for (int i = 0; i < nl; i++)
                    if (i != rep && dsu_find(&ld, i) == r) dsu_unite(&d, repitem, roots_local.a[i]);
                comp_item[dsu_find(&d, repitem)] = total + nid;
            }
            for (int i = 0; i < nl; i++) ivec_free(&mem[i]);
            free(mem); free(lsum);
            dsu_free(&ld);
        }
        free(lmap);
        ivec_free(&roots_local);
        p = q;
    }

    /* ---- 2) 展平为只以状态为叶的层级 ---- */
    int *item_node = (int *)xmalloc((size_t)(total ? total : 1) * sizeof(int));
    for (int i = 0; i < nstate; i++) {
        int nd = hier_new_node(H, 1);
        H->state_of_leaf[nd] = i;
        item_node[i] = nd;
    }
    for (int i = nstate; i < total; i++) item_node[i] = -1;
    int *mnode = (int *)xmalloc((size_t)(nm ? nm : 1) * sizeof(int));
    for (int i = 0; i < nm; i++) mnode[i] = -1;

    H->extra = (len_t *)xmalloc((size_t)(total * 2 + 16) * sizeof(len_t));
    for (int i = 0; i < total * 2 + 16; i++) H->extra[i] = len_zero();
    int extra_cap = total * 2 + 16;
    int *charr = NULL; int nch = 0, capch = 0;
    for (int i = 0; i < nm; i++) {
        ivec kids; ivec_init(&kids);
        for (int k = 0; k < mchild[i].n; k++) {
            int it = mchild[i].a[k];
            int nd = (it < total) ? item_node[it] : mnode[it - total];
            if (nd < 0) continue;
            int dup = 0;
            for (int z = 0; z < kids.n; z++) if (kids.a[z] == nd) { dup = 1; break; }
            if (!dup) ivec_push(&kids, nd);
        }
        if (kids.n == 0) mnode[i] = -1;
        else if (kids.n == 1) mnode[i] = kids.a[0];
        else {
            int nd = hier_new_node(H, 0);
            if (nd >= extra_cap) {
                int old = extra_cap;
                extra_cap = nd * 2 + 16;
                H->extra = (len_t *)xrealloc(H->extra, (size_t)extra_cap * sizeof(len_t));
                for (int t = old; t < extra_cap; t++) H->extra[t] = len_zero();
            }
            H->h[nd] = mh[i];
            H->nu[nd] = mh[i];
            H->extra[nd] = msum[i];
            H->child_cnt[nd] = kids.n;
            H->nchild[nd] = kids.n;
            H->child_start[nd] = nch;
            for (int k = 0; k < kids.n; k++) {
                if (nch == capch) { capch = capch ? capch * 2 : 16; charr = (int *)xrealloc(charr, (size_t)capch * sizeof(int)); }
                charr[nch++] = kids.a[k];
                H->parent[kids.a[k]] = nd;
            }
            mnode[i] = nd;
        }
        ivec_free(&kids);
    }
    H->children = (int *)xmalloc((size_t)(nch ? nch : 1) * sizeof(int));
    memcpy(H->children, charr, (size_t)nch * sizeof(int));
    free(charr);

    /* ---- 3) 根（先给任意顺序，供第一次布局） ---- */
    ivec roots; ivec_init(&roots);
    for (int i = 0; i < H->nnode; i++) if (H->parent[i] < 0) ivec_push(&roots, i);
    H->root_cnt = roots.n;
    H->roots = (int *)xmalloc((size_t)(roots.n ? roots.n : 1) * sizeof(int));
    for (int i = 0; i < roots.n; i++) H->roots[i] = roots.a[i];
    H->nstate = nstate;
    hier_layout(H);

    /* ---- 4) 连接森林（Kruskal 取连接记录的子集） ---- */
    Forest F;
    {
        forest_init(&F, nstate, njoin);
        DSU fd; dsu_init(&fd, nstate);
        for (int i = 0; i < njoin; i++) {
            int u = joins[i].su, v = joins[i].sv;
            if (u < 0 || v < 0 || u >= nstate || v >= nstate || u == v) continue;
            if (dsu_unite(&fd, u, v)) forest_add(&F, u, v, joins[i].w);
        }
        dsu_free(&fd);
    }
    Tree FT;
    {
        int nv = nstate + 1;
        int cap = F.narc / 2 + nstate + 4;
        int *tu = (int *)xmalloc((size_t)cap * sizeof(int));
        int *tv = (int *)xmalloc((size_t)cap * sizeof(int));
        len_t *tw = (len_t *)xmalloc((size_t)cap * sizeof(len_t));
        int nt = 0;
        char *vis = (char *)xcalloc((size_t)(nstate ? nstate : 1), 1);
        int *queue = (int *)xmalloc((size_t)(nstate ? nstate : 1) * sizeof(int));
        for (int u = 0; u < nstate; u++) {
            if (vis[u]) continue;
            int qh = 0, qt = 0;
            queue[qt++] = u; vis[u] = 1;
            while (qh < qt) {
                int x = queue[qh++];
                for (int e = F.head[x]; e != -1; e = F.nxt[e]) {
                    int y = F.to[e];
                    if (vis[y]) continue;
                    vis[y] = 1;
                    tu[nt] = x; tv[nt] = y; tw[nt] = F.w[e]; nt++;
                    queue[qt++] = y;
                }
            }
            if (nt < cap) { tu[nt] = u; tv[nt] = nstate; tw[nt] = len_zero(); nt++; }
        }
        tree_build(&FT, nv, nstate, tu, tv, tw, nt);
        free(tu); free(tv); free(tw); free(vis); free(queue);
    }

    /* ---- 5) 弧 -> 最近公共祖先；子节点按轻弧拓扑排序 ---- */
    {
        int *leaf_of_state = (int *)xmalloc((size_t)(nstate ? nstate : 1) * sizeof(int));
        for (int i = 0; i < nstate; i++) leaf_of_state[i] = -1;
        for (int i = 0; i < H->nnode; i++) if (H->is_leaf[i]) leaf_of_state[H->state_of_leaf[i]] = i;
        for (int i = 0; i < H->nnode; i++) {
            if (H->is_leaf[i]) continue;
            int d = H->child_cnt[i];
            if (d < 2) continue;
            int *adj = (int *)xcalloc((size_t)(d * d), sizeof(int));
            int *indeg = (int *)xcalloc((size_t)d, sizeof(int));
            for (int a = 0; a < S->narc; a++) {
                if (len_ge(S->arc_len[a], H->nu[i])) continue;
                int lu = leaf_of_state[S->arc_head[a]], lv = leaf_of_state[S->arc_to[a]];
                if (lu < 0 || lv < 0) continue;
                int pu = H->leaf_node[lu], pv = H->leaf_node[lv];
                int cu = -1, cv = -1;
                for (int k = 0; k < d; k++) {
                    int c = H->children[H->child_start[i] + k];
                    if (H->first_leaf[c] <= pu && pu <= H->last_leaf[c]) cu = k;
                    if (H->first_leaf[c] <= pv && pv <= H->last_leaf[c]) cv = k;
                }
                if (cu < 0 || cv < 0 || cu == cv) continue;
                if (!adj[cu * d + cv]) { adj[cu * d + cv] = 1; indeg[cv]++; }
            }
            int *ord = (int *)xmalloc((size_t)d * sizeof(int));
            int qh = 0, qt = 0;
            for (int k = 0; k < d; k++) if (indeg[k] == 0) ord[qt++] = k;
            while (qh < qt) {
                int u = ord[qh++];
                for (int k = 0; k < d; k++) if (adj[u * d + k]) { if (--indeg[k] == 0) ord[qt++] = k; }
            }
            if (qt < d) for (int k = 0; k < d; k++) ord[k] = k;   /* 有环：退化为原顺序 */
            int *nk = (int *)xmalloc((size_t)d * sizeof(int));
            for (int k = 0; k < d; k++) nk[k] = H->children[H->child_start[i] + ord[k]];
            for (int k = 0; k < d; k++) H->children[H->child_start[i] + k] = nk[k];
            free(nk); free(adj); free(indeg); free(ord);
        }
        free(leaf_of_state);
    }

    /* 子节点顺序已变：重算布局，并写回秩次 */
    hier_layout(H);
    for (int i = 0; i < H->nnode; i++) {
        if (H->is_leaf[i]) continue;
        for (int k = 0; k < H->child_cnt[i]; k++) H->child_rank[H->children[H->child_start[i] + k]] = k;
    }

    /* ---- 6) 根的凝聚序（根之间的弧非循环） ---- */
    {
        int nr = H->root_cnt;
        int *adj = (int *)xcalloc((size_t)(nr * nr > 0 ? nr * nr : 1), sizeof(int));
        int *indeg = (int *)xcalloc((size_t)(nr ? nr : 1), sizeof(int));
        int *leaf_root = (int *)xmalloc((size_t)(H->nleaf ? H->nleaf : 1) * sizeof(int));
        for (int r = 0; r < nr; r++)
            for (int lp = H->first_leaf[H->roots[r]]; lp <= H->last_leaf[H->roots[r]]; lp++) leaf_root[lp] = r;
        for (int a = 0; a < S->narc; a++) {
            int lu = -1, lv = -1;
            for (int i = 0; i < H->nnode; i++) { if (H->is_leaf[i] && H->state_of_leaf[i] == S->arc_head[a]) { lu = i; break; } }
            for (int i = 0; i < H->nnode; i++) { if (H->is_leaf[i] && H->state_of_leaf[i] == S->arc_to[a]) { lv = i; break; } }
            if (lu < 0 || lv < 0) continue;
            int ru = leaf_root[H->leaf_node[lu]], rv = leaf_root[H->leaf_node[lv]];
            if (ru == rv) continue;
            if (!adj[ru * nr + rv]) { adj[ru * nr + rv] = 1; indeg[rv]++; }
        }
        int *ord = (int *)xmalloc((size_t)(nr ? nr : 1) * sizeof(int));
        int qh = 0, qt = 0;
        for (int k = 0; k < nr; k++) if (indeg[k] == 0) ord[qt++] = k;
        while (qh < qt) { int u = ord[qh++]; for (int k = 0; k < nr; k++) if (adj[u * nr + k]) { if (--indeg[k] == 0) ord[qt++] = k; } }
        if (qt == nr) {
            int *tmp = (int *)xmalloc((size_t)(nr ? nr : 1) * sizeof(int));
            for (int k = 0; k < nr; k++) tmp[k] = H->roots[ord[k]];
            memcpy(H->roots, tmp, (size_t)nr * sizeof(int));
            free(tmp);
        }
        free(adj); free(indeg); free(ord); free(leaf_root);
    }
    hier_layout(H);   /* 根次序确定后最终布局 */

    /* ---- 7) 相邻子节点代表点之间的物理间距 ---- */
    {
        H->gap = (len_t *)xrealloc(H->gap, (size_t)(nch > 0 ? nch : 1) * sizeof(len_t));
        for (int i = 0; i < nch; i++) H->gap[i] = len_zero();
        for (int i = 0; i < H->nnode; i++) {
            if (H->is_leaf[i]) continue;
            for (int k = 1; k < H->child_cnt[i]; k++) {
                int c0 = H->children[H->child_start[i] + k - 1];
                int c1 = H->children[H->child_start[i] + k];
                int s0 = H->leaf_state[H->first_leaf[c0]];
                int s1 = H->leaf_state[H->first_leaf[c1]];
                H->gap[H->child_start[i] + k] = tree_distance(&FT, s0, s1);
            }
        }
    }

    free(F.head); free(F.nxt); free(F.to); free(F.w);
    tree_free(&FT);
    for (int i = 0; i < nm; i++) ivec_free(&mchild[i]);
    free(mchild); free(mh); free(msum);
    free(item_node); free(mnode); free(comp_item);
    dsu_free(&d);
    ivec_free(&roots);
    (void)nroot_out; (void)roots_out; (void)aux_state;
    return 0;
}
/* 一致性检查：子节点区间必须不重不漏地铺满父区间 */
void hier_validate(const Hierarchy *H)
{
    int bad = 0;
    for (int i = 0; i < H->nnode; i++) {
        if (H->is_leaf[i]) continue;
        int expect = H->first_leaf[i];
        for (int k = 0; k < H->child_cnt[i]; k++) {
            int c = H->children[H->child_start[i] + k];
            if (c < 0 || c >= H->nnode) { bad++; continue; }
            if (H->first_leaf[c] != expect) {
                if (bad < 4) fprintf(stderr, "[val] 节点 %d 第 %d 子 %d 起点 %d 期望 %d\n", i, k, c, H->first_leaf[c], expect);
                bad++;
            }
            expect = H->last_leaf[c] + 1;
        }
        if (expect != H->last_leaf[i] + 1) {
            if (bad < 4) fprintf(stderr, "[val] 节点 %d 子区间未覆盖到 %d\n", i, H->last_leaf[i]);
            bad++;
        }
    }
    if (bad) fprintf(stderr, "[val] 共 %d 处区间问题\n", bad);
    else fprintf(stderr, "[val] 层级一致性检查通过\n");
}