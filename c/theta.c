/* theta.c -- 阈值分量（块与 theta）与连接记录的生成（§5.3、§8.1、引理 7.3、8.3）
 *
 * 块 = 状态图里成对反向出现的近邻弧的连通分量。
 * 用"最小最大"多源 Dijkstra 在去掉宏观弧后的状态图上精确求四个量：
 *   I_0(z)：能到达 z 的宏观边的最小权重（忽略径向权重）
 *   O_0(z)：从 z 出发能穿越的宏观边的最小权重
 *   I_*(z)、O_*(z)：同一批路径、但取整条访问路径的瓶颈
 * B_0 = max(I_0,O_0)，theta = max(I_*,O_*)；引理 8.3 给出 B_0 <= theta <= 4K^2 h B_0。
 *
 * 连接记录（引理 7.3 / 8.3）：块 MST 边、每个状态按 theta 连到其归属枢纽的
 * 辅助顶点、枢纽走廊按 sigma 相连。对这些记录按倍增网格的类同时合并，
 * 就恰好得到各阈值下的强连通分量。
 */
#include "sssp.h"

typedef struct { int n; int *head, *nxt, *to; len_t *w; int narc; } AL;
static void al_init(AL *a, int n, int narc)
{
    a->n = n;
    a->head = (int *)xmalloc((size_t)(n ? n : 1) * sizeof(int));
    for (int i = 0; i < n; i++) a->head[i] = -1;
    a->nxt = (int *)xmalloc((size_t)(narc ? narc : 1) * sizeof(int));
    a->to = (int *)xmalloc((size_t)(narc ? narc : 1) * sizeof(int));
    a->w = (len_t *)xmalloc((size_t)(narc ? narc : 1) * sizeof(len_t));
    a->narc = 0;
}
static void al_add(AL *a, int u, int v, len_t w)
{
    int i = a->narc++;
    a->to[i] = v; a->w[i] = w; a->nxt[i] = a->head[u]; a->head[u] = i;
}
static void al_free(AL *a) { free(a->head); free(a->nxt); free(a->to); free(a->w); memset(a, 0, sizeof(*a)); }

typedef struct { len_t *key; int *it; int *pos; int n, cap; } MH;
static void mh_init(MH *h, int n)
{
    h->n = 0; h->cap = n + 4;
    h->key = (len_t *)xmalloc((size_t)h->cap * sizeof(len_t));
    h->it = (int *)xmalloc((size_t)h->cap * sizeof(int));
    h->pos = (int *)xmalloc((size_t)n * sizeof(int));
    for (int i = 0; i < n; i++) h->pos[i] = -1;
}
static void mh_free(MH *h) { free(h->key); free(h->it); free(h->pos); }
static void mh_push(MH *h, int v, len_t k)
{
    int i;
    if (h->pos[v] >= 0) {
        i = h->pos[v];
        if (!len_lt(k, h->key[i])) return;
        h->key[i] = k;
    } else {
        i = h->n++;
        h->it[i] = v; h->key[i] = k; h->pos[v] = i;
    }
    while (i > 0) {
        int p = (i - 1) / 2;
        if (len_lt(h->key[i], h->key[p])) {
            len_t t = h->key[i]; h->key[i] = h->key[p]; h->key[p] = t;
            int ti = h->it[i]; h->it[i] = h->it[p]; h->it[p] = ti;
            h->pos[h->it[i]] = i; h->pos[h->it[p]] = p; i = p;
        } else break;
    }
}
static int mh_pop(MH *h, len_t *out)
{
    if (h->n == 0) return -1;
    int v = h->it[0];
    *out = h->key[0];
    h->pos[v] = -1;
    h->n--;
    if (h->n > 0) {
        h->it[0] = h->it[h->n]; h->key[0] = h->key[h->n]; h->pos[h->it[0]] = 0;
        int i = 0;
        for (;;) {
            int l = 2 * i + 1, r = l + 1, m = i;
            if (l < h->n && len_lt(h->key[l], h->key[m])) m = l;
            if (r < h->n && len_lt(h->key[r], h->key[m])) m = r;
            if (m == i) break;
            len_t t = h->key[i]; h->key[i] = h->key[m]; h->key[m] = t;
            int ti = h->it[i]; h->it[i] = h->it[m]; h->it[m] = ti;
            h->pos[h->it[i]] = i; h->pos[h->it[m]] = m; i = m;
        }
    }
    return v;
}

static void bdij(int n, const AL *A, const int *src, const len_t *sc, int nsrc, int use_w, len_t *out)
{
    for (int i = 0; i < n; i++) out[i] = len_inf();
    MH h; mh_init(&h, n);
    for (int s = 0; s < nsrc; s++)
        if (len_lt(sc[s], out[src[s]])) { out[src[s]] = sc[s]; mh_push(&h, src[s], sc[s]); }
    char *done = (char *)xcalloc((size_t)(n ? n : 1), 1);
    for (;;) {
        len_t du; int u = mh_pop(&h, &du);
        if (u < 0) break;
        if (done[u]) continue;
        done[u] = 1;
        if (len_gt(du, out[u])) continue;
        for (int e = A->head[u]; e != -1; e = A->nxt[e]) {
            int v = A->to[e];
            len_t nd = use_w ? len_max(du, A->w[e]) : du;
            if (len_lt(nd, out[v])) { out[v] = nd; mh_push(&h, v, nd); }
        }
    }
    free(done); mh_free(&h);
}

int outer_build_hierarchy(Hierarchy *H, const SGraph *S, const Owner *O,
                          int nstates, const int *aux_state, int naux)
{
    int total = nstates + naux;

    /* ---- 块：成对反向弧 ---- */
    typedef struct { int u, v; len_t w; } AR;
    AR *ar = (AR *)xmalloc((size_t)(S->narc ? S->narc : 1) * sizeof(AR));
    for (int i = 0; i < S->narc; i++) { ar[i].u = S->arc_head[i]; ar[i].v = S->arc_to[i]; ar[i].w = S->arc_len[i]; }
    /* 以 (u,v) 排序，再按权重配对 */
    for (int i = 1; i < S->narc; i++) {
        AR k = ar[i]; int j = i - 1;
        while (j >= 0 && (ar[j].u > k.u || (ar[j].u == k.u && ar[j].v > k.v))) { ar[j + 1] = ar[j]; j--; }
        ar[j + 1] = k;
    }
    DSU d; dsu_init(&d, nstates);
    {
        int lo = 0;
        while (lo < S->narc) {
            int hi = lo;
            while (hi < S->narc && ar[hi].u == ar[lo].u && ar[hi].v == ar[lo].v) hi++;
            int a = ar[lo].u, b = ar[lo].v;
            /* 反向 (b,a) 是否存在 */
            int l = 0, r = S->narc - 1, found = 0;
            while (l <= r) {
                int mid = (l + r) / 2;
                if (ar[mid].u == b && ar[mid].v == a) { found = 1; break; }
                if (ar[mid].u < b || (ar[mid].u == b && ar[mid].v < a)) l = mid + 1; else r = mid - 1;
            }
            if (found && a != b) dsu_unite(&d, a, b);
            lo = hi;
        }
        for (int i = 0; i < S->narc; i++) if (ar[i].u == ar[i].v) { /* 自环不影响块 */ }
    }
    int *blk = (int *)xmalloc((size_t)(nstates ? nstates : 1) * sizeof(int));
    {
        int *map = (int *)xmalloc((size_t)(nstates ? nstates : 1) * sizeof(int));
        for (int i = 0; i < nstates; i++) map[i] = -1;
        int nb = 0;
        for (int i = 0; i < nstates; i++) {
            int r = dsu_find(&d, i);
            if (map[r] < 0) map[r] = nb++;
            blk[i] = map[r];
        }
        free(map);
    }
    dsu_free(&d);

    /* ---- 宏观弧源 ---- */
    int nmacro = 0;
    for (int i = 0; i < S->narc; i++) if (S->arc_macro[i]) nmacro++;
    int *msrc = (int *)xmalloc((size_t)(nmacro ? nmacro : 1) * sizeof(int));
    int *mdst = (int *)xmalloc((size_t)(nmacro ? nmacro : 1) * sizeof(int));
    len_t *mw = (len_t *)xmalloc((size_t)(nmacro ? nmacro : 1) * sizeof(len_t));
    {
        int k = 0;
        for (int i = 0; i < S->narc; i++) if (S->arc_macro[i]) {
            msrc[k] = S->arc_head[i]; mdst[k] = S->arc_to[i]; mw[k] = S->arc_len[i]; k++;
        }
    }

    /* ---- 四个阈值量 ---- */
    len_t *I0 = (len_t *)xmalloc((size_t)nstates * sizeof(len_t));
    len_t *O0 = (len_t *)xmalloc((size_t)nstates * sizeof(len_t));
    len_t *Is = (len_t *)xmalloc((size_t)nstates * sizeof(len_t));
    len_t *Os = (len_t *)xmalloc((size_t)nstates * sizeof(len_t));
    AL fwd, rev;
    al_init(&fwd, nstates, S->narc);
    al_init(&rev, nstates, S->narc);
    for (int i = 0; i < S->narc; i++) {
        al_add(&fwd, S->arc_head[i], S->arc_to[i], S->arc_len[i]);
        al_add(&rev, S->arc_to[i], S->arc_head[i], S->arc_len[i]);
    }
    /* 把宏观弧从"无宏观"图里去掉：这里直接用完整图但把宏观弧权重当作 0 处理？
       引理要求宏-free 路径；宏观弧本身只作为种子/终点，不参与内部转移。
       实现上：把宏观弧的弧长置为 inf 会让它不可用，而种子已经单独给出。 */
    {
        AL fwd2, rev2;
        al_init(&fwd2, nstates, S->narc - nmacro);
        al_init(&rev2, nstates, S->narc - nmacro);
        for (int i = 0; i < S->narc; i++) {
            if (S->arc_macro[i]) continue;
            al_add(&fwd2, S->arc_head[i], S->arc_to[i], S->arc_len[i]);
            al_add(&rev2, S->arc_to[i], S->arc_head[i], S->arc_len[i]);
        }
        bdij(nstates, &fwd2, mdst, mw, nmacro, 0, I0);
        bdij(nstates, &rev2, msrc, mw, nmacro, 0, O0);
        bdij(nstates, &fwd2, mdst, mw, nmacro, 1, Is);
        bdij(nstates, &rev2, msrc, mw, nmacro, 1, Os);
        al_free(&fwd2); al_free(&rev2);
    }
    al_free(&fwd); al_free(&rev);

    /* ---- 连接记录 ---- */
    ivec ju, jv, jsu, jsv;
    lvec jw;
    ivec_init(&ju); ivec_init(&jv); ivec_init(&jsu); ivec_init(&jsv); lvec_init(&jw);

    /* (a) 块 MST 边：块内经 Kruskal 取最小生成树的边 */
    {
        typedef struct { len_t w; int u, v; } E2;
        int cap = S->narc;
        E2 *es = (E2 *)xmalloc((size_t)(cap ? cap : 1) * sizeof(E2));
        int ne = 0;
        for (int i = 0; i < S->narc; i++) {
            int u = S->arc_head[i], v = S->arc_to[i];
            if (u >= v) continue;                    /* 每对只取一次 */
            if (blk[u] != blk[v]) continue;
            es[ne].u = u; es[ne].v = v; es[ne].w = S->arc_len[i]; ne++;
        }
        for (int i = 1; i < ne; i++) {
            E2 k = es[i]; int j = i - 1;
            while (j >= 0 && len_gt(es[j].w, k.w)) { es[j + 1] = es[j]; j--; }
            es[j + 1] = k;
        }
        DSU bd; dsu_init(&bd, nstates);
        for (int i = 0; i < ne; i++) {
            if (dsu_unite(&bd, es[i].u, es[i].v)) {
                ivec_push(&ju, es[i].u); ivec_push(&jv, es[i].v);
                lvec_push(&jw, es[i].w);
                ivec_push(&jsu, es[i].u); ivec_push(&jsv, es[i].v);
            }
        }
        dsu_free(&bd);
        free(es);
    }
    /* (b) 状态 -> 归属枢纽的辅助顶点，权重 theta */
    for (int z = 0; z < nstates; z++) {
        len_t th = len_max(len_max(Is[z], Os[z]), len_zero());
        if (!len_finite(th)) continue;
        int v = S->proj[z];
        int p = O->owner[v];
        int aux = nstates + p;
        if (aux >= total) continue;
        ivec_push(&ju, z); ivec_push(&jv, aux);
        lvec_push(&jw, th);
        ivec_push(&jsu, z);
        ivec_push(&jsv, aux_state && aux_state[p] >= 0 ? aux_state[p] : z);
    }
    /* (c) 枢纽走廊 */
    for (int i = 0; i < O->npiv_edge; i++) {
        int a = nstates + O->piv_edge_u[i], b = nstates + O->piv_edge_v[i];
        ivec_push(&ju, a); ivec_push(&jv, b);
        lvec_push(&jw, O->piv_edge_w[i]);
        ivec_push(&jsu, aux_state ? aux_state[O->piv_edge_u[i]] : a);
        ivec_push(&jsv, aux_state ? aux_state[O->piv_edge_v[i]] : b);
    }

    int njoin = ju.n;
    JoinRec *joins = (JoinRec *)xmalloc((size_t)(njoin ? njoin : 1) * sizeof(JoinRec));
    for (int i = 0; i < njoin; i++) {
        joins[i].u = ju.a[i]; joins[i].v = jv.a[i]; joins[i].w = jw.a[i];
        joins[i].su = jsu.a[i]; joins[i].sv = jsv.a[i];
    }
    hier_build_from_joins(H, S, nstates, naux, aux_state, njoin, joins, 0, NULL);
    H->nstate = nstates;

    free(joins); ivec_free(&ju); ivec_free(&jv); ivec_free(&jsu); ivec_free(&jsv); lvec_free(&jw);
    free(msrc); free(mdst); free(mw); free(blk); free(ar);
    free(I0); free(O0); free(Is); free(Os);
    return 0;
}
