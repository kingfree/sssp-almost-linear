/* route.c -- 第 6 节：由随机样本得到短的枢纽路由证书（算法 2、3）
 *
 * 算法 3 Estimate（样本距离估计器）：
 *   1. 按权重排序样本边，用 Kruskal 扫描建"合并森林"：每次成功合并生成一个
 *      二元父节点，记录合并权重与权重带；
 *   2. 需求 (x,y) 的合并点 = 森林中两叶子的最近公共祖先；
 *   3. 对权重带锚点 a：收缩 < a/T^2 的边，保留 [a/T^2, T^2 a] 的边得商图 Q_a；
 *   4. 尺度 Delta = a*2^s 上做随机中心分解：随机置换 + 半径 Delta(1+j/H) 的
 *      剪枝 Dijkstra，顶点取置换序中第一个覆盖它的中心；
 *   5. 由粗到细串成划分树（每个划分节点存值 4*Delta + a），需求取两端叶子
 *      的最近公共划分 -> 上界 Z。
 *
 * 算法 2 NameRoutes（证书分配）：反复抽样估计，按 Z 分组并用 P*a_* >= Z 的
 * 阈值找出最小可用输入权重；至少一半边达标则提交并从 U 删除；残留 <= 32t 条
 * 各以自己的权重为证书；最后排序去重得到词典。
 */
#include "sssp.h"

typedef struct { int u, v; len_t w; int eid; } SampleEdge;

static int cmp_se(const void *pa, const void *pb)
{
    const SampleEdge *a = (const SampleEdge *)pa, *b = (const SampleEdge *)pb;
    int c = len_cmp(a->w, b->w);
    if (c) return c;
    return a->eid - b->eid;
}

static int log2_ceil_int(long long x)
{
    int k = 0; long long p = 1;
    while (p < x) { p <<= 1; k++; }
    return k;
}
static double pow2_ceil_d(double x)
{
    double p = 1.0;
    while (p < x) p *= 2.0;
    return p;
}

/* ------------------------------------------------------------------ */
/* 算法 3：样本中的距离估计                                            */
/* ------------------------------------------------------------------ */

typedef struct {
    int n;
    int *head, *nxt, *to, *_pad;
    len_t *w;
    int narc;
} Small;

static void small_init(Small *s, int n)
{
    s->n = n;
    s->head = (int *)xmalloc((size_t)(n ? n : 1) * sizeof(int));
    for (int i = 0; i < n; i++) s->head[i] = -1;
    s->nxt = NULL; s->to = NULL; s->w = NULL; s->narc = 0; s->_pad = NULL;
}
static void small_add(Small *s, int u, int v, len_t w)
{
    s->to = (int *)xrealloc(s->to, (size_t)(s->narc + 1) * sizeof(int));
    s->nxt = (int *)xrealloc(s->nxt, (size_t)(s->narc + 1) * sizeof(int));
    s->w = (len_t *)xrealloc(s->w, (size_t)(s->narc + 1) * sizeof(len_t));
    s->to[s->narc] = v; s->w[s->narc] = w; s->nxt[s->narc] = s->head[u]; s->head[u] = s->narc;
    s->narc++;
}
static void small_free(Small *s)
{
    free(s->head); free(s->nxt); free(s->to); free(s->w);
    memset(s, 0, sizeof(*s));
}

/* 小顶堆（商图上用） */
typedef struct { int *it; len_t *key; int *pos; int n, cap; } MiniHeap;
static void mh_init(MiniHeap *h, int n)
{
    h->cap = n + 4; h->n = 0;
    h->it = (int *)xmalloc((size_t)h->cap * sizeof(int));
    h->key = (len_t *)xmalloc((size_t)h->cap * sizeof(len_t));
    h->pos = (int *)xmalloc((size_t)n * sizeof(int));
    for (int i = 0; i < n; i++) h->pos[i] = -1;
}
static void mh_free(MiniHeap *h) { free(h->it); free(h->key); free(h->pos); }
static void mh_push(MiniHeap *h, int v, len_t k)
{
    if (h->pos[v] >= 0) {
        int i = h->pos[v];
        if (!len_lt(k, h->key[i])) return;
        h->key[i] = k;
        while (i > 0) {
            int p = (i - 1) / 2;
            if (len_lt(h->key[i], h->key[p])) {
                len_t t = h->key[i]; h->key[i] = h->key[p]; h->key[p] = t;
                int ti = h->it[i]; h->it[i] = h->it[p]; h->it[p] = ti;
                h->pos[h->it[i]] = i; h->pos[h->it[p]] = p;
                i = p;
            } else break;
        }
        return;
    }
    int i = h->n++;
    h->it[i] = v; h->key[i] = k; h->pos[v] = i;
    while (i > 0) {
        int p = (i - 1) / 2;
        if (len_lt(h->key[i], h->key[p])) {
            len_t t = h->key[i]; h->key[i] = h->key[p]; h->key[p] = t;
            int ti = h->it[i]; h->it[i] = h->it[p]; h->it[p] = ti;
            h->pos[h->it[i]] = i; h->pos[h->it[p]] = p;
            i = p;
        } else break;
    }
}
static int mh_pop(MiniHeap *h, len_t *out)
{
    if (h->n == 0) return -1;
    int v = h->it[0];
    if (out) *out = h->key[0];
    h->pos[v] = -1;
    h->n--;
    if (h->n > 0) {
        h->it[0] = h->it[h->n]; h->key[0] = h->key[h->n];
        h->pos[h->it[0]] = 0;
        int i = 0;
        for (;;) {
            int l = 2 * i + 1, r = l + 1, m = i;
            if (l < h->n && len_lt(h->key[l], h->key[m])) m = l;
            if (r < h->n && len_lt(h->key[r], h->key[m])) m = r;
            if (m == i) break;
            len_t t = h->key[i]; h->key[i] = h->key[m]; h->key[m] = t;
            int ti = h->it[i]; h->it[i] = h->it[m]; h->it[m] = ti;
            h->pos[h->it[i]] = i; h->pos[h->it[m]] = m;
            i = m;
        }
    }
    return v;
}

void estimate_sample(int t, const int *su, const int *sv, const len_t *sw, int ns,
                     int h, const int *qu, const int *qv,
                     int *out_index, int *out_disc, int *out_ntab, len_t **out_table)
{
    for (int i = 0; i < h; i++) { out_index[i] = -1; out_disc[i] = 0; }
    *out_ntab = 0; *out_table = NULL;
    if (t < 2 || ns == 0) { for (int i = 0; i < h; i++) out_disc[i] = 1; return; }

    SampleEdge *es = (SampleEdge *)xmalloc((size_t)ns * sizeof(SampleEdge));
    for (int i = 0; i < ns; i++) { es[i].u = su[i]; es[i].v = sv[i]; es[i].w = sw[i]; es[i].eid = i; }
    qsort(es, (size_t)ns, sizeof(SampleEdge), cmp_se);

    long long Tq = 4LL * t;
    int logT = log2_ceil_int(Tq);
    double T = pow2_ceil_d((double)Tq);

    /* ---- 权重带 ---- */
    int *band_of_edge = (int *)xmalloc((size_t)ns * sizeof(int));
    lvec anchors; lvec_init(&anchors);
    {
        int i = 0;
        int nb = 0;
        while (i < ns) {
            len_t a = es[i].w;
            len_t lim = len_mul_pow2(a, logT);
            lvec_push(&anchors, a);
            while (i < ns && len_lt(es[i].w, lim)) { band_of_edge[i] = nb; i++; }
            nb++;
        }
    }
    int nband = anchors.n;

    /* ---- Kruskal 合并森林 ---- */
    int maxnode = 2 * t + 2;
    int *mleft = (int *)xmalloc((size_t)maxnode * sizeof(int));
    int *mright = (int *)xmalloc((size_t)maxnode * sizeof(int));
    int *mband = (int *)xmalloc((size_t)maxnode * sizeof(int));
    int *mpar = (int *)xmalloc((size_t)maxnode * sizeof(int));
    int *mdepth = (int *)xcalloc((size_t)maxnode, sizeof(int));
    for (int i = 0; i < maxnode; i++) { mleft[i] = mright[i] = -1; mpar[i] = -1; mband[i] = -1; }
    int *compnode = (int *)xmalloc((size_t)t * sizeof(int));
    for (int i = 0; i < t; i++) compnode[i] = i;
    DSU d; dsu_init(&d, t);
    int nextnode = t;
    for (int i = 0; i < ns; i++) {
        int ra = dsu_find(&d, es[i].u), rb = dsu_find(&d, es[i].v);
        if (ra == rb) continue;
        if (nextnode >= maxnode) break;
        int k = nextnode++;
        mleft[k] = compnode[ra];
        mright[k] = compnode[rb];
        mband[k] = band_of_edge[i];
        mpar[mleft[k]] = k;
        mpar[mright[k]] = k;
        int nr = dsu_find(&d, dsu_unite(&d, ra, rb));
        compnode[nr] = k;
    }
    {
        /* 自顶向下推深度（父节点编号大于子节点） */
        for (int k = nextnode - 1; k >= t; k--) {
            if (mleft[k] >= 0) mdepth[mleft[k]] = mdepth[k] + 1;
            if (mright[k] >= 0) mdepth[mright[k]] = mdepth[k] + 1;
        }
    }
    int *root_of = (int *)xmalloc((size_t)t * sizeof(int));
    for (int i = 0; i < t; i++) root_of[i] = dsu_find(&d, i);

    int *dem_band = (int *)xmalloc((size_t)(h ? h : 1) * sizeof(int));
    for (int i = 0; i < h; i++) {
        dem_band[i] = -1;
        if (root_of[qu[i]] != root_of[qv[i]]) { out_disc[i] = 1; continue; }
        int a = qu[i], b = qv[i];
        while (mdepth[a] > mdepth[b]) a = mpar[a];
        while (mdepth[b] > mdepth[a]) b = mpar[b];
        while (a != b) { a = mpar[a]; b = mpar[b]; }
        dem_band[i] = mband[a];
    }

    /* ---- 划分树 ---- */
    int maxs = log2_ceil_int((long long)((double)t * T * T) + 1);
    int capn = (t + 1) * (maxs + 3) + 4;
    int *pn_par = (int *)xmalloc((size_t)capn * sizeof(int));
    int *pn_dep = (int *)xmalloc((size_t)capn * sizeof(int));
    len_t *pn_val = (len_t *)xmalloc((size_t)capn * sizeof(len_t));
    int npn = 0;
    for (int i = 0; i < capn; i++) { pn_par[i] = -1; pn_dep[i] = 0; pn_val[i] = len_zero(); }

    int *qpair_u = (int *)xmalloc((size_t)(h ? h : 1) * sizeof(int));
    int *qpair_v = (int *)xmalloc((size_t)(h ? h : 1) * sizeof(int));

    /* 按 band 分组的需求顺序 */
    int *ord = (int *)xmalloc((size_t)(h ? h : 1) * sizeof(int));
    {
        int *cnt = (int *)xcalloc((size_t)nband + 2, sizeof(int));
        for (int i = 0; i < h; i++) if (dem_band[i] >= 0) cnt[dem_band[i] + 1]++;
        for (int i = 0; i < nband + 1; i++) cnt[i + 1] += cnt[i];
        int *pos = (int *)xmalloc((size_t)(nband + 2) * sizeof(int));
        for (int i = 0; i < nband + 2; i++) pos[i] = cnt[i];
        for (int i = 0; i < h; i++)
            if (dem_band[i] >= 0) ord[pos[dem_band[i]]++] = i;
        free(cnt); free(pos);
    }

    int p = 0;
    while (p < h) {
        int band = dem_band[ord[p]];
        int q = p;
        while (q < h && dem_band[ord[q]] == band) q++;
        len_t a = anchors.a[band];
        len_t lo = len_mul_pow2(a, -2 * logT);
        len_t hi = len_mul_pow2(a, 2 * logT);

        /* 收缩 + 保留边 */
        DSU cd; dsu_init(&cd, t);
        int i = 0;
        while (i < ns && len_lt(es[i].w, lo)) { dsu_unite(&cd, es[i].u, es[i].v); i++; }
        int *ku = (int *)xmalloc((size_t)(ns + 1) * sizeof(int));
        int *kv = (int *)xmalloc((size_t)(ns + 1) * sizeof(int));
        len_t *kw = (len_t *)xmalloc((size_t)(ns + 1) * sizeof(len_t));
        int nk = 0;
        int *map = (int *)xmalloc((size_t)t * sizeof(int));
        for (int j = 0; j < t; j++) map[j] = -1;
        int nq = 0;
        for (int j = i; j < ns && len_le(es[j].w, hi); j++) {
            int cu = dsu_find(&cd, es[j].u), cv = dsu_find(&cd, es[j].v);
            if (cu == cv) continue;
            ku[nk] = cu; kv[nk] = cv; kw[nk] = es[j].w; nk++;
            if (map[cu] < 0) map[cu] = nq++;
            if (map[cv] < 0) map[cv] = nq++;
        }
        for (int j = p; j < q; j++) {
            int idx = ord[j];
            int cu = dsu_find(&cd, qu[idx]), cv = dsu_find(&cd, qv[idx]);
            if (cu == cv || map[cu] < 0 || map[cv] < 0) {
                /* 理论上不会发生（引理 6.4）；防御性地记为不连通 */
                out_disc[idx] = 1;
                qpair_u[j - p] = -1; qpair_v[j - p] = -1;
            } else {
                qpair_u[j - p] = map[cu];
                qpair_v[j - p] = map[cv];
            }
        }
        dsu_free(&cd);

        if (nq > 0) {
            Small sm; small_init(&sm, nq);
            for (int j = 0; j < nk; j++) {
                small_add(&sm, map[ku[j]], map[kv[j]], kw[j]);
                small_add(&sm, map[kv[j]], map[ku[j]], kw[j]);
            }
            int *perm = (int *)xmalloc((size_t)nq * sizeof(int));
            len_t *best = (len_t *)xmalloc((size_t)nq * sizeof(len_t));
            int *label = (int *)xmalloc((size_t)nq * sizeof(int));
            len_t *dist = (len_t *)xmalloc((size_t)nq * sizeof(len_t));
            char *done = (char *)xmalloc((size_t)nq);
            int *prev_part = (int *)xmalloc((size_t)nq * sizeof(int));
            int *cur_part = (int *)xmalloc((size_t)nq * sizeof(int));
            int *sorti = (int *)xmalloc((size_t)nq * sizeof(int));
            int *tmp = (int *)xmalloc((size_t)nq * sizeof(int));
            int *tmp2 = (int *)xmalloc((size_t)nq * sizeof(int));
            int *leafn = (int *)xmalloc((size_t)nq * sizeof(int));
            MiniHeap mh; mh_init(&mh, nq);

            /* 根：整个商图分量，值 4a + a */
            int root = npn++;
            pn_par[root] = -1; pn_dep[root] = 0;
            pn_val[root] = len_add(len_mul_pow2(a, 2), a);
            for (int j = 0; j < nq; j++) prev_part[j] = root;

            for (int s = 0; s <= maxs; s++) {
                len_t Delta = len_mul_pow2(a, s);
                for (int j = 0; j < nq; j++) perm[j] = j;
                for (int j = nq - 1; j > 0; j--) { int k = rng_int(j + 1); int tp = perm[j]; perm[j] = perm[k]; perm[k] = tp; }
                len_t R = len_make(Delta.v * (1.0 + rng_double()));
                for (int j = 0; j < nq; j++) { best[j] = len_inf(); label[j] = -1; }
                for (int ci = 0; ci < nq; ci++) {
                    int c = perm[ci];
                    for (int j = 0; j < nq; j++) { dist[j] = len_inf(); done[j] = 0; }
                    for (int j = 0; j < nq; j++) mh.pos[j] = -1;
                    mh.n = 0;
                    dist[c] = len_zero();
                    mh_push(&mh, c, len_zero());
                    for (;;) {
                        len_t du;
                        int u = mh_pop(&mh, &du);
                        if (u < 0) break;
                        if (done[u]) continue;
                        done[u] = 1;
                        if (len_gt(du, R)) break;
                        if (len_le(best[u], du)) continue;   /* 剪枝：已有更早中心的更短覆盖 */
                        best[u] = du;
                        if (label[u] < 0) label[u] = ci;
                        for (int e = sm.head[u]; e != -1; e = sm.nxt[e]) {
                            int v = sm.to[e];
                            len_t nd = len_add(du, sm.w[e]);
                            if (len_le(nd, R) && len_lt(nd, dist[v])) { dist[v] = nd; mh_push(&mh, v, nd); }
                        }
                    }
                }
                /* 本尺度划分 = (上一尺度划分, 中心标号) */
                for (int j = 0; j < nq; j++) sorti[j] = j;
                {
                    int *cnt = (int *)xcalloc((size_t)nq + 1, sizeof(int));
                    for (int j = 0; j < nq; j++) cnt[(label[j] < 0 ? 0 : label[j]) + 1]++;
                    for (int j = 0; j < nq; j++) cnt[j + 1] += cnt[j];
                    for (int j = 0; j < nq; j++) tmp[cnt[label[j] < 0 ? 0 : label[j]]++] = j;
                    int *cnt2 = (int *)xcalloc((size_t)nq + 1, sizeof(int));
                    for (int j = 0; j < nq; j++) cnt2[(prev_part[tmp[j]] < 0 ? 0 : prev_part[tmp[j]]) + 1]++;
                    for (int j = 0; j < nq; j++) cnt2[j + 1] += cnt2[j];
                    for (int j = 0; j < nq; j++) tmp2[cnt2[prev_part[tmp[j]] < 0 ? 0 : prev_part[tmp[j]]]++] = tmp[j];
                    for (int j = 0; j < nq; j++) sorti[j] = tmp2[j];
                    free(cnt); free(cnt2);
                }
                int last = -1;
                for (int j = 0; j < nq; j++) {
                    int v = sorti[j];
                    int same = (j > 0) && (prev_part[v] == prev_part[sorti[j - 1]]) && (label[v] == label[sorti[j - 1]]);
                    if (!same) {
                        last = npn++;
                        if (last >= capn) { fprintf(stderr, "estimate: 划分树节点溢出\n"); exit(1); }
                        pn_par[last] = prev_part[v];
                        pn_dep[last] = pn_dep[pn_par[last]] + 1;
                        pn_val[last] = len_add(len_mul_pow2(Delta, 2), a);
                    }
                    cur_part[v] = last;
                }
                for (int j = 0; j < nq; j++) prev_part[j] = cur_part[j];
            }
            /* 每个顶点一片叶子 */
            for (int j = 0; j < nq; j++) {
                int nd = npn++;
                if (nd >= capn) { fprintf(stderr, "estimate: 叶子溢出\n"); exit(1); }
                pn_par[nd] = cur_part[j];
                pn_dep[nd] = pn_dep[pn_par[nd]] + 1;
                pn_val[nd] = len_zero();
                leafn[j] = nd;
            }
            /* 需求：两端叶子的最近公共祖先 */
            for (int j = p; j < q; j++) {
                int idx = ord[j];
                if (qpair_u[j - p] < 0) continue;
                int x = leafn[qpair_u[j - p]], y = leafn[qpair_v[j - p]];
                while (pn_dep[x] > pn_dep[y]) x = pn_par[x];
                while (pn_dep[y] > pn_dep[x]) y = pn_par[y];
                while (x != y) { x = pn_par[x]; y = pn_par[y]; }
                out_index[idx] = x;
            }
            mh_free(&mh);
            free(perm); free(best); free(label); free(dist); free(done);
            free(prev_part); free(cur_part); free(sorti); free(tmp); free(tmp2); free(leafn);
            small_free(&sm);
        }
        free(ku); free(kv); free(kw); free(map);
        p = q;
    }

    *out_ntab = npn;
    *out_table = pn_val;
    free(pn_par); free(pn_dep);
    free(es); free(band_of_edge); free(anchors.a);
    free(mleft); free(mright); free(mband); free(mpar); free(mdepth);
    free(compnode); free(root_of); free(dem_band); free(ord);
    free(qpair_u); free(qpair_v);
    dsu_free(&d);
}

/* ------------------------------------------------------------------ */
/* 算法 2：证书分配                                                    */
/* ------------------------------------------------------------------ */

typedef struct { int u, v; len_t w; int rid; } ReqKey;
static int cmp_reqkey(const void *pa, const void *pb)
{
    const ReqKey *a = (const ReqKey *)pa, *b = (const ReqKey *)pb;
    if (a->u != b->u) return a->u - b->u;
    if (a->v != b->v) return a->v - b->v;
    return a->rid - b->rid;
}

void routes_build(Routes *R, int npiv, int nreq, const int *ru, const int *rv,
                  const len_t *rw, uint64_t seed)
{
    memset(R, 0, sizeof(*R));
    rng_seed(seed);
    int t = npiv;
    R->npiv = t;
    R->nreq = nreq;
    R->ru = (int *)xmalloc((size_t)(nreq ? nreq : 1) * sizeof(int));
    R->rv = (int *)xmalloc((size_t)(nreq ? nreq : 1) * sizeof(int));
    R->rw = (len_t *)xmalloc((size_t)(nreq ? nreq : 1) * sizeof(len_t));
    R->name = (len_t *)xmalloc((size_t)(nreq ? nreq : 1) * sizeof(len_t));
    R->dict_of_req = (int *)xmalloc((size_t)(nreq ? nreq : 1) * sizeof(int));
    for (int i = 0; i < nreq; i++) {
        R->ru[i] = ru[i]; R->rv[i] = rv[i]; R->rw[i] = rw[i];
        R->name[i] = len_inf(); R->dict_of_req[i] = -1;
    }
    if (t < 2 || nreq == 0) { R->P = len_make(1.0); return; }

    int u = log2_ceil_int((long long)t + 2);
    int g = 4 * u + 3;
    double Pneed = 8.0 * 256.0 * (double)g * (double)u * (double)u;
    double Pv = 1.0;
    int Plog = 0;
    while (Pv < Pneed) { Pv *= 2.0; Plog++; }
    R->P = len_make(Pv);

    /* 每条枢纽对只留最轻的投影边作代表 */
    ReqKey *rk = (ReqKey *)xmalloc((size_t)nreq * sizeof(ReqKey));
    for (int i = 0; i < nreq; i++) {
        int a = ru[i], b = rv[i];
        if (a > b) { int tmp = a; a = b; b = tmp; }
        rk[i].u = a; rk[i].v = b; rk[i].w = rw[i]; rk[i].rid = i;
    }
    qsort(rk, (size_t)nreq, sizeof(ReqKey), cmp_reqkey);

    int nrep = 0;
    int *rep_u = (int *)xmalloc((size_t)nreq * sizeof(int));
    int *rep_v = (int *)xmalloc((size_t)nreq * sizeof(int));
    len_t *rep_w = (len_t *)xmalloc((size_t)nreq * sizeof(len_t));
    int *rep_lo = (int *)xmalloc((size_t)nreq * sizeof(int));
    int *rep_hi = (int *)xmalloc((size_t)nreq * sizeof(int));
    {
        int i = 0;
        while (i < nreq) {
            int j = i;
            len_t best = rk[i].w;
            while (j < nreq && rk[j].u == rk[i].u && rk[j].v == rk[i].v) {
                if (len_lt(rk[j].w, best)) best = rk[j].w;
                j++;
            }
            rep_u[nrep] = rk[i].u; rep_v[nrep] = rk[i].v; rep_w[nrep] = best;
            rep_lo[nrep] = i; rep_hi[nrep] = j;
            nrep++;
            i = j;
        }
    }

    /* 活的代表集合 U；提交时写入 name */
    int *live = (int *)xmalloc((size_t)nrep * sizeof(int));
    int nlive = nrep;
    for (int i = 0; i < nrep; i++) live[i] = i;
    len_t *lambda_of_rep = (len_t *)xmalloc((size_t)nrep * sizeof(len_t));
    for (int i = 0; i < nrep; i++) lambda_of_rep[i] = len_inf();

    /* 词典候选（提交过的 a_*） */
    lvec dictcand; lvec_init(&dictcand);

    int sample_cap = 16 * t;
    if (sample_cap < 2) sample_cap = 2;
    int *su = (int *)xmalloc((size_t)sample_cap * sizeof(int));
    int *sv = (int *)xmalloc((size_t)sample_cap * sizeof(int));
    len_t *sw = (len_t *)xmalloc((size_t)sample_cap * sizeof(len_t));
    int *du_arr = (int *)xmalloc((size_t)nrep * sizeof(int));
    int *dv_arr = (int *)xmalloc((size_t)nrep * sizeof(int));
    int *didx = (int *)xmalloc((size_t)nrep * sizeof(int));
    int *ddisc = (int *)xmalloc((size_t)nrep * sizeof(int));

    while (nlive > 32 * t) {
        int ns = 0;
        /* 从 live 中无放回抽 sample_cap 条 */
        {
            int *tmp = (int *)xmalloc((size_t)nlive * sizeof(int));
            memcpy(tmp, live, (size_t)nlive * sizeof(int));
            int want = sample_cap < nlive ? sample_cap : nlive;
            for (int i = 0; i < want; i++) {
                int j = i + rng_int(nlive - i);
                int t2 = tmp[i]; tmp[i] = tmp[j]; tmp[j] = t2;
                su[ns] = rep_u[tmp[i]]; sv[ns] = rep_v[tmp[i]]; sw[ns] = rep_w[tmp[i]];
                ns++;
            }
            free(tmp);
        }
        for (int i = 0; i < nlive; i++) { du_arr[i] = rep_u[live[i]]; dv_arr[i] = rep_v[live[i]]; }
        int ntab = 0;
        len_t *table = NULL;
        estimate_sample(t, su, sv, sw, ns, nlive, du_arr, dv_arr, didx, ddisc, &ntab, &table);

        /* 按 Z 值分组：把 (Z, 代表) 排序 */
        typedef struct { len_t z; int rep; int ok; } ZG;
        ZG *zg = (ZG *)xmalloc((size_t)nlive * sizeof(ZG));
        int nzg = 0;
        for (int i = 0; i < nlive; i++) {
            if (ddisc[i] || didx[i] < 0 || !table) { continue; }
            zg[nzg].z = table[didx[i]];
            zg[nzg].rep = live[i];
            zg[nzg].ok = 1;
            nzg++;
        }
        /* 组内按权重排序以做阈值二分 */
        int *grep = (int *)xmalloc((size_t)nzg * sizeof(int));
        int *gwo = (int *)xmalloc((size_t)nzg * sizeof(int));
        for (int i = 0; i < nzg; i++) { grep[i] = zg[i].rep; gwo[i] = i; }
        /* 用简单排序：按 (Z, 权重) */
        for (int i = 1; i < nzg; i++) {
            int key = gwo[i], j = i - 1;
            while (j >= 0) {
                int a2 = gwo[j];
                int c1 = len_cmp(zg[a2].z, zg[key].z);
                int c2 = c1 ? c1 : len_cmp(rep_w[zg[a2].rep], rep_w[zg[key].rep]);
                if (c2 <= 0) break;
                gwo[j + 1] = gwo[j];
                j--;
            }
            gwo[j + 1] = key;
        }
        int nqual = 0;
        int commit_start = dictcand.n;
        int *commit_rep = (int *)xmalloc((size_t)nzg * sizeof(int));
        len_t *commit_lam = (len_t *)xmalloc((size_t)nzg * sizeof(len_t));
        int ncommit = 0;
        {
            int i = 0;
            while (i < nzg) {
                int j = i;
                while (j < nzg && len_eq(zg[gwo[j]].z, zg[gwo[i]].z)) j++;
                len_t Z = zg[gwo[i]].z;
                /* 该组按权重升序：找第一个使 P*w >= Z 的 w */
                int found = -1;
                {
                    int lo = i, hi = j - 1, ans = -1;
                    while (lo <= hi) {
                        int mid = (lo + hi) / 2;
                        len_t Pw = len_mul_pow2(rep_w[zg[gwo[mid]].rep], Plog);
                        if (len_ge(Pw, Z)) { ans = mid; hi = mid - 1; }
                        else lo = mid + 1;
                    }
                    found = ans;
                }
                if (found >= 0) {
                    len_t lam = rep_w[zg[gwo[found]].rep];
                    for (int k = found; k < j; k++) {
                        int rp = zg[gwo[k]].rep;
                        if (len_ge(rep_w[rp], lam)) {
                            commit_rep[ncommit] = rp;
                            commit_lam[ncommit] = lam;
                            ncommit++;
                        }
                    }
                }
                i = j;
            }
        }
        nqual = ncommit;
        if (nqual * 2 >= nlive) {
            /* 提交 */
            char *dead = (char *)xcalloc((size_t)nrep, 1);
            for (int i = 0; i < ncommit; i++) {
                int rp = commit_rep[i];
                lambda_of_rep[rp] = commit_lam[i];
                lvec_push(&dictcand, commit_lam[i]);
                for (int k = rep_lo[rp]; k < rep_hi[rp]; k++) {
                    int req = rk[k].rid;
                    R->name[req] = commit_lam[i];
                }
                dead[rp] = 1;
            }
            int w2 = 0;
            for (int i = 0; i < nlive; i++) if (!dead[live[i]]) live[w2++] = live[i];
            nlive = w2;
            free(dead);
        } else {
            R->rejects++;
            dictcand.n = commit_start;
        }
        R->trials++;
        free(zg); free(grep); free(gwo); free(commit_rep); free(commit_lam);
        free(table);
    }

    /* 残留代表：各自以自己的权重为证书 */
    for (int i = 0; i < nlive; i++) {
        int rp = live[i];
        lambda_of_rep[rp] = rep_w[rp];
        lvec_push(&dictcand, rep_w[rp]);
        for (int k = rep_lo[rp]; k < rep_hi[rp]; k++) R->name[rk[k].rid] = rep_w[rp];
    }

    /* 词典：排序去重 */
    int nd = dictcand.n;
    len_t *dv = (len_t *)xmalloc((size_t)(nd ? nd : 1) * sizeof(len_t));
    for (int i = 0; i < nd; i++) dv[i] = dictcand.a[i];
    /* 长度排序（只比较，不转换成整数） */
    for (int i = 1; i < nd; i++) {
        len_t k = dv[i];
        int j = i - 1;
        while (j >= 0 && len_gt(dv[j], k)) { dv[j + 1] = dv[j]; j--; }
        dv[j + 1] = k;
    }
    int cnt = 0;
    for (int i = 0; i < nd; i++) if (cnt == 0 || !len_eq(dv[i], dv[cnt - 1])) dv[cnt++] = dv[i];
    R->ndict = cnt;
    R->dict = dv;
    for (int i = 0; i < nreq; i++) {
        R->dict_of_req[i] = -1;
        if (!len_finite(R->name[i])) continue;
        for (int k = 0; k < cnt; k++) if (len_eq(R->dict[k], R->name[i])) { R->dict_of_req[i] = k; break; }
    }
    free(rep_u); free(rep_v); free(rep_w); free(rep_lo); free(rep_hi);
    free(live); free(lambda_of_rep); free(rk);
    free(su); free(sv); free(sw); free(du_arr); free(dv_arr); free(didx); free(ddisc);
    lvec_free(&dictcand);
}

void routes_free(Routes *R)
{
    free(R->ru); free(R->rv); free(R->rw); free(R->reps_of_edge);
    free(R->name); free(R->dict); free(R->dict_of_req);
    memset(R, 0, sizeof(*R));
}