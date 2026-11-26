/* graph.c -- 图、Dijkstra、树上查询、MST（算法 1） */
#include "sssp.h"

void graph_init(Graph *g)
{
    memset(g, 0, sizeof(*g));
}
void graph_free(Graph *g)
{
    free(g->eu); free(g->ev); free(g->ew); free(g->eid);
    free(g->head); free(g->nxt); free(g->to); free(g->arc_e);
    memset(g, 0, sizeof(*g));
}
int graph_add_edge(Graph *g, int u, int v, len_t w, int eid)
{
    int e = g->m++;
    g->eu = (int *)xrealloc(g->eu, (size_t)g->m * sizeof(int));
    g->ev = (int *)xrealloc(g->ev, (size_t)g->m * sizeof(int));
    g->ew = (len_t *)xrealloc(g->ew, (size_t)g->m * sizeof(len_t));
    g->eid = (int *)xrealloc(g->eid, (size_t)g->m * sizeof(int));
    g->eu[e] = u; g->ev[e] = v; g->ew[e] = w; g->eid[e] = eid;
    return e;
}
void graph_finalize(Graph *g)
{
    int n = g->n;
    g->head = (int *)xmalloc((size_t)(n ? n : 1) * sizeof(int));
    for (int i = 0; i < n; i++) g->head[i] = -1;
    g->narc = 2 * g->m;
    g->to = (int *)xmalloc((size_t)(g->narc ? g->narc : 1) * sizeof(int));
    g->nxt = (int *)xmalloc((size_t)(g->narc ? g->narc : 1) * sizeof(int));
    g->arc_e = (int *)xmalloc((size_t)(g->narc ? g->narc : 1) * sizeof(int));
    for (int e = 0; e < g->m; e++) {
        int a = 2 * e, b = 2 * e + 1;
        g->to[a] = g->ev[e]; g->arc_e[a] = e; g->nxt[a] = g->head[g->eu[e]]; g->head[g->eu[e]] = a;
        g->to[b] = g->eu[e]; g->arc_e[b] = e; g->nxt[b] = g->head[g->ev[e]]; g->head[g->ev[e]] = b;
    }
}

void dist_init(Dist *D, int n)
{
    D->d = (len_t *)xmalloc((size_t)n * sizeof(len_t));
    D->par_vertex = (int *)xmalloc((size_t)n * sizeof(int));
    D->par_edge = (int *)xmalloc((size_t)n * sizeof(int));
    for (int i = 0; i < n; i++) { D->d[i] = len_inf(); D->par_vertex[i] = -1; D->par_edge[i] = -1; }
}
void dist_free(Dist *D) { free(D->d); free(D->par_vertex); free(D->par_edge); D->d = NULL; }

/* ---------------- Dijkstra ---------------- */

void dijkstra_heap(const Graph *g, int s, Dist *D)
{
    Heap h;
    heap_init(&h, g->n);
    for (int i = 0; i < g->n; i++) { D->d[i] = len_inf(); D->par_vertex[i] = -1; D->par_edge[i] = -1; }
    D->d[s] = len_zero();
    heap_push(&h, len_zero(), s);
    char *done = (char *)xcalloc((size_t)g->n, 1);
    for (;;) {
        len_t du;
        int u = heap_pop(&h, &du);
        if (u < 0) break;
        if (done[u]) continue;
        done[u] = 1;
        for (int a = g->head[u]; a != -1; a = g->nxt[a]) {
            int v = g->to[a];
            len_t nd = len_add(du, g->ew[g->arc_e[a]]);
            if (len_lt(nd, D->d[v])) {
                D->d[v] = nd;
                D->par_vertex[v] = u;
                D->par_edge[v] = g->arc_e[a];
                heap_decrease(&h, v, nd);
            }
        }
    }
    free(done);
    heap_free(&h);
}

void dijkstra_scan(const Graph *g, int s, Dist *D)
{
    int n = g->n;
    for (int i = 0; i < n; i++) { D->d[i] = len_inf(); D->par_vertex[i] = -1; D->par_edge[i] = -1; }
    char *done = (char *)xcalloc((size_t)n, 1);
    D->d[s] = len_zero();
    for (int it = 0; it < n; it++) {
        int u = -1;
        for (int v = 0; v < n; v++)
            if (!done[v] && len_finite(D->d[v]))
                if (u < 0 || len_lt(D->d[v], D->d[u])) u = v;
        if (u < 0) break;
        done[u] = 1;
        for (int a = g->head[u]; a != -1; a = g->nxt[a]) {
            int v = g->to[a];
            len_t nd = len_add(D->d[u], g->ew[g->arc_e[a]]);
            if (len_lt(nd, D->d[v])) { D->d[v] = nd; D->par_vertex[v] = u; D->par_edge[v] = g->arc_e[a]; }
        }
    }
    free(done);
}

void dijkstra_radius(const Graph *g, int s, len_t R, Dist *D)
{
    Heap h;
    heap_init(&h, g->n);
    for (int i = 0; i < g->n; i++) { D->d[i] = len_inf(); D->par_vertex[i] = -1; D->par_edge[i] = -1; }
    D->d[s] = len_zero();
    heap_push(&h, len_zero(), s);
    char *done = (char *)xcalloc((size_t)g->n, 1);
    for (;;) {
        len_t du;
        int u = heap_pop(&h, &du);
        if (u < 0) break;
        if (len_gt(du, R)) break;
        if (done[u]) continue;
        done[u] = 1;
        for (int a = g->head[u]; a != -1; a = g->nxt[a]) {
            int v = g->to[a];
            len_t nd = len_add(du, g->ew[g->arc_e[a]]);
            if (len_lt(nd, D->d[v]) && len_le(nd, R)) {
                D->d[v] = nd;
                D->par_vertex[v] = u;
                D->par_edge[v] = g->arc_e[a];
                heap_decrease(&h, v, nd);
            }
        }
    }
    free(done);
    heap_free(&h);
}

/* ---------------- 树查询 ---------------- */

void tree_build(Tree *T, int n, int root, const int *tu, const int *tv, const len_t *tw, int nt)
{
    memset(T, 0, sizeof(*T));
    T->n = n; T->root = root;
    int lg = 1;
    while ((1 << lg) < n) lg++;
    T->lg = lg + 1;
    T->depth = (int *)xmalloc((size_t)n * sizeof(int));
    T->dist_root = (len_t *)xmalloc((size_t)n * sizeof(len_t));
    T->up = (int *)xmalloc((size_t)T->lg * (size_t)n * sizeof(int));
    T->up_max = (len_t *)xmalloc((size_t)T->lg * (size_t)n * sizeof(len_t));

    int *head = (int *)xmalloc((size_t)n * sizeof(int));
    for (int i = 0; i < n; i++) head[i] = -1;
    int *nxt = (int *)xmalloc((size_t)(2 * nt + 1) * sizeof(int));
    int *to = (int *)xmalloc((size_t)(2 * nt + 1) * sizeof(int));
    len_t *ww = (len_t *)xmalloc((size_t)(2 * nt + 1) * sizeof(len_t));
    for (int i = 0; i < nt; i++) {
        int a = 2 * i, b = 2 * i + 1;
        to[a] = tv[i]; ww[a] = tw[i]; nxt[a] = head[tu[i]]; head[tu[i]] = a;
        to[b] = tu[i]; ww[b] = tw[i]; nxt[b] = head[tv[i]]; head[tv[i]] = b;
    }
    int *queue = (int *)xmalloc((size_t)n * sizeof(int));
    char *vis = (char *)xcalloc((size_t)n, 1);
    int qh = 0, qt = 0;
    queue[qt++] = root; vis[root] = 1;
    T->depth[root] = 0; T->dist_root[root] = len_zero();
    T->up[(size_t)0 * n + root] = root;
    T->up_max[(size_t)0 * n + root] = len_zero();
    while (qh < qt) {
        int u = queue[qh++];
        for (int a = head[u]; a != -1; a = nxt[a]) {
            int v = to[a];
            if (vis[v]) continue;
            vis[v] = 1;
            T->depth[v] = T->depth[u] + 1;
            T->dist_root[v] = len_add(T->dist_root[u], ww[a]);
            T->up[(size_t)0 * n + v] = u;
            T->up_max[(size_t)0 * n + v] = ww[a];
            queue[qt++] = v;
        }
    }
    for (int k = 1; k < T->lg; k++) {
        for (int v = 0; v < n; v++) {
            int mid = T->up[(size_t)(k - 1) * n + v];
            int anc = T->up[(size_t)(k - 1) * n + mid];
            T->up[(size_t)k * n + v] = anc;
            len_t m1 = T->up_max[(size_t)(k - 1) * n + v];
            len_t m2 = T->up_max[(size_t)(k - 1) * n + mid];
            T->up_max[(size_t)k * n + v] = len_max(m1, m2);
        }
    }
    free(head); free(nxt); free(to); free(ww); free(queue); free(vis);
}
void tree_free(Tree *T)
{
    free(T->depth); free(T->dist_root); free(T->up); free(T->up_max);
    memset(T, 0, sizeof(*T));
}
int tree_lca(Tree *T, int u, int v)
{
    int n = T->n;
    if (T->depth[u] < T->depth[v]) { int t = u; u = v; v = t; }
    int diff = T->depth[u] - T->depth[v];
    for (int k = 0; diff; k++, diff >>= 1)
        if (diff & 1) u = T->up[(size_t)k * n + u];
    if (u == v) return u;
    for (int k = T->lg - 1; k >= 0; k--) {
        if (T->up[(size_t)k * n + u] != T->up[(size_t)k * n + v]) {
            u = T->up[(size_t)k * n + u];
            v = T->up[(size_t)k * n + v];
        }
    }
    return T->up[(size_t)0 * n + u];
}
len_t tree_distance(Tree *T, int u, int v)
{
    if (u == v) return len_zero();
    {
        int w = tree_lca(T, u, v);
        /* D(u,v) = dist(u) + dist(v) - 2 dist(w) 用带符号对表示 */
        sgn_t s = sgn_from(len_add(T->dist_root[u], T->dist_root[v]));
        sgn_t t = sgn_from(len_mul_pow2(T->dist_root[w], 1));
        sgn_t r = sgn_sub(s, t);
        return len_max(r.pos, len_zero());
    }
}
len_t tree_bottleneck(Tree *T, int u, int v)
{
    int n = T->n;
    len_t best = len_zero();
    if (u == v) return best;
    if (T->depth[u] < T->depth[v]) { int t = u; u = v; v = t; }
    int diff = T->depth[u] - T->depth[v];
    for (int k = 0; diff; k++, diff >>= 1) {
        if (diff & 1) {
            best = len_max(best, T->up_max[(size_t)k * n + u]);
            u = T->up[(size_t)k * n + u];
        }
    }
    if (u == v) return best;
    for (int k = T->lg - 1; k >= 0; k--) {
        if (T->up[(size_t)k * n + u] != T->up[(size_t)k * n + v]) {
            best = len_max(best, T->up_max[(size_t)k * n + u]);
            best = len_max(best, T->up_max[(size_t)k * n + v]);
            u = T->up[(size_t)k * n + u];
            v = T->up[(size_t)k * n + v];
        }
    }
    best = len_max(best, T->up_max[(size_t)0 * n + u]);
    best = len_max(best, T->up_max[(size_t)0 * n + v]);
    return best;
}

/* ---------------- 拓扑排序 ---------------- */

int topo_sort(int n, const ivec *adj, int *order)
{
    int *indeg = (int *)xcalloc((size_t)(n ? n : 1), sizeof(int));
    for (int u = 0; u < n; u++)
        for (int i = 0; i < adj[u].n; i++) indeg[adj[u].a[i]]++;
    int *q = (int *)xmalloc((size_t)(n ? n : 1) * sizeof(int));
    int qh = 0, qt = 0;
    for (int u = 0; u < n; u++) if (indeg[u] == 0) q[qt++] = u;
    int cnt = 0;
    while (qh < qt) {
        int u = q[qh++];
        order[cnt++] = u;
        for (int i = 0; i < adj[u].n; i++) {
            int v = adj[u].a[i];
            if (--indeg[v] == 0) q[qt++] = v;
        }
    }
    free(indeg); free(q);
    return cnt == n;
}

/* ---------------- Kruskal ---------------- */

typedef struct { len_t w; int e; } WE;
static int cmp_we(const void *a, const void *b)
{
    const WE *x = (const WE *)a, *y = (const WE *)b;
    int c = len_cmp(x->w, y->w);
    if (c) return c;
    return x->e - y->e;
}
int mst_kruskal(const Graph *g, int *out_edges)
{
    int m = g->m;
    WE *arr = (WE *)xmalloc((size_t)(m ? m : 1) * sizeof(WE));
    for (int i = 0; i < m; i++) { arr[i].w = g->ew[i]; arr[i].e = i; }
    qsort(arr, (size_t)m, sizeof(WE), cmp_we);
    DSU d; dsu_init(&d, g->n);
    int cnt = 0;
    for (int i = 0; i < m; i++) {
        if (dsu_unite(&d, g->eu[arr[i].e], g->ev[arr[i].e])) out_edges[cnt++] = arr[i].e;
    }
    dsu_free(&d);
    free(arr);
    return cnt;
}

/* ---------------- 算法 1：有界生成树生长趟 ---------------- */

typedef struct { len_t w; int e; int v; } Bound;

/* 以 (w, e) 为键的小顶堆，元素是"外部顶点" */
typedef struct {
    Bound *b;
    int n, cap;
    int *pos;   /* 顶点 -> 堆下标 */
} BHeap;

static void bh_init(BHeap *h, int nver)
{
    h->b = NULL; h->n = 0; h->cap = 0;
    h->pos = (int *)xmalloc((size_t)(nver ? nver : 1) * sizeof(int));
    for (int i = 0; i < nver; i++) h->pos[i] = -1;
}
static void bh_free(BHeap *h) { free(h->b); free(h->pos); h->b = NULL; h->pos = NULL; }
static int bkey_lt(Bound a, Bound c)
{
    int r = len_cmp(a.w, c.w);
    if (r) return r < 0;
    return a.e < c.e;
}
static void bh_up(BHeap *h, int i)
{
    while (i > 0) {
        int p = (i - 1) / 2;
        if (bkey_lt(h->b[i], h->b[p])) {
            Bound t = h->b[i]; h->b[i] = h->b[p]; h->b[p] = t;
            h->pos[h->b[i].v] = i; h->pos[h->b[p].v] = p;
            i = p;
        } else break;
    }
}
static void bh_down(BHeap *h, int i)
{
    for (;;) {
        int l = 2 * i + 1, r = l + 1, m = i;
        if (l < h->n && bkey_lt(h->b[l], h->b[m])) m = l;
        if (r < h->n && bkey_lt(h->b[r], h->b[m])) m = r;
        if (m == i) break;
        Bound t = h->b[i]; h->b[i] = h->b[m]; h->b[m] = t;
        h->pos[h->b[i].v] = i; h->pos[h->b[m].v] = m;
        i = m;
    }
}
static void bh_offer(BHeap *h, Bound x)
{
    int i = h->pos[x.v];
    if (i >= 0) {
        if (bkey_lt(x, h->b[i])) { h->b[i] = x; bh_up(h, i); }
        return;
    }
    if (h->n == h->cap) {
        h->cap = h->cap ? h->cap * 2 : 16;
        h->b = (Bound *)xrealloc(h->b, (size_t)h->cap * sizeof(Bound));
    }
    i = h->n++;
    h->b[i] = x; h->pos[x.v] = i;
    bh_up(h, i);
}
static int bh_pop(BHeap *h, Bound *out)
{
    if (h->n == 0) return 0;
    *out = h->b[0];
    h->pos[out->v] = -1;
    h->n--;
    if (h->n > 0) {
        h->b[0] = h->b[h->n];
        h->pos[h->b[0].v] = 0;
        bh_down(h, 0);
    }
    return 1;
}

/* 一趟生长：quotient 用邻接表（qu 顶点集，边带原始下标） */
typedef struct {
    int nq;
    int *head, *nxt, *to, *eid;
    len_t *w;
    int narc;
} QGraph;

static void qg_build(QGraph *q, const Graph *g, int nq, const int *comp)
{
    memset(q, 0, sizeof(*q));
    q->nq = nq;
    q->head = (int *)xmalloc((size_t)(nq ? nq : 1) * sizeof(int));
    for (int i = 0; i < nq; i++) q->head[i] = -1;
    q->nxt = (int *)xmalloc((size_t)(2 * g->m + 1) * sizeof(int));
    q->to = (int *)xmalloc((size_t)(2 * g->m + 1) * sizeof(int));
    q->eid = (int *)xmalloc((size_t)(2 * g->m + 1) * sizeof(int));
    q->w = (len_t *)xmalloc((size_t)(2 * g->m + 1) * sizeof(len_t));
    q->narc = 0;
    for (int e = 0; e < g->m; e++) {
        int a = comp[g->eu[e]], b = comp[g->ev[e]];
        if (a == b) continue;
        int i = q->narc++;
        q->to[i] = b; q->eid[i] = e; q->w[i] = g->ew[e]; q->nxt[i] = q->head[a]; q->head[a] = i;
        int j = q->narc++;
        q->to[j] = a; q->eid[j] = e; q->w[j] = g->ew[e]; q->nxt[j] = q->head[b]; q->head[b] = j;
    }
}
static void qg_free(QGraph *q)
{
    free(q->head); free(q->nxt); free(q->to); free(q->eid); free(q->w);
    memset(q, 0, sizeof(*q));
}

static int mst_one_pass(const Graph *g, int nq, QGraph *q, int limit, int *visited_out, ivec *chosen)
{
    /* visited: 0 = 未访问；>0 = 本趟生长的编号（-1 表示更早的生长） */
    int *gid = (int *)xcalloc((size_t)(nq ? nq : 1), sizeof(int));
    int growth = 0;
    BHeap h;
    bh_init(&h, nq);
    for (int start = 0; start < nq; start++) {
        if (gid[start]) continue;
        growth++;
        bh_free(&h);
        bh_init(&h, nq);
        int u = start;
        for (;;) {
            gid[u] = growth;
            for (int a = q->head[u]; a != -1; a = q->nxt[a]) {
                int v = q->to[a];
                if (gid[v] == growth) continue;         /* 本生长内部 */
                Bound b; b.w = q->w[a]; b.e = g->eid[q->eid[a]]; b.v = v;
                bh_offer(&h, b);
            }
            if (h.n == 0 || (limit >= 0 && h.n > limit)) break;
            Bound best;
            bh_pop(&h, &best);
            ivec_push(chosen, best.e);
            if (gid[best.v] && gid[best.v] < growth) break;  /* 连到更早的生长 */
            u = best.v;
        }
    }
    bh_free(&h);
    int nvis = growth;
    if (visited_out) *visited_out = nvis;
    free(gid);
    return growth;
}

int mst_growth_passes(const Graph *g, int *out_edges)
{
    int n = g->n, m = g->m;
    if (n <= 1) return 0;
    int Lk = 1;
    while ((1 << (Lk - 1)) < (n + 2)) Lk++;         /* ceil(log2(n+2)) */
    int B1 = (int)ceil(log2((double)Lk));
    if (B1 < 2) B1 = 2;
    int B2 = Lk;
    int limits[3] = { B1, B2, -1 };                  /* -1 = 无穷 */

    DSU d; dsu_init(&d, n);
    ivec chosen; ivec_init(&chosen);
    for (int pass = 0; pass < 3; pass++) {
        /* 构造 quotient：comp[] 为紧凑分量编号 */
        int *raw = (int *)xmalloc((size_t)n * sizeof(int));
        int *comp = (int *)xmalloc((size_t)n * sizeof(int));
        int nq = 0;
        for (int v = 0; v < n; v++) raw[v] = dsu_find(&d, v);
        int *map = (int *)xmalloc((size_t)n * sizeof(int));
        for (int v = 0; v < n; v++) map[v] = -1;
        for (int v = 0; v < n; v++) {
            if (map[raw[v]] < 0) map[raw[v]] = nq++;
            comp[v] = map[raw[v]];
        }
        QGraph q;
        qg_build(&q, g, nq, comp);
        int before = chosen.n;
        mst_one_pass(g, nq, &q, limits[pass], NULL, &chosen);
        qg_free(&q);
        /* 收缩本趟选择的边 */
        for (int i = before; i < chosen.n; i++) {
            int e = chosen.a[i];
            dsu_unite(&d, g->eu[e], g->ev[e]);
        }
        free(raw); free(comp); free(map);
        if (nq == 1) break;
    }
    int cnt = 0;
    /* 去重（同一原始边可能被多次选择？论文保证不会；这里保险） */
    char *seen = (char *)xcalloc((size_t)(m ? m : 1), 1);
    for (int i = 0; i < chosen.n; i++) {
        int e = chosen.a[i];
        if (!seen[e]) { seen[e] = 1; out_edges[cnt++] = e; }
    }
    free(seen);
    ivec_free(&chosen);
    dsu_free(&d);
    return cnt;
}
