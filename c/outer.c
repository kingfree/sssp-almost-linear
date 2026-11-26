/* outer.c -- 第 8.1 节 外层状态图、阈值分量与层级组装
 *
 * 状态图：以归属坐标 J 为坐标、壳层比 K = 16 的双壳构造。
 *   - 顶点的固有壳层 i(v) 由 a*K^i <= J_v < a*K^(i+1) 确定（a 为近邻分量的最小正 J）；
 *   - 内向 j 号窗口含固有壳层 j 与 j+1，外向含 j-1 与 j；每个顶点至多两个内向、
 *     两个外向副本；
 *   - 窗口内的近邻边双向；向内降到更低固有壳层 = "近邻内转移"；远距径向边做
 *     "远距转移"（从高坐标端点的每个内向副本降到低端点的固有内向副本，反向则
 *     从低端点的外向副本升到高端点的固有外向副本）；相位切换为零代价；
 *   - 宏观边从尾的任意普通副本重置到头的固有内向副本；
 *   - 额外的归属树状态 v_in / v_out 让 J 不必在归属树上单调。
 *
 * 阈值分量：块 = 窗口内双向近邻边的连通分量；用块 DAG 上的两趟传播得到
 * I_0 与 O_0（忽略径向权重的最小宏观权重）、以及 I_ 与 O_ 的星号版本（整条路径的瓶颈），
 * θ = max(星号版 I, 星号版 O)；再由"块 MST 边 + 状态到枢纽的 θ 连接 + 枢纽走廊连接"
 * 给出阈值强连通分量，见引理 7.3。
 *
 * 层级：把连接记录按倍增网格的类 [a, 2a) 同时合并（引理 5.2 的阈值 SCC 刻画），
 * 每个类生成一个高为 a 的内部节点；辅助顶点只用于识别核，最后被压缩掉。
 */
#include "sssp.h"

void sgraph_init(SGraph *S) { memset(S, 0, sizeof(*S)); }
void sgraph_free(SGraph *S)
{
    free(S->proj); free(S->phase); free(S->shell);
    free(S->native_in); free(S->native_out);
    free(S->arc_head); free(S->arc_next); free(S->arc_to); free(S->arc_len); free(S->arc_macro);
    free(S->vhead);
    memset(S, 0, sizeof(*S));
}
int sgraph_add_state(SGraph *S, int proj, int phase, int shell)
{
    int id = S->nstates++;
    S->proj = (int *)xrealloc(S->proj, (size_t)S->nstates * sizeof(int));
    S->phase = (int *)xrealloc(S->phase, (size_t)S->nstates * sizeof(int));
    S->shell = (int *)xrealloc(S->shell, (size_t)S->nstates * sizeof(int));
    S->proj[id] = proj; S->phase[id] = phase; S->shell[id] = shell;
    return id;
}
void sgraph_add_arc(SGraph *S, int u, int v, len_t w)
{
    int a = S->narc++;
    S->arc_head = (int *)xrealloc(S->arc_head, (size_t)S->narc * sizeof(int));
    S->arc_next = (int *)xrealloc(S->arc_next, (size_t)S->narc * sizeof(int));
    S->arc_to = (int *)xrealloc(S->arc_to, (size_t)S->narc * sizeof(int));
    S->arc_len = (len_t *)xrealloc(S->arc_len, (size_t)S->narc * sizeof(len_t));
    S->arc_head[a] = u; S->arc_to[a] = v; S->arc_len[a] = w;
    S->arc_next[a] = -1;
    S->arc_macro = (int *)xrealloc(S->arc_macro, (size_t)S->narc * sizeof(int));
    S->arc_macro[a] = S->cur_macro;
}
void sgraph_finalize(SGraph *S)
{
    S->vhead = (int *)xmalloc((size_t)(S->nstates ? S->nstates : 1) * sizeof(int));
    for (int i = 0; i < S->nstates; i++) S->vhead[i] = -1;
    /* 把弧按尾分组（保持插入序的逆序无妨） */
    for (int a = 0; a < S->narc; a++) {
        int u = S->arc_head[a];
        S->arc_next[a] = S->vhead[u];
        S->vhead[u] = a;
    }
}

/* ---------------- 外层状态图构造 ---------------- */

#define K_SHELL 16

static int find_or_make(SGraph *S, int *cache, int proj, int phase, int shell)
{
    if (*cache >= 0) return *cache;
    *cache = sgraph_add_state(S, proj, phase, shell);
    return *cache;
}

/* 构造外层状态图。nears[] 标记那些"两端坐标均为正且比值 <= K^2"的径向边 */
void build_outer_states(SGraph *S, StateMap *M, const Graph *g, const Owner *O,
                        const len_t *J, const len_t *de, const int *radial,
                        const int *near, const int *mst_edges, int nmst)
{
    int n = g->n;
    memset(S, 0, sizeof(*S));
    memset(M, 0, sizeof(*M));
    M->v_in = (int *)xmalloc((size_t)n * sizeof(int));
    M->v_inlo = (int *)xmalloc((size_t)n * sizeof(int));
    M->v_out = (int *)xmalloc((size_t)n * sizeof(int));
    M->v_outhi = (int *)xmalloc((size_t)n * sizeof(int));
    M->special_in = (int *)xmalloc((size_t)n * sizeof(int));
    M->special_out = (int *)xmalloc((size_t)n * sizeof(int));
    M->shell_of = (int *)xmalloc((size_t)n * sizeof(int));
    M->base_of = (double *)xmalloc((size_t)n * sizeof(double));
    for (int v = 0; v < n; v++) {
        M->v_in[v] = M->v_inlo[v] = M->v_out[v] = M->v_outhi[v] = -1;
        M->special_in[v] = M->special_out[v] = -1;
        M->shell_of[v] = -1;
        M->base_of[v] = 0.0;
    }
    char *ispiv = (char *)xcalloc((size_t)n, 1);
    for (int i = 0; i < O->npiv; i++) ispiv[O->piv_vertex[i]] = 1;

    /* 1) 近邻分量：径向边两端坐标均为正且比值 <= K^2 */
    DSU nd; dsu_init(&nd, n);
    for (int e = 0; e < g->m; e++) {
        if (!radial[e] || !len_finite(de[e])) continue;
        if (!near[e]) continue;
        dsu_unite(&nd, g->eu[e], g->ev[e]);
    }
    /* 每个近邻分量的最小正 J 作为壳层底 */
    double *cbase = (double *)xmalloc((size_t)n * sizeof(double));
    for (int v = 0; v < n; v++) cbase[v] = 0.0;
    for (int v = 0; v < n; v++) {
        if (isinf(J[v].v) || J[v].v <= 0.0) continue;
        int r = dsu_find(&nd, v);
        if (cbase[r] == 0.0 || J[v].v < cbase[r]) cbase[r] = J[v].v;
    }
    for (int v = 0; v < n; v++) M->base_of[v] = (isinf(J[v].v) || J[v].v <= 0.0) ? 0.0 : cbase[dsu_find(&nd, v)];
    dsu_free(&nd);
    free(cbase);

    /* 2) 固有壳层 i(v) */
    for (int v = 0; v < n; v++) {
        if (ispiv[v]) continue;
        if (isinf(J[v].v) || J[v].v <= 0.0) continue;
        double a = M->base_of[v];
        if (a <= 0.0) continue;
        double r = J[v].v / a;
        int i = 0;
        double p = 1.0;
        while (p * K_SHELL <= r + 1e-12) { p *= K_SHELL; i++; }
        /* a*K^i <= J_v < a*K^(i+1) */
        while (a * p > J[v].v * (1.0 + 1e-12) && i > 0) { p /= K_SHELL; i--; }
        M->shell_of[v] = i;
    }

    /* 3) 建状态：枢纽 = 一内向一外向；正坐标 = 至多两个内向 + 两个外向 */
    int *piv_in = (int *)xmalloc((size_t)(O->npiv ? O->npiv : 1) * sizeof(int));
    int *piv_out = (int *)xmalloc((size_t)(O->npiv ? O->npiv : 1) * sizeof(int));
    for (int i = 0; i < O->npiv; i++) {
        int v = O->piv_vertex[i];
        piv_in[i] = sgraph_add_state(S, v, PH_IN, -1);
        piv_out[i] = sgraph_add_state(S, v, PH_OUT, -1);
        M->v_in[v] = piv_in[i];
        M->v_out[v] = piv_out[i];
    }
    for (int v = 0; v < n; v++) {
        if (ispiv[v]) continue;
        int i = M->shell_of[v];
        if (i < 0) {
            /* 零坐标或不可达：当作枢纽式单副本 */
            M->v_in[v] = sgraph_add_state(S, v, PH_IN, -1);
            M->v_out[v] = sgraph_add_state(S, v, PH_OUT, -1);
            continue;
        }
        M->v_in[v] = sgraph_add_state(S, v, PH_IN, i);
        M->v_out[v] = sgraph_add_state(S, v, PH_OUT, i);
        if (i - 1 >= 0) M->v_inlo[v] = sgraph_add_state(S, v, PH_IN, i - 1);
        M->v_outhi[v] = sgraph_add_state(S, v, PH_OUT, i + 1);
    }
    S->native_in = (int *)xmalloc((size_t)n * sizeof(int));
    S->native_out = (int *)xmalloc((size_t)n * sizeof(int));
    for (int v = 0; v < n; v++) { S->native_in[v] = M->v_in[v]; S->native_out[v] = M->v_out[v]; }
    /* 归属树特殊状态 */
    for (int v = 0; v < n; v++) {
        if (ispiv[v]) continue;
        M->special_in[v] = sgraph_add_state(S, v, PH_TIN, -1);
        M->special_out[v] = sgraph_add_state(S, v, PH_TOUT, -1);
    }

    /* 4) 窗口内近邻边（双向） */
    for (int e = 0; e < g->m; e++) {
        if (!radial[e] || !len_finite(de[e]) || !near[e]) continue;
        int x = g->eu[e], y = g->ev[e];
        if (ispiv[x] || ispiv[y]) continue;
        int ix = M->shell_of[x], iy = M->shell_of[y];
        if (ix < 0 || iy < 0) continue;
        len_t w = g->ew[e];
        /* 内向：窗口 j 同时覆盖 x 与 y 的固有壳层 */
        for (int j = (ix < iy ? ix : iy) - 1; j <= (ix > iy ? ix : iy); j++) {
            if (j < 0) continue;
            int sx = (j == ix) ? M->v_in[x] : (j == ix - 1 ? M->v_inlo[x] : -1);
            int sy = (j == iy) ? M->v_in[y] : (j == iy - 1 ? M->v_inlo[y] : -1);
            if (sx >= 0 && sy >= 0) { sgraph_add_arc(S, sx, sy, w); sgraph_add_arc(S, sy, sx, w); }
        }
        /* 外向：窗口 j 覆盖固有壳层 j-1 与 j */
        for (int j = (ix < iy ? ix : iy); j <= (ix > iy ? ix : iy) + 1; j++) {
            if (j < 0) continue;
            int sx = (j == ix) ? M->v_out[x] : (j == ix + 1 ? M->v_outhi[x] : -1);
            int sy = (j == iy) ? M->v_out[y] : (j == iy + 1 ? M->v_outhi[y] : -1);
            if (sx >= 0 && sy >= 0) { sgraph_add_arc(S, sx, sy, w); sgraph_add_arc(S, sy, sx, w); }
        }
        /* 近邻内转移：内向降到更低固有壳层；外向升到更高固有壳层 */
        if (iy < ix) sgraph_add_arc(S, M->v_in[x], M->v_in[y], w);
        else sgraph_add_arc(S, M->v_in[y], M->v_in[x], w);
        if (iy > ix) sgraph_add_arc(S, M->v_out[x], M->v_out[y], w);
        else if (ix > iy) sgraph_add_arc(S, M->v_out[y], M->v_out[x], w);
    }

    /* 5) 远距径向边与宏观边 */
    for (int e = 0; e < g->m; e++) {
        if (!radial[e] || !len_finite(de[e])) continue;
        if (near[e]) continue;
        int x = g->eu[e], y = g->ev[e];
        len_t w = g->ew[e];
        double Jx = J[x].v, Jy = J[y].v;
        if (ispiv[x] || ispiv[y]) {
            /* 枢纽端点：径向弧指向枢纽内向副本，外向弧离开枢纽外向副本 */
            int p = ispiv[x] ? x : y;
            int v = ispiv[x] ? y : x;
            int pin = S->native_in[p], pout = S->native_out[p];
            if (M->v_in[v] >= 0) {
                sgraph_add_arc(S, M->v_in[v], pin, w);
                if (M->v_inlo[v] >= 0) sgraph_add_arc(S, M->v_inlo[v], pin, w);
            }
            sgraph_add_arc(S, pout, S->native_out[v], w);
            continue;
        }
        int hi = (Jx >= Jy) ? x : y;
        int lo = (Jx >= Jy) ? y : x;
        /* 远距内向：高端的每个内向副本 -> 低端的固有内向副本 */
        sgraph_add_arc(S, M->v_in[hi], S->native_in[lo], len_add(w, de[e]));
        if (M->v_inlo[hi] >= 0) sgraph_add_arc(S, M->v_inlo[hi], S->native_in[lo], len_add(w, de[e]));
        /* 远距外向：低端的每个外向副本 -> 高端的固有外向副本 */
        sgraph_add_arc(S, M->v_out[lo], S->native_out[hi], len_add(w, de[e]));
        if (M->v_outhi[lo] >= 0) sgraph_add_arc(S, M->v_outhi[lo], S->native_out[hi], len_add(w, de[e]));
    }
    S->cur_macro = 1;
    for (int e = 0; e < g->m; e++) {
        if (radial[e]) continue;
        int x = g->eu[e], y = g->ev[e];
        len_t w = g->ew[e];
        /* 宏观边：尾的每个普通副本 -> 头的固有内向副本 */
        for (int dir = 0; dir < 2; dir++) {
            int t = dir ? y : x, h = dir ? x : y;
            int tin[6], ntin = 0;
            if (M->v_in[t] >= 0) tin[ntin++] = M->v_in[t];
            if (M->v_inlo[t] >= 0) tin[ntin++] = M->v_inlo[t];
            if (M->v_out[t] >= 0) tin[ntin++] = M->v_out[t];
            if (M->v_outhi[t] >= 0) tin[ntin++] = M->v_outhi[t];
            for (int k = 0; k < ntin; k++)
                if (S->native_in[h] >= 0) sgraph_add_arc(S, tin[k], S->native_in[h], w);
        }
    }

    S->cur_macro = 0;

    /* 6) 相位切换（内向副本 -> 固有外向副本） */
    for (int v = 0; v < n; v++) {
        if (M->v_in[v] >= 0 && M->v_out[v] >= 0) sgraph_add_arc(S, M->v_in[v], M->v_out[v], len_zero());
        if (M->v_inlo[v] >= 0 && M->v_out[v] >= 0) sgraph_add_arc(S, M->v_inlo[v], M->v_out[v], len_zero());
    }

    /* 7) 归属树特殊状态 */
    for (int i = 0; i < nmst; i++) {
        int e = mst_edges[i];
        int a = g->eu[e], b = g->ev[e];
        len_t w = g->ew[e];
        /* 沿归属树（以本 MST 为树、根为 0）向内走：a -> b 与 b -> a 都要加，
           方向由树结构决定；用两侧度数法简单处理：把 v_in 连到相邻顶点的 v_in/pivot_in */
        if (M->special_in[a] >= 0) {
            if (M->special_in[b] >= 0) sgraph_add_arc(S, M->special_in[a], M->special_in[b], w);
            else if (S->native_in[b] >= 0) sgraph_add_arc(S, M->special_in[a], S->native_in[b], w);
        }
        if (M->special_in[b] >= 0) {
            if (M->special_in[a] >= 0) sgraph_add_arc(S, M->special_in[b], M->special_in[a], w);
            else if (S->native_in[a] >= 0) sgraph_add_arc(S, M->special_in[b], S->native_in[a], w);
        }
        if (M->special_out[a] >= 0) {
            if (M->special_out[b] >= 0) sgraph_add_arc(S, M->special_out[a], M->special_out[b], w);
            else if (S->native_out[b] >= 0) sgraph_add_arc(S, S->native_out[b], M->special_out[a], w);
        }
        if (M->special_out[b] >= 0) {
            if (M->special_out[a] >= 0) sgraph_add_arc(S, M->special_out[b], M->special_out[a], w);
            else if (S->native_out[a] >= 0) sgraph_add_arc(S, S->native_out[a], M->special_out[b], w);
        }
    }
    for (int v = 0; v < n; v++) {
        if (M->special_in[v] < 0) continue;
        /* 普通内向副本 -> 特殊内向（零）；特殊外向 -> 普通固有外向（零） */
        if (M->v_in[v] >= 0) sgraph_add_arc(S, M->v_in[v], M->special_in[v], len_zero());
        if (M->v_inlo[v] >= 0) sgraph_add_arc(S, M->v_inlo[v], M->special_in[v], len_zero());
        if (M->v_out[v] >= 0) sgraph_add_arc(S, M->special_out[v], M->v_out[v], len_zero());
        if (M->v_outhi[v] >= 0) sgraph_add_arc(S, M->special_out[v], M->v_outhi[v], len_zero());
    }

    sgraph_finalize(S);
    free(piv_in); free(piv_out); free(ispiv);
}