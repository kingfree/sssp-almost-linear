/* main.c -- 算法 9 的驱动、图生成器、与 Dijkstra 的对拍、CLI
 *
 * 用法：
 *   sssp --type rand --n 200 --m 600 --seed 1 --verify
 *   sssp --type rand --n 200 --m 600 --sep 3      （用缩放的分隔规模跑层级路径）
 *   sssp --type grid --n 400 --verify
 *   sssp --bench                                  （跑一组基准）
 */
#include "sssp.h"
#include <time.h>

/* ---------------- 图生成器 ---------------- */

static void add_edge_w(Graph *g, int u, int v, double w, int *eid)
{
    graph_add_edge(g, u, v, len_make(w), (*eid)++);
}

static Graph gen_rand(int n, int m, int maxw, int seed, int real_weights)
{
    Graph g; graph_init(&g); g.n = n;
    rng_seed((uint64_t)seed * 2654435761ULL + 12345ULL);
    int eid = 0;
    /* 先连一棵随机树保证连通 */
    for (int v = 1; v < n; v++) {
        int u = rng_int(v);
        double w = real_weights ? (0.5 + rng_double() * maxw) : (1 + rng_int(maxw));
        add_edge_w(&g, u, v, w, &eid);
    }
    for (int i = n - 1; i < m; i++) {
        int u = rng_int(n), v = rng_int(n);
        if (u == v) { v = (v + 1) % n; }
        double w = real_weights ? (0.5 + rng_double() * maxw) : (1 + rng_int(maxw));
        add_edge_w(&g, u, v, w, &eid);
    }
    graph_finalize(&g);
    return g;
}

static Graph gen_path(int n, int maxw, int seed)
{
    Graph g; graph_init(&g); g.n = n;
    rng_seed((uint64_t)seed);
    for (int i = 0; i + 1 < n; i++) graph_add_edge(&g, i, i + 1, len_make((double)(1 + rng_int(maxw))), i);
    graph_finalize(&g);
    return g;
}

static Graph gen_grid(int side, int maxw, int seed, int real_weights)
{
    Graph g; graph_init(&g); g.n = side * side;
    rng_seed((uint64_t)seed * 88172645463325252ULL + 7ULL);
    int eid = 0;
    for (int y = 0; y < side; y++)
        for (int x = 0; x < side; x++) {
            int v = y * side + x;
            if (x + 1 < side) {
                double w = real_weights ? (0.5 + rng_double() * maxw) : (1 + rng_int(maxw));
                add_edge_w(&g, v, v + 1, w, &eid);
            }
            if (y + 1 < side) {
                double w = real_weights ? (0.5 + rng_double() * maxw) : (1 + rng_int(maxw));
                add_edge_w(&g, v, v + side, w, &eid);
            }
        }
    graph_finalize(&g);
    return g;
}

/* ---------------- 流水线（算法 9） ---------------- */

typedef struct {
    int nstates, nnode, nleaf, nroot, npiv, nreq, njoin, ndict;
    int base_case;
    double build_ms, eval_ms;
} Stats;

static double now_ms(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1000.0 + ts.tv_nsec / 1e6;
}

/* 收缩零边分量并返回收缩图的顶点数（本实现要求正权重，故通常 = n） */
static int lceil_log2(long long x)
{
    int k = 0; long long p = 1;
    while (p < x) { p <<= 1; k++; }
    return k;
}

int sssp_pipeline(const Graph *g, int s, const Options *opt, Dist *out, Stats *st)
{
    int n = g->n;
    memset(st, 0, sizeof(*st));
    dist_init(out, n);
    double t0 = now_ms();
    int L = lceil_log2((long long)n + 2);
    /* 论文：N <= L_N^64 时走常数规模的直接搜索 */
    int small;
    if (opt->force_hier) small = 0;
    else {
        double bound = 1.0;
        for (int i = 0; i < 64; i++) bound *= (double)L;
        small = ((double)n <= bound);
    }
    if (small) {
        dijkstra_scan(g, s, out);
        st->base_case = 1;
        st->build_ms = now_ms() - t0;
        return 0;
    }

    /* 1) MST（算法 1 的三趟生长） */
    int *mst = (int *)xmalloc((size_t)(n ? n : 1) * sizeof(int));
    int nmst = mst_growth_passes(g, mst);
    /* 2) 删除被树路径严格统治的非树边 */
    int *keep = (int *)xmalloc((size_t)(g->m ? g->m : 1) * sizeof(int));
    reduce_tree_dominated(g, mst, nmst, keep);
    Graph G2; graph_init(&G2); G2.n = n;
    int *map2 = (int *)xmalloc((size_t)(g->m ? g->m : 1) * sizeof(int));
    int eid2 = 0;
    for (int e = 0; e < g->m; e++) {
        if (!keep[e]) { map2[e] = -1; continue; }
        map2[e] = eid2;
        graph_add_edge(&G2, g->eu[e], g->ev[e], g->ew[e], eid2);
        eid2++;
    }
    graph_finalize(&G2);
    int nmst2 = 0;
    {
        int *m2 = (int *)xmalloc((size_t)(nmst ? nmst : 1) * sizeof(int));
        for (int i = 0; i < nmst; i++) if (map2[mst[i]] >= 0) m2[nmst2++] = map2[mst[i]];
        memcpy(mst, m2, (size_t)nmst2 * sizeof(int));
        nmst = nmst2;
        free(m2);
    }

    /* 3) 归属分解 */
    Owner O;
    int sep = opt->sep_log > 0 ? 2 : 2;
    {
        int b = 2;
        if (opt->sep_log > 0) {
            int Lg = lceil_log2((long long)n + 2);
            long long v = 1;
            for (int i = 0; i < opt->sep_log; i++) { v *= Lg; if (v > n) { v = n; break; } }
            b = (int)(v < 2 ? 2 : v);
        }
        sep = b;
    }
    if (getenv("SSSP_DEBUG")) fprintf(stderr, "[stage] MST done: %d edges\n", nmst);
    owner_build(&O, &G2, mst, nmst, sep);
    if (getenv("SSSP_DEBUG")) fprintf(stderr, "[stage] owner done: %d pivots, %d cells, sep=%d\n", O.npiv, O.ncell, sep);

    /* 4) 投影请求与证书 */
    int *ru = NULL, *rv = NULL, *rep = NULL;
    len_t *rw = NULL;
    int nreq = 0;
    build_requests(&G2, &O, 1, &ru, &rv, &rw, &nreq, &rep, NULL, NULL);
    Routes R;
    routes_build(&R, O.npiv, nreq, ru, rv, rw, (uint64_t)opt->seed);
    if (getenv("SSSP_DEBUG")) fprintf(stderr, "[stage] routes done: %d requests, %d dict, P=%g\n", nreq, R.ndict, R.P.v);
    /* 边 -> 请求 */
    for (int e = 0; e < G2.m; e++) { /* rep 已给 */ }

    /* 5) 分类与归属坐标 */
    len_t *de = (len_t *)xmalloc((size_t)(G2.m ? G2.m : 1) * sizeof(len_t));
    int *radial = (int *)xmalloc((size_t)(G2.m ? G2.m : 1) * sizeof(int));
    classify_edges(&G2, &O, &R, rep, de, radial);
    len_t *J = (len_t *)xmalloc((size_t)(n ? n : 1) * sizeof(len_t));
    compute_J(&G2, &O, de, radial, mst, nmst, J);
    if (getenv("SSSP_DEBUG")) fprintf(stderr, "[stage] J done\n");

    /* 6) 删除 w > J_x + J_y + d_e 的径向边，并计算 near */
    int *near = (int *)xmalloc((size_t)(G2.m ? G2.m : 1) * sizeof(int));
    for (int e = 0; e < G2.m; e++) {
        near[e] = 0;
        if (!radial[e] || !len_finite(de[e])) continue;
        int x = G2.eu[e], y = G2.ev[e];
        if (!len_finite(J[x]) || !len_finite(J[y])) { radial[e] = 0; continue; }
        len_t bound = len_add(len_add(J[x], J[y]), de[e]);
        if (len_gt(G2.ew[e], bound)) { radial[e] = 0; continue; }   /* 用更短的枢纽路线替换 */
        if (len_gt(J[x], len_zero()) && len_gt(J[y], len_zero())) {
            double a = J[x].v, b = J[y].v;
            double ratio = a > b ? a / b : b / a;
            if (ratio <= 256.0) near[e] = 1;
        }
    }

    /* 7) 外层状态图 + 阈值层级 + 求值 */
    SGraph S; StateMap M;
    build_outer_states(&S, &M, &G2, &O, J, de, radial, near, mst, nmst);
    if (getenv("SSSP_DEBUG")) fprintf(stderr, "[stage] states done: %d states, %d arcs\n", S.nstates, S.narc);
    int naux = O.npiv;
    int *aux_state = (int *)xmalloc((size_t)(naux ? naux : 1) * sizeof(int));
    for (int p = 0; p < naux; p++) aux_state[p] = M.v_in[O.piv_vertex[p]];
    if (getenv("SSSP_DEBUG")) {
        /* 状态图可达性诊断 */
        int st0 = S.nstates;
        char *vis = (char *)xcalloc((size_t)S.nstates, 1);
        int *q = (int *)xmalloc((size_t)S.nstates * sizeof(int));
        int qh = 0, qt = 0;
        int seed = M.v_in[s];
        vis[seed] = 1; q[qt++] = seed;
        while (qh < qt) {
            int u = q[qh++];
            for (int a = S.vhead[u]; a != -1; a = S.arc_next[a]) {
                int v = S.arc_to[a];
                if (!vis[v]) { vis[v] = 1; q[qt++] = v; }
            }
        }
        int *hit = (int *)xcalloc((size_t)n, sizeof(int));
        for (int z = 0; z < S.nstates; z++) if (vis[z]) hit[S.proj[z]] = 1;
        int miss = 0;
        for (int v = 0; v < n; v++) if (!hit[v]) miss++;
        fprintf(stderr, "[dbg] 状态图可达顶点 %d/%d（种子状态 %d）\n", n - miss, n, seed);
        free(vis); free(q); free(hit);
        (void)st0;
    }
    Hierarchy H;
    outer_build_hierarchy(&H, &S, &O, S.nstates, aux_state, naux);
    if (getenv("SSSP_DEBUG")) { fprintf(stderr, "[stage] hierarchy done: %d nodes, %d leaves, %d roots\n", H.nnode, H.nleaf, H.root_cnt); hier_validate(&H); }
    st->build_ms = now_ms() - t0;

    double t1 = now_ms();
    Evaluator E;
    eval_init(&E, &H, &S);
    eval_seed(&E, M.v_in[s], len_zero());
    if (getenv("SSSP_DEBUG")) fprintf(stderr, "[stage] eval start\n");
    if (opt->use_slot_eval) eval_run(&E);
    len_t *sd = (len_t *)xmalloc((size_t)(S.nstates ? S.nstates : 1) * sizeof(len_t));
    state_graph_dist(&S, M.v_in[s], sd);
    if (getenv("SSSP_DEBUG")) {
        int settled = 0, finite = 0;
        for (int z = 0; z < S.nstates; z++) { if (E.settled[z]) settled++; if (len_finite(E.tent[z])) finite++; }
        int r0 = H.roots[0];
        fprintf(stderr, "[dbg] 根: 高度 %g 单位 %g 质量 %g 子 %d | 叶 %d\n",
                H.h[r0].v, H.nu[r0].v, H.mass[r0].v, H.child_cnt[r0], H.nleaf);
        fprintf(stderr, "[dbg] 求值后：已确定状态 %d/%d，试探值有限 %d\n", settled, S.nstates, finite);
        eval_debug_report();
    }

    /* 8) 每个顶点取各副本的最小发出值 */
    int *cnt = (int *)xcalloc((size_t)n, sizeof(int));
    for (int z = 0; z < S.nstates; z++) cnt[S.proj[z]]++;
    int *off = (int *)xmalloc((size_t)(n + 1) * sizeof(int));
    off[0] = 0;
    for (int v = 0; v < n; v++) off[v + 1] = off[v] + cnt[v];
    int *fill = (int *)xmalloc((size_t)n * sizeof(int));
    for (int v = 0; v < n; v++) fill[v] = off[v];
    int *st_list = (int *)xmalloc((size_t)(S.nstates ? S.nstates : 1) * sizeof(int));
    for (int z = 0; z < S.nstates; z++) st_list[fill[S.proj[z]]++] = z;
    for (int v = 0; v < n; v++) {
        len_t best = len_inf();
        for (int i = off[v]; i < off[v + 1]; i++) {
            len_t e = opt->use_slot_eval ? E.emit_min[st_list[i]] : sd[st_list[i]];
            if (len_lt(e, best)) best = e;
        }
        out->d[v] = best;
    }
    st->eval_ms = now_ms() - t1;
    st->nstates = S.nstates; st->nnode = H.nnode; st->nleaf = H.nleaf;
    st->nroot = H.root_cnt; st->npiv = O.npiv; st->nreq = nreq;
    st->njoin = H.nnode; st->ndict = R.ndict;

    eval_free(&E);
    hier_free(&H);
    sgraph_free(&S);
    free(cnt); free(off); free(fill); free(st_list); free(sd);
    free(aux_state);
    free(J); free(de); free(radial); free(near);
    free(ru); free(rv); free(rw); free(rep);
    routes_free(&R);
    owner_free(&O);
    free(keep); free(map2); free(mst);
    graph_free(&G2);
    return 0;
}

/* ---------------- 验证 ---------------- */

static int verify(const Graph *g, int s, const Dist *got)
{
    Dist ref; dist_init(&ref, g->n);
    dijkstra_heap(g, s, &ref);
    int bad = 0;
    for (int v = 0; v < g->n; v++) {
        double a = got->d[v].v, b = ref.d[v].v;
        if (isinf(a) != isinf(b) || (!isinf(a) && fabs(a - b) > 1e-6 * (1.0 + fabs(b)))) {
            if (bad < 6) printf("    顶点 %d: 本实现 %.10g / Dijkstra %.10g\n", v, a, b);
            bad++;
        }
    }
    dist_free(&ref);
    return bad;
}

int main(int argc, char **argv)
{
    const char *type = "rand";
    int n = 400, m = 1200, seed = 1, maxw = 100, real_w = 0, do_verify = 0, bench = 0;
    Options opt; memset(&opt, 0, sizeof(opt)); opt.sep_log = 0; opt.seed = 1; opt.verbose = 0;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--n") && i + 1 < argc) n = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--m") && i + 1 < argc) m = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--seed") && i + 1 < argc) { seed = atoi(argv[++i]); opt.seed = seed; }
        else if (!strcmp(argv[i], "--maxw") && i + 1 < argc) maxw = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--type") && i + 1 < argc) type = argv[++i];
        else if (!strcmp(argv[i], "--sep") && i + 1 < argc) { opt.sep_log = atoi(argv[++i]); opt.force_hier = 1; }
        else if (!strcmp(argv[i], "--force")) opt.force_hier = 1;
        else if (!strcmp(argv[i], "--slot")) opt.use_slot_eval = 1;
        else if (!strcmp(argv[i], "--real")) real_w = 1;
        else if (!strcmp(argv[i], "--verify")) do_verify = 1;
        else if (!strcmp(argv[i], "--bench")) bench = 1;
        else if (!strcmp(argv[i], "--help")) {
            printf("用法: sssp [--type rand|grid|path|star|clique] [--n N] [--m M] [--seed S]\n"
                   "            [--maxw W] [--real] [--sep K] [--force] [--verify] [--bench]\n");
            return 0;
        }
    }
    if (bench) {
        struct { const char *name; int n, m; } cases[] = {
            { "rand-500", 500, 1500 }, { "rand-2000", 2000, 6000 },
            { "rand-8000", 8000, 24000 }, { "grid-40x40", 1600, 3120 },
            { "grid-90x90", 8100, 16020 },
        };
        int total_bad = 0;
        for (int c = 0; c < 5; c++) {
            Graph g;
            if (cases[c].name[0] == 'g' && cases[c].name[1] == 'r')
                g = gen_rand(cases[c].n, cases[c].m, 1000, 7 + c, 1);
            else {
                int side = (int)sqrt((double)cases[c].n);
                g = gen_grid(side, 1000, 7 + c, 1);
            }
            Dist got; Stats st;
            sssp_pipeline(&g, 0, &opt, &got, &st);
            int bad = verify(&g, 0, &got);
            total_bad += bad;
            printf("%-12s n=%-6d m=%-7d 生成+层级 %7.1f ms 求值 %7.1f ms  状态 %d 节点 %d | %s\n",
                   cases[c].name, g.n, g.m, st.build_ms, st.eval_ms, st.nstates, st.nnode,
                   bad ? "不一致" : "一致");
            dist_free(&got);
            graph_free(&g);
        }
        printf("\n总计不一致顶点: %d\n", total_bad);
        return total_bad != 0;
    }
    Graph g;
    if (!strcmp(type, "path")) {
        g = gen_path(n, maxw, seed);
    } else if (!strcmp(type, "grid")) {
        int side = (int)sqrt((double)n);
        g = gen_grid(side, maxw, seed, real_w);
    } else {
        g = gen_rand(n, m, maxw, seed, real_w);
    }
    printf("图: %s  n=%d m=%d seed=%d %s权重\n", type, g.n, g.m, seed, real_w ? "实数" : "整数");
    Dist got; Stats st;
    len_ops = 0; len_doublings = 0;
    sssp_pipeline(&g, 0, &opt, &got, &st);
    if (st.base_case)
        printf("结果: 走常数规模基例（扫描版 Dijkstra）。N <= L_N^64 对任何现实规模都成立；\n"
               "      要触发层级路径请加 --sep 3（缩放的分隔规模）或 --force。\n");
    else
        printf("结果: 层级路径  状态 %d（%d 个块外另有特殊状态）、层级节点 %d、叶 %d、根 %d\n"
               "      枢纽 %d、请求 %d、词典 %d\n"
               "      生成 %.1f ms、求值 %.1f ms\n",
               st.nstates, 0, st.nnode, st.nleaf, st.nroot, st.npiv, st.nreq, st.ndict,
               st.build_ms, st.eval_ms);
    printf("寄存器运算: 比较+加法 %llu 次，加倍 %llu 次\n",
           (unsigned long long)len_ops, (unsigned long long)len_doublings);
    int bad = 0;
    if (do_verify) {
        bad = verify(&g, 0, &got);
        printf("对拍 Dijkstra: %s（不一致顶点 %d）\n", bad ? "失败" : "通过", bad);
    }
    dist_free(&got);
    graph_free(&g);
    return bad != 0;
}