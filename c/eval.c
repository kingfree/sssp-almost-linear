/* eval.c -- 第 4 节：有序层级的求值（算法 7 Scan、算法 8 Evaluate）
 *
 * 数据结构：叶布局（DFS 序）+ 段树维护"未开启区间的最小试探值"（对应引理 3.9 的
 * split--findmin 接口；这里用 O(log N) 的段树实现同一组操作：区间最小、减键、
 * 取最小叶）。每个内部节点维护自己的 piece 双向链表（未注册连续段 RUN 与已注册
 * 子节点 CHILD），以及时钟（槽下标）、到期时刻与游标。
 */
#include "sssp.h"

static void seg_pull(Evaluator *E, int i);

static void seg_init(Evaluator *E, int nleaf)
{
    int sz = 1;
    while (sz < nleaf) sz <<= 1;
    E->seg_n = sz;
    E->seg_arg = (int *)xmalloc((size_t)(2 * sz) * sizeof(int));
    E->seg_val = (len_t *)xmalloc((size_t)(2 * sz) * sizeof(len_t));
    for (int i = 0; i < 2 * sz; i++) { E->seg_arg[i] = -1; E->seg_val[i] = len_inf(); }
}
static void seg_set(Evaluator *E, int leafpos, len_t v)
{
    int i = leafpos + E->seg_n;
    E->seg_val[i] = v; E->seg_arg[i] = (len_finite(v) ? leafpos : -1);
    for (i >>= 1; i >= 1; i >>= 1) seg_pull(E, i);
    if (i == 0) { }
}
static void seg_pull(Evaluator *E, int i)
{
    int l = 2 * i, r = l + 1;
    len_t vl = E->seg_val[l], vr = E->seg_val[r];
    int c = len_cmp(vl, vr);
    if (c < 0 || (c == 0 && E->seg_arg[l] >= 0 && (E->seg_arg[r] < 0 || E->seg_arg[l] <= E->seg_arg[r]))) {
        E->seg_val[i] = vl; E->seg_arg[i] = E->seg_arg[l];
    } else {
        E->seg_val[i] = vr; E->seg_arg[i] = E->seg_arg[r];
    }
}
/* 区间 [lo,hi] 内的最小叶 */
static int seg_query(Evaluator *E, int lo, int hi, len_t *out)
{
    int res = -1;
    len_t best = len_inf();
    int l = lo + E->seg_n, r = hi + E->seg_n + 1;
    while (l < r) {
        if (l & 1) {
            if (len_lt(E->seg_val[l], best) || (len_eq(E->seg_val[l], best) && E->seg_arg[l] >= 0 && (res < 0 || E->seg_arg[l] < res))) {
                best = E->seg_val[l]; res = E->seg_arg[l];
            }
            l++;
        }
        if (r & 1) {
            r--;
            if (len_lt(E->seg_val[r], best) || (len_eq(E->seg_val[r], best) && E->seg_arg[r] >= 0 && (res < 0 || E->seg_arg[r] < res))) {
                best = E->seg_val[r]; res = E->seg_arg[r];
            }
        }
        l >>= 1; r >>= 1;
    }
    if (out) *out = best;
    return res;
}

/* ---------------- piece 池 ---------------- */

static int pc_new(Evaluator *E)
{
    if (E->npiece == E->cap_piece) {
        E->cap_piece = E->cap_piece ? E->cap_piece * 2 : 64;
        E->pc_kind = (int *)xrealloc(E->pc_kind, (size_t)E->cap_piece * sizeof(int));
        E->pc_lo = (int *)xrealloc(E->pc_lo, (size_t)E->cap_piece * sizeof(int));
        E->pc_hi = (int *)xrealloc(E->pc_hi, (size_t)E->cap_piece * sizeof(int));
        E->pc_child = (int *)xrealloc(E->pc_child, (size_t)E->cap_piece * sizeof(int));
        E->pc_next = (int *)xrealloc(E->pc_next, (size_t)E->cap_piece * sizeof(int));
        E->pc_prev = (int *)xrealloc(E->pc_prev, (size_t)E->cap_piece * sizeof(int));
        E->pc_owner = (int *)xrealloc(E->pc_owner, (size_t)E->cap_piece * sizeof(int));
    }
    int i = E->npiece++;
    E->pc_kind[i] = 0; E->pc_lo[i] = E->pc_hi[i] = -1; E->pc_child[i] = -1;
    E->pc_next[i] = E->pc_prev[i] = -1; E->pc_owner[i] = -1;
    return i;
}

void eval_init(Evaluator *E, Hierarchy *H, const SGraph *S)
{
    memset(E, 0, sizeof(*E));
    E->H = H; E->S = S;
    int ns = S->nstates;
    E->tent = (len_t *)xmalloc((size_t)(ns ? ns : 1) * sizeof(len_t));
    E->settled = (int *)xcalloc((size_t)(ns ? ns : 1), sizeof(int));
    E->retired = (int *)xcalloc((size_t)(ns ? ns : 1), sizeof(int));
    E->emit_min = (len_t *)xmalloc((size_t)(ns ? ns : 1) * sizeof(len_t));
    E->leaf_pos_of_state = (int *)xmalloc((size_t)(ns ? ns : 1) * sizeof(int));
    for (int i = 0; i < ns; i++) { E->tent[i] = len_inf(); E->emit_min[i] = len_inf(); E->leaf_pos_of_state[i] = -1; }
    for (int i = 0; i < H->nleaf; i++) E->leaf_pos_of_state[H->leaf_state[i]] = i;
    seg_init(E, H->nleaf);
    int NN = H->nnode;
    E->opened = (int *)xcalloc((size_t)(NN ? NN : 1), sizeof(int));
    E->clock_origin = (len_t *)xmalloc((size_t)(NN ? NN : 1) * sizeof(len_t));
    E->slot_k = (int *)xcalloc((size_t)(NN ? NN : 1), sizeof(int));
    E->expiry = (len_t *)xmalloc((size_t)(NN ? NN : 1) * sizeof(len_t));
    E->piece_head = (int *)xmalloc((size_t)(NN ? NN : 1) * sizeof(int));
    E->cursor_piece = (int *)xmalloc((size_t)(NN ? NN : 1) * sizeof(int));
    for (int i = 0; i < NN; i++) { E->clock_origin[i] = len_zero(); E->expiry[i] = len_inf(); E->piece_head[i] = -1; E->cursor_piece[i] = -1; }
    /* 每个内部节点初始一个 RUN 覆盖其叶区间 */
    for (int i = 0; i < NN; i++) {
        if (H->is_leaf[i]) continue;
        int pc = pc_new(E);
        E->pc_kind[pc] = 0;
        E->pc_lo[pc] = H->first_leaf[i];
        E->pc_hi[pc] = H->last_leaf[i];
        E->pc_owner[pc] = i;
        E->piece_head[i] = pc;
        E->cursor_piece[i] = pc;
    }
}

void eval_free(Evaluator *E)
{
    free(E->tent); free(E->settled); free(E->retired); free(E->emit_min);
    free(E->leaf_pos_of_state); free(E->seg_arg); free(E->seg_val);
    free(E->opened); free(E->clock_origin); free(E->slot_k); free(E->expiry);
    free(E->piece_head); free(E->cursor_piece);
    free(E->pc_kind); free(E->pc_lo); free(E->pc_hi); free(E->pc_child);
    free(E->pc_next); free(E->pc_prev); free(E->pc_owner);
    memset(E, 0, sizeof(*E));
}

void eval_seed(Evaluator *E, int state, len_t cost)
{
    if (state < 0) return;
    if (len_lt(cost, E->tent[state])) {
        E->tent[state] = cost;
        int lp = E->leaf_pos_of_state[state];
        if (lp >= 0) seg_set(E, lp, cost);
    }
}

static void eval_relax(Evaluator *E, int state, len_t val)
{
    for (int a = E->S->vhead[state]; a != -1; a = E->S->arc_next[a]) {
        int v = E->S->arc_to[a];
        if (E->settled[v] || E->retired[v]) continue;
        len_t nd = len_add(val, E->S->arc_len[a]);
        if (len_lt(nd, E->tent[v])) {
            E->tent[v] = nd;
            int lp = E->leaf_pos_of_state[v];
            if (lp >= 0) seg_set(E, lp, nd);
        }
    }
}

/* 子节点 -> 其在父中的 piece 下标（按需查找） */
static int find_child_piece(Evaluator *E, int parent, int child)
{
    for (int p = E->piece_head[parent]; p != -1; p = E->pc_next[p])
        if (E->pc_kind[p] == 1 && E->pc_child[p] == child) return p;
    return -1;
}

static void unlink_piece(Evaluator *E, int parent, int pc)
{
    if (E->pc_prev[pc] >= 0) E->pc_next[E->pc_prev[pc]] = E->pc_next[pc];
    else E->piece_head[parent] = E->pc_next[pc];
    if (E->pc_next[pc] >= 0) E->pc_prev[E->pc_next[pc]] = E->pc_prev[pc];
    E->pc_next[pc] = E->pc_prev[pc] = -2;   /* 已摘除 */
}

/* 在父节点里把覆盖 leafpos 的 RUN 拆成 [左 RUN, CHILD, 右 RUN] */
static int split_run(Evaluator *E, int parent, int pc, int child)
{
    int H0 = E->H->first_leaf[child], H1 = E->H->last_leaf[child];
    int lo = E->pc_lo[pc], hi = E->pc_hi[pc];
    int left = -1, right = -1;
    int mid = pc_new(E);
    E->pc_kind[mid] = 1; E->pc_child[mid] = child; E->pc_owner[mid] = parent;
    if (lo < H0) {
        left = pc_new(E);
        E->pc_kind[left] = 0; E->pc_lo[left] = lo; E->pc_hi[left] = H0 - 1; E->pc_owner[left] = parent;
    }
    if (H1 < hi) {
        right = pc_new(E);
        E->pc_kind[right] = 0; E->pc_lo[right] = H1 + 1; E->pc_hi[right] = hi; E->pc_owner[right] = parent;
    }
    int prev = E->pc_prev[pc], next = E->pc_next[pc];
    /* 组装 */
    int first = (left >= 0) ? left : mid;
    int lastp = (right >= 0) ? right : mid;
    if (left >= 0) { E->pc_next[left] = mid; E->pc_prev[mid] = left; }
    if (right >= 0) { E->pc_next[mid] = right; E->pc_prev[right] = mid; }
    E->pc_prev[first] = prev;
    if (prev >= 0) E->pc_next[prev] = first; else E->piece_head[parent] = first;
    E->pc_next[lastp] = next;
    if (next >= 0) E->pc_prev[next] = lastp;
    E->pc_next[pc] = E->pc_prev[pc] = -2;
    return first;
}

static int child_of_leafpos(Hierarchy *H, int node, int leafpos)
{
    for (int k = 0; k < H->child_cnt[node]; k++) {
        int c = H->children[H->child_start[node] + k];
        if (H->first_leaf[c] <= leafpos && leafpos <= H->last_leaf[c]) return c;
    }
    return -1;
}

static void visit_node(Evaluator *E, int X, len_t lo, int has_lo, len_t hi, int has_hi, len_t cap);

/* 诊断计数 */
static long dbg_split, dbg_leafvisit, dbg_emit, dbg_childvisit, dbg_slots, dbg_rejectwin;
static long dbg_visitnode, dbg_immediate_break, dbg_open;

/* 扫描一个节点的 piece 链表：处理一个高度为 nu 的时间槽片段 */
static void scan_slot(Evaluator *E, int X, len_t slo, len_t shi, len_t cap, len_t pad)
{
    Hierarchy *H = E->H;
    int p = E->piece_head[X];
    while (p != -1) {
        int nx = E->pc_next[p];
        if (E->pc_kind[p] == 0) {
            len_t z; 
            int lp = seg_query(E, E->pc_lo[p], E->pc_hi[p], &z);
            if (lp >= 0 && len_finite(z) && len_ge(z, slo) && len_lt(z, shi) && len_le(z, cap)) {
                int c = child_of_leafpos(H, X, lp);
                if (c < 0) { p = nx; continue; }
                int firstpc = split_run(E, X, p, c);
                dbg_split++;
                /* 双侧子段的补充候选值：在相邻子节点代表点上发布 u = max{b, z + S(i) + g(i,j)} */
                int cpos = H->child_rank[c];
                /* 左邻居 */
                if (cpos > 0) {
                    int cj = H->children[H->child_start[X] + cpos - 1];
                    int s1 = H->leaf_state[H->first_leaf[cj]];
                    len_t g = H->gap[H->child_start[X] + cpos];
                    len_t u = len_add(len_add(z, H->mass[c]), g);
                    if (len_lt(u, pad)) u = pad;       /* 填充到名义槽右端点 b：本轮不再到期 */
                    if (len_lt(u, E->tent[s1])) {
                        /* 只有当它在 cj 的许可窗口内且不超过其期限时才发布 */
                        E->tent[s1] = u;
                        int lp2 = E->leaf_pos_of_state[s1];
                        if (lp2 >= 0) seg_set(E, lp2, u);
                    }
                }
                /* 右邻居 */
                if (cpos + 1 < H->child_cnt[X]) {
                    int cj = H->children[H->child_start[X] + cpos + 1];
                    int s1 = H->leaf_state[H->first_leaf[cj]];
                    len_t g = H->gap[H->child_start[X] + cpos + 1];
                    len_t u = len_add(len_add(z, H->mass[c]), g);
                    if (len_lt(u, pad)) u = pad;
                    if (len_lt(u, E->tent[s1])) {
                        E->tent[s1] = u;
                        int lp2 = E->leaf_pos_of_state[s1];
                        if (lp2 >= 0) seg_set(E, lp2, u);
                    }
                }
                p = firstpc;                        /* 游标回到第一个新片段 */
                nx = E->pc_next[p];
                continue;
            }
            /* 最小值为无穷：本槽内没有到期项。可能是尚未收到候选值，
               因此不能摘除（后续松弛仍可能让它变小）；到节点过期后自然停止。 */
            p = nx;
        } else {
            int c = E->pc_child[p];
            dbg_childvisit++;
            visit_node(E, c, slo, 1, shi, 1, len_min(cap, E->expiry[c]));
            /* 子节点若已无 piece 则摘掉 */
            if (E->piece_head[c] == -1) unlink_piece(E, X, p);
            p = nx;
        }
    }
}

static void open_node(Evaluator *E, int X, len_t D0)
{
    E->opened[X] = 1;
    E->clock_origin[X] = D0;
    E->slot_k[X] = 0;
    E->expiry[X] = len_add(D0, E->H->mass[X]);
    E->cursor_piece[X] = E->piece_head[X];
}

static void visit_node(Evaluator *E, int X, len_t lo, int has_lo, len_t hi, int has_hi, len_t cap)
{
    Hierarchy *H = E->H;
    if (H->is_leaf[X]) {
        int st = H->state_of_leaf[X];
        if (E->settled[st] || E->retired[st]) return;
        len_t v = E->tent[st];
        if (!len_finite(v)) return;
        dbg_leafvisit++;
        if (has_lo && len_lt(v, lo)) { dbg_rejectwin++; return; }
        if (has_hi && len_ge(v, hi)) { dbg_rejectwin++; return; }
        if (len_gt(v, cap)) { dbg_rejectwin++; return; }
        /* 确定并松弛 */
        E->settled[st] = 1;
        dbg_emit++;
        if (len_lt(v, E->emit_min[st])) E->emit_min[st] = v;
        int lp = E->leaf_pos_of_state[st];
        if (lp >= 0) seg_set(E, lp, len_inf());
        eval_relax(E, st, v);
        return;
    }
    dbg_visitnode++;
    if (!E->opened[X]) {
        len_t z;
        int lp = seg_query(E, H->first_leaf[X], H->last_leaf[X], &z);
        if (lp < 0 || !len_finite(z)) return;
        dbg_open++;
        open_node(E, X, z);
    }
    len_t nu = H->nu[X];
    if (!len_gt(nu, len_zero())) nu = len_make(1.0);
    for (;;) {
        len_t a = E->clock_origin[X];
        if (E->slot_k[X] > 0) a = len_mul_pow2(nu, 0);   /* 占位，下面重算 */
        /* a = origin + slot_k * nu（用倍增表示 k*nu） */
        {
            len_t off = len_zero();
            int k = E->slot_k[X];
            for (int b = 0; (1 << b) <= k; b++) {
                if (k & (1 << b)) off = len_add(off, len_mul_pow2(nu, b));
            }
            a = len_add(E->clock_origin[X], off);
        }
        len_t b = len_add(a, nu);
        if (has_hi && len_ge(a, hi)) { dbg_immediate_break++; break; }
        if (len_gt(a, E->expiry[X])) break;          /* 节点到期：覆盖到 D0 + S(C)（含端点） */
        if (has_lo && len_le(b, lo)) { E->slot_k[X]++; continue; }
        if (len_gt(a, cap)) break;
        if (E->slot_k[X] > 100000000) break;          /* 安全阀 */
        len_t slo = (has_lo && len_gt(lo, a)) ? lo : a;
        len_t shi = (has_hi && len_lt(hi, b)) ? hi : b;
        len_t cap2 = len_min(len_min(b, cap), E->expiry[X]);
        dbg_slots++;
        scan_slot(E, X, slo, shi, cap2, b);
        if (E->piece_head[X] == -1) break;
        E->slot_k[X]++;
    }
}

/* 处理一个根（叶子或内部节点） */
static void eval_root(Evaluator *E, int X)
{
    Hierarchy *H = E->H;
    if (H->is_leaf[X]) {
        int st = H->state_of_leaf[X];
        if (!E->settled[st] && len_finite(E->tent[st])) {
            E->settled[st] = 1;
            if (len_lt(E->tent[st], E->emit_min[st])) E->emit_min[st] = E->tent[st];
            int lp = E->leaf_pos_of_state[st];
            if (lp >= 0) seg_set(E, lp, len_inf());
            eval_relax(E, st, E->tent[st]);
        }
        return;
    }
    len_t z;
    int lp = seg_query(E, H->first_leaf[X], H->last_leaf[X], &z);
    if (lp < 0 || !len_finite(z)) return;
    if (!E->opened[X]) open_node(E, X, z);
    visit_node(E, X, len_zero(), 0, len_zero(), 0, len_inf());
}

/* 引理 4.6 的求值：按根的凝聚序处理；这里再加一层"直到没有新进展"的
   外循环，使后到的候选值也能推进（等价于把根序中的所有依赖都处理完）。*/
void eval_run(Evaluator *E)
{
    Hierarchy *H = E->H;
    if (getenv("SSSP_DEBUG")) {
        dbg_split = dbg_leafvisit = dbg_emit = dbg_childvisit = dbg_slots = dbg_rejectwin = 0;
    }
    int passes = 0, maxpass = H->root_cnt + 8;
    for (;;) {
        int progress = 0;
        for (int r = 0; r < H->root_cnt; r++) {
            int before = 0;
            if (getenv("SSSP_DEBUG")) for (int z = 0; z < E->S->nstates; z++) if (E->settled[z]) before++;
            eval_root(E, H->roots[r]);
            if (getenv("SSSP_DEBUG")) {
                int after = 0;
                for (int z = 0; z < E->S->nstates; z++) if (E->settled[z]) after++;
                if (after > before) progress = 1;
            } else progress = 1;
        }
        if (!progress || ++passes >= maxpass) break;
    }
}

void eval_debug_report(void)
{
    fprintf(stderr, "[dbg] 槽扫描 %ld，拆段 %ld，子访问 %ld，叶访问 %ld，确定 %ld，窗口拒绝 %ld\n",
            dbg_slots, dbg_split, dbg_childvisit, dbg_leafvisit, dbg_emit, dbg_rejectwin);
    fprintf(stderr, "[dbg] visit_node %ld，开启 %ld，立即跳出 %ld\n", dbg_visitnode, dbg_open, dbg_immediate_break);
}

len_t eval_distance(const Evaluator *E, int vertex, const int *states_of_vertex, int cnt)
{
    (void)vertex;
    len_t best = len_inf();
    for (int i = 0; i < cnt; i++) {
        int st = states_of_vertex[i];
        if (len_lt(E->emit_min[st], best)) best = E->emit_min[st];
    }
    return best;
}
/* 状态图上的精确 Dijkstra：独立验证 §8.1 的构造确实以原长度实现了每条物理测地线
   （引理 8.1 / 5.5 的结论）。 */
void state_graph_dist(const SGraph *S, int seed_state, len_t *out_state_dist)
{
    int n = S->nstates;
    for (int i = 0; i < n; i++) out_state_dist[i] = len_inf();
    Heap h; heap_init(&h, n);
    out_state_dist[seed_state] = len_zero();
    heap_push(&h, len_zero(), seed_state);
    char *done = (char *)xcalloc((size_t)(n ? n : 1), 1);
    for (;;) {
        len_t du; int u = heap_pop(&h, &du);
        if (u < 0) break;
        if (done[u]) continue;
        done[u] = 1;
        for (int a = S->vhead[u]; a != -1; a = S->arc_next[a]) {
            int v = S->arc_to[a];
            len_t nd = len_add(du, S->arc_len[a]);
            if (len_lt(nd, out_state_dist[v])) { out_state_dist[v] = nd; heap_decrease(&h, v, nd); }
        }
    }
    free(done);
    heap_free(&h);
}
