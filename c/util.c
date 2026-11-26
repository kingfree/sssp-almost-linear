/* util.c -- 基础容器、随机数、并查集、二叉堆 */
#include "sssp.h"

uint64_t len_ops = 0;
uint64_t len_doublings = 0;

void *xmalloc(size_t n)
{
    void *p = malloc(n ? n : 1);
    if (!p) { fprintf(stderr, "oom\n"); exit(1); }
    return p;
}
void *xcalloc(size_t n, size_t s)
{
    void *p = calloc(n ? n : 1, s ? s : 1);
    if (!p) { fprintf(stderr, "oom\n"); exit(1); }
    return p;
}
void *xrealloc(void *p, size_t n)
{
    void *q = realloc(p, n ? n : 1);
    if (!q) { fprintf(stderr, "oom\n"); exit(1); }
    return q;
}

void ivec_init(ivec *v) { v->a = NULL; v->n = 0; v->cap = 0; }
void ivec_reserve(ivec *v, int cap)
{
    if (cap <= v->cap) return;
    v->cap = cap < 8 ? 8 : cap;
    v->a = (int *)xrealloc(v->a, (size_t)v->cap * sizeof(int));
}
void ivec_push(ivec *v, int x)
{
    if (v->n == v->cap) ivec_reserve(v, v->cap ? v->cap * 2 : 8);
    v->a[v->n++] = x;
}
void ivec_clear(ivec *v) { v->n = 0; }
void ivec_free(ivec *v) { free(v->a); v->a = NULL; v->n = v->cap = 0; }

void lvec_init(lvec *v) { v->a = NULL; v->n = 0; v->cap = 0; }
void lvec_push(lvec *v, len_t x)
{
    if (v->n == v->cap) {
        v->cap = v->cap ? v->cap * 2 : 8;
        v->a = (len_t *)xrealloc(v->a, (size_t)v->cap * sizeof(len_t));
    }
    v->a[v->n++] = x;
}
void lvec_free(lvec *v) { free(v->a); v->a = NULL; v->n = v->cap = 0; }

/* splitmix64 + xorshift */
static uint64_t rng_s = 0x2545F4914F6CDD1DULL;
void rng_seed(uint64_t s) { rng_s = s ? s : 0x9E3779B97F4A7C15ULL; }
uint64_t rng_next(void)
{
    uint64_t z = (rng_s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
int rng_int(int n)
{
    if (n <= 1) return 0;
    return (int)(rng_next() % (uint64_t)n);
}
double rng_double(void) { return (double)(rng_next() >> 11) / 9007199254740992.0; }

void dsu_init(DSU *d, int n)
{
    d->n = n;
    d->p = (int *)xmalloc((size_t)n * sizeof(int));
    d->rk = (int *)xcalloc((size_t)n, sizeof(int));
    for (int i = 0; i < n; i++) d->p[i] = i;
}
int dsu_find(DSU *d, int x)
{
    int r = x;
    while (d->p[r] != r) r = d->p[r];
    while (d->p[x] != r) { int nx = d->p[x]; d->p[x] = r; x = nx; }
    return r;
}
int dsu_unite(DSU *d, int a, int b)
{
    int ra = dsu_find(d, a), rb = dsu_find(d, b);
    if (ra == rb) return 0;
    if (d->rk[ra] < d->rk[rb]) { int t = ra; ra = rb; rb = t; }
    d->p[rb] = ra;
    if (d->rk[ra] == d->rk[rb]) d->rk[ra]++;
    return 1;
}
void dsu_free(DSU *d) { free(d->p); free(d->rk); d->p = d->rk = NULL; d->n = 0; }

void heap_init(Heap *h, int maxn)
{
    h->cap = maxn > 0 ? maxn : 1;
    h->key = (len_t *)xmalloc((size_t)h->cap * sizeof(len_t));
    h->item = (int *)xmalloc((size_t)h->cap * sizeof(int));
    h->pos = (int *)xmalloc((size_t)h->cap * sizeof(int));
    for (int i = 0; i < h->cap; i++) h->pos[i] = -1;
    h->n = 0;
}
void heap_free(Heap *h) { free(h->key); free(h->item); free(h->pos); h->key = NULL; h->item = NULL; h->pos = NULL; h->n = h->cap = 0; }
void heap_clear(Heap *h)
{
    for (int i = 0; i < h->n; i++) h->pos[h->item[i]] = -1;
    h->n = 0;
}
int heap_empty(Heap *h) { return h->n == 0; }

static void heap_up(Heap *h, int i)
{
    while (i > 0) {
        int p = (i - 1) / 2;
        if (len_lt(h->key[i], h->key[p])) {
            len_t tk = h->key[i]; h->key[i] = h->key[p]; h->key[p] = tk;
            int ti = h->item[i]; h->item[i] = h->item[p]; h->item[p] = ti;
            h->pos[h->item[i]] = i; h->pos[h->item[p]] = p;
            i = p;
        } else break;
    }
}
static void heap_down(Heap *h, int i)
{
    for (;;) {
        int l = 2 * i + 1, r = l + 1, m = i;
        if (l < h->n && len_lt(h->key[l], h->key[m])) m = l;
        if (r < h->n && len_lt(h->key[r], h->key[m])) m = r;
        if (m == i) break;
        len_t tk = h->key[i]; h->key[i] = h->key[m]; h->key[m] = tk;
        int ti = h->item[i]; h->item[i] = h->item[m]; h->item[m] = ti;
        h->pos[h->item[i]] = i; h->pos[h->item[m]] = m;
        i = m;
    }
}
int heap_push(Heap *h, len_t k, int item)
{
    if (h->pos[item] >= 0) { heap_decrease(h, item, k); return h->pos[item]; }
    if (h->n == h->cap) {
        int old = h->cap;
        h->cap = h->cap * 2 + 1;
        h->key = (len_t *)xrealloc(h->key, (size_t)h->cap * sizeof(len_t));
        h->item = (int *)xrealloc(h->item, (size_t)h->cap * sizeof(int));
        h->pos = (int *)xrealloc(h->pos, (size_t)h->cap * sizeof(int));
        for (int i = old; i < h->cap; i++) h->pos[i] = -1;
    }
    int i = h->n++;
    h->key[i] = k; h->item[i] = item; h->pos[item] = i;
    heap_up(h, i);
    return i;
}
void heap_decrease(Heap *h, int item, len_t k)
{
    int i = h->pos[item];
    if (i < 0) { heap_push(h, k, item); return; }
    if (!len_lt(k, h->key[i])) return;
    h->key[i] = k;
    heap_up(h, i);
}
int heap_pop(Heap *h, len_t *out)
{
    if (h->n == 0) return -1;
    int it = h->item[0];
    if (out) *out = h->key[0];
    h->pos[it] = -1;
    h->n--;
    if (h->n > 0) {
        h->key[0] = h->key[h->n]; h->item[0] = h->item[h->n];
        h->pos[h->item[0]] = 0;
        heap_down(h, 0);
    }
    return it;
}
