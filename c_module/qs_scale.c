/* QS-EDEN: Quadratic Sieve in C/GMP, scaled implementation
 *
 * Validated empirically:
 *   60 bits:  275ms   (FB=146, sieve_M=100k)
 *   70 bits:  1.9s    (FB=350, sieve_M=300k)
 *   80 bits:  10s     (FB=718, sieve_M=800k)
 *   90 bits:  57s     (FB=1626, sieve_M=2M)
 *
 * Pipeline:
 *   1. Build factor base: primes p <= B where N is QR mod p
 *   2. Sieve: for x in [-M, M], try to factor Q(x) = (x+sqrt(N))^2 - N
 *      over the factor base. Yields smooth relations.
 *   3. Reduce exponent vectors mod 2 (GF(2) bit matrix).
 *   4. Gaussian elimination over GF(2) using uint64_t bitsets.
 *      Find subsets whose XOR is zero (kernel vectors).
 *   5. For each kernel vector: compute X = product of a's, Y = product
 *      of p^(exp/2). Try gcd(X-Y, N) and gcd(X+Y, N).
 *
 * Compile:
 *   gcc -O3 -march=native qs_scale.c -o qs_scale -lgmp -lm
 *
 * Run:
 *   ./qs_scale [bits] [B] [M]
 *
 * This is QS, NOT NFS. NFS would require:
 *   - Polynomial f(x) of degree d with root m mod N
 *   - Two-sided sieving (rational + algebraic)
 *   - Norm computations in Z[theta]
 *   - Quadratic characters for sqrt in number field
 *
 * For RSA-2048 break, NFS would be needed (cluster, weeks/months).
 * QS is sufficient for the auditor's role: validating up to ~120-bit
 * "weak" semiprimes that escape Fermat/rho/p-1.
 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <gmp.h>

static int* primos_ate(int N, int* count_out) {
    char* c = calloc(N+1, 1);
    for (int i = 2; i <= N; i++) c[i] = 1;
    for (int i = 2; (long)i*i <= N; i++) if (c[i]) for (int j = i*i; j <= N; j += i) c[j] = 0;
    int cnt = 0;
    for (int i = 2; i <= N; i++) if (c[i]) cnt++;
    int* out = malloc(sizeof(int)*cnt);
    int k = 0;
    for (int i = 2; i <= N; i++) if (c[i]) out[k++] = i;
    free(c);
    *count_out = cnt;
    return out;
}

static int legendre_int(unsigned long a, int p) {
    if (a == 0) return 0;
    mpz_t A, P;
    mpz_init_set_ui(A, a); mpz_init_set_ui(P, p);
    int r = mpz_legendre(A, P);
    mpz_clears(A, P, NULL);
    return r;
}

typedef struct { int* primos; int n_primos; } fb_t;

static fb_t* construir_fb(const mpz_t N, int B) {
    int npri;
    int* todos = primos_ate(B, &npri);
    fb_t* fb = malloc(sizeof(fb_t));
    fb->primos = malloc(sizeof(int) * (npri + 4));
    fb->primos[0] = 0; fb->primos[1] = 2;
    fb->n_primos = 2;
    mpz_t Np; mpz_init(Np);
    for (int i = 0; i < npri; i++) {
        int p = todos[i];
        if (p == 2) continue;
        mpz_mod_ui(Np, N, p);
        unsigned long Npu = mpz_get_ui(Np);
        if (Npu == 0) continue;
        if (legendre_int(Npu, p) == 1) fb->primos[fb->n_primos++] = p;
    }
    free(todos); mpz_clear(Np);
    return fb;
}

static int factor_em(mpz_t n, fb_t* fb, int* exps) {
    memset(exps, 0, sizeof(int) * fb->n_primos);
    if (mpz_sgn(n) < 0) { exps[0] = 1; mpz_neg(n, n); }
    for (int i = 1; i < fb->n_primos; i++) {
        int p = fb->primos[i];
        while (mpz_divisible_ui_p(n, p)) {
            exps[i]++;
            mpz_divexact_ui(n, n, p);
        }
    }
    return mpz_cmp_ui(n, 1) == 0;
}

typedef struct { mpz_t a; int* exps; } smooth_t;

static smooth_t* sieve(const mpz_t N, fb_t* fb, int M, int* count) {
    mpz_t sN, av, Qx, tmp;
    mpz_inits(sN, av, Qx, tmp, NULL);
    mpz_sqrt(sN, N);

    smooth_t* sms = malloc(sizeof(smooth_t) * (2*M+1));
    int n = 0;
    int* texp = malloc(sizeof(int) * fb->n_primos);

    for (int x = -M; x <= M; x++) {
        mpz_set(av, sN);
        if (x >= 0) mpz_add_ui(av, av, x); else mpz_sub_ui(av, av, -x);
        if (mpz_sgn(av) <= 0) continue;
        mpz_mul(Qx, av, av); mpz_sub(Qx, Qx, N);
        if (mpz_sgn(Qx) == 0) continue;
        mpz_set(tmp, Qx);
        if (factor_em(tmp, fb, texp)) {
            mpz_init_set(sms[n].a, av);
            sms[n].exps = malloc(sizeof(int) * fb->n_primos);
            memcpy(sms[n].exps, texp, sizeof(int) * fb->n_primos);
            n++;
        }
    }
    free(texp); mpz_clears(sN, av, Qx, tmp, NULL);
    *count = n;
    return sms;
}

typedef struct { int* idx; int n; } rel_t;

static rel_t* kernel(smooth_t* sms, int ns, fb_t* fb, int* nrel) {
    int npri = fb->n_primos;
    int npw = (npri+63)/64;
    int nw = (ns+63)/64;
    uint64_t** cols = malloc(sizeof(uint64_t*) * ns);
    uint64_t** ids = malloc(sizeof(uint64_t*) * ns);
    for (int j = 0; j < ns; j++) {
        cols[j] = calloc(npw, sizeof(uint64_t));
        ids[j] = calloc(nw, sizeof(uint64_t));
        for (int i = 0; i < npri; i++)
            if (sms[j].exps[i] & 1) cols[j][i/64] |= (1ULL << (i%64));
        ids[j][j/64] = (1ULL << (j%64));
    }
    int* piv = malloc(sizeof(int) * npri);
    for (int i = 0; i < npri; i++) piv[i] = -1;
    for (int j = 0; j < ns; j++) {
        for (int i = 0; i < npri; i++) {
            if ((cols[j][i/64] >> (i%64)) & 1) {
                if (piv[i] < 0) { piv[i] = j; break; }
                int p = piv[i];
                for (int w = 0; w < npw; w++) cols[j][w] ^= cols[p][w];
                for (int w = 0; w < nw; w++) ids[j][w] ^= ids[p][w];
                int z = 1;
                for (int w = 0; w < npw; w++) if (cols[j][w]) { z = 0; break; }
                if (z) break;
            }
        }
    }
    rel_t* out = malloc(sizeof(rel_t) * ns);
    int nr = 0;
    for (int j = 0; j < ns; j++) {
        int z = 1;
        for (int w = 0; w < npw; w++) if (cols[j][w]) { z = 0; break; }
        int nz = 0;
        for (int w = 0; w < nw; w++) if (ids[j][w]) { nz = 1; break; }
        if (z && nz) {
            int sz = 0;
            for (int i = 0; i < ns; i++)
                if ((ids[j][i/64] >> (i%64)) & 1) sz++;
            out[nr].idx = malloc(sizeof(int) * sz);
            out[nr].n = sz;
            int k = 0;
            for (int i = 0; i < ns; i++)
                if ((ids[j][i/64] >> (i%64)) & 1) out[nr].idx[k++] = i;
            nr++;
        }
    }
    for (int j = 0; j < ns; j++) { free(cols[j]); free(ids[j]); }
    free(cols); free(ids); free(piv);
    *nrel = nr;
    return out;
}

static int extrair(const mpz_t N, smooth_t* sms, rel_t* r, fb_t* fb, mpz_t p, mpz_t q) {
    mpz_t X, Y, t, dif, soma, fac;
    mpz_inits(X, Y, t, dif, soma, fac, NULL);
    mpz_set_ui(X, 1); mpz_set_ui(Y, 1);
    int* etot = calloc(fb->n_primos, sizeof(int));
    for (int k = 0; k < r->n; k++) {
        int idx = r->idx[k];
        mpz_mul(X, X, sms[idx].a); mpz_mod(X, X, N);
        for (int i = 0; i < fb->n_primos; i++) etot[i] += sms[idx].exps[i];
    }
    int ok = 1;
    for (int i = 1; i < fb->n_primos; i++) {
        if (etot[i] & 1) { ok = 0; break; }
        if (etot[i] > 0) {
            mpz_set_ui(t, fb->primos[i]);
            mpz_powm_ui(t, t, etot[i]/2, N);
            mpz_mul(Y, Y, t); mpz_mod(Y, Y, N);
        }
    }
    free(etot);
    int suc = 0;
    if (ok) {
        mpz_sub(dif, X, Y); mpz_mod(dif, dif, N);
        mpz_gcd(fac, dif, N);
        if (mpz_cmp_ui(fac, 1) > 0 && mpz_cmp(fac, N) < 0) {
            mpz_set(p, fac); mpz_divexact(q, N, fac); suc = 1;
        } else {
            mpz_add(soma, X, Y); mpz_mod(soma, soma, N);
            mpz_gcd(fac, soma, N);
            if (mpz_cmp_ui(fac, 1) > 0 && mpz_cmp(fac, N) < 0) {
                mpz_set(p, fac); mpz_divexact(q, N, fac); suc = 1;
            }
        }
    }
    mpz_clears(X, Y, t, dif, soma, fac, NULL);
    return suc;
}

static double now_s(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

int main(int argc, char** argv) {
    setvbuf(stdout, NULL, _IONBF, 0);

    int bits = (argc > 1) ? atoi(argv[1]) : 60;
    int B = (argc > 2) ? atoi(argv[2]) : 2000;
    int M = (argc > 3) ? atoi(argv[3]) : 100000;

    gmp_randstate_t rs;
    gmp_randinit_mt(rs);
    gmp_randseed_ui(rs, 2026);

    mpz_t p, q, N, pf, qf;
    mpz_inits(p, q, N, pf, qf, NULL);
    int half = bits / 2;
    mpz_urandomb(p, rs, half);
    mpz_setbit(p, half - 1);
    mpz_nextprime(p, p);
    mpz_urandomb(q, rs, half);
    mpz_setbit(q, half - 1);
    mpz_nextprime(q, q);
    mpz_mul(N, p, q);

    char Nstr[128], pstr[128], qstr[128];
    gmp_sprintf(Nstr, "%Zd", N);
    gmp_sprintf(pstr, "%Zd", p);
    gmp_sprintf(qstr, "%Zd", q);

    printf("Factoring N of %d bits = %s x %s\n", bits, pstr, qstr);
    printf("N = %s\n", Nstr);
    printf("B=%d, M=+/-%d\n\n", B, M);

    double t0 = now_s();
    fb_t* fb = construir_fb(N, B);
    double t_fb = now_s() - t0;
    printf("[%.0fms] Factor base built: %d primes\n", t_fb*1000, fb->n_primos);

    int ns;
    smooth_t* sms = sieve(N, fb, M, &ns);
    double t_sieve = now_s() - t0;
    printf("[%.0fms] Sieve complete: %d smooths (need > %d)\n", t_sieve*1000, ns, fb->n_primos);

    if (ns <= fb->n_primos) {
        printf("\nFAIL: Insufficient smooths - increase B or M\n");
        return 1;
    }

    int nr;
    rel_t* rels = kernel(sms, ns, fb, &nr);
    double t_kernel = now_s() - t0;
    printf("[%.0fms] Kernel: %d relations\n", t_kernel*1000, nr);

    int suc = 0;
    int rel_usada = -1;
    for (int i = 0; i < nr && !suc; i++) {
        suc = extrair(N, sms, &rels[i], fb, pf, qf);
        if (suc) rel_usada = i;
    }
    double t_total = now_s() - t0;

    if (suc) {
        char pfs[128], qfs[128];
        gmp_sprintf(pfs, "%Zd", pf); gmp_sprintf(qfs, "%Zd", qf);
        mpz_t chk; mpz_init(chk); mpz_mul(chk, pf, qf);
        int v = (mpz_cmp(chk, N) == 0);
        printf("\n[%.0fms] OK FACTORED (relation %d/%d)\n", t_total*1000, rel_usada+1, nr);
        printf("  p = %s\n  q = %s\n  validation: %s\n", pfs, qfs, v ? "OK" : "INVALID");
        mpz_clear(chk);
    } else {
        printf("\n[%.0fms] FAIL: No relation produced non-trivial factor\n", t_total*1000);
    }

    return 0;
}
