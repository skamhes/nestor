#include <immintrin.h>

void intrinsic_grad(double *q, int ci, int ck, double *cxyz, double *ccgradq){

// const u_int64_t bmt = 0xFFFFFFFFFFFFFFFF; // bit mask true
// const u_int64_t bmf = 0x0000000000000000; // bit mask false

const __mmask8 bitmask3 = 0b00000111; 
const __mmask8 bitmask5 = 0b00011111; 
// const int switch_l128   = 0b11111111;

__m512d qi = _mm512_maskz_loadu_pd (bitmask5, &q[ci*5]);
__m512d qk = _mm512_maskz_loadu_pd (bitmask5, &q[ck*5]);

qi = _mm512_sub_pd (qk, qi);   // dq = qk - qi

__m256d rc = _mm256_maskz_loadu_pd (bitmask3, &cxyz[0]);

__m256d grad[5]; // zero initialize

grad[0] = _mm256_maskz_loadu_pd (bitmask3, &ccgradq[0]);
grad[1] = _mm256_maskz_loadu_pd (bitmask3, &ccgradq[3]);
grad[2] = _mm256_maskz_loadu_pd (bitmask3, &ccgradq[6]);
grad[3] = _mm256_maskz_loadu_pd (bitmask3, &ccgradq[9]);
grad[4] = _mm256_maskz_loadu_pd (bitmask3, &ccgradq[12]);

__m128d q01 = _mm512_extractf64x2_pd (qi, 0); // move lowest 128 bits of qi to q01
__m256d q0 = _mm256_maskz_broadcastsd_pd (bitmask3, q01); // broadcast lower 64 bits of q01 to all threads of q0
grad[0] = _mm256_maskz_fmadd_pd (bitmask3, rc, q0, grad[0]); // calculate first column of grad

q01 = _mm_maskz_permute_pd (bitmask3, q01, 3); // swap lower two qwords
__m256d q1 = _mm256_maskz_broadcastsd_pd (bitmask3, q01);
grad[1] = _mm256_maskz_fmadd_pd (bitmask3, rc, q1, grad[1]);

__m128d q23 = _mm512_extractf64x2_pd (qi, 1); // 0 = lowest 128 bits
__m256d q2 = _mm256_maskz_broadcastsd_pd (bitmask3, q23);
grad[2] = _mm256_maskz_fmadd_pd (bitmask3, rc, q2, grad[2]);

q23 = _mm_maskz_permute_pd (bitmask3, q23, 3); // swap lower two qwords
__m256d q3 = _mm256_maskz_broadcastsd_pd (bitmask3, q23);
grad[3] = _mm256_maskz_fmadd_pd (bitmask3, rc, q3, grad[3]);

__m128d q4x = _mm512_extractf64x2_pd (qi, 2); // 0 = lowest 128 bits
__m256d q4 = _mm256_maskz_broadcastsd_pd (bitmask3, q4x);
grad[4] = _mm256_maskz_fmadd_pd (bitmask3, rc, q4, grad[4]);

_mm256_mask_storeu_pd (&ccgradq[0], bitmask3, grad[0]);
_mm256_mask_storeu_pd (&ccgradq[3], bitmask3, grad[1]);
_mm256_mask_storeu_pd (&ccgradq[6], bitmask3, grad[2]);
_mm256_mask_storeu_pd (&ccgradq[9], bitmask3, grad[3]);
_mm256_mask_storeu_pd (&ccgradq[12], bitmask3, grad[4]);

}