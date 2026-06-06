#include <immintrin.h>


void old_intrinsic_grad(double *q, int ci, int ck, double *cxyz, double *ccgradq){

// const u_int64_t bmt = 0xFFFFFFFFFFFFFFFF; // bit mask true
// const u_int64_t bmf = 0x0000000000000000; // bit mask false

const __mmask8 bitmask3 = 0b00000111; 
const __mmask8 bitmask5 = 0b00011111; 
// const int switch_l128   = 0b11111111;

__m512d qi = _mm512_maskz_loadu_pd (bitmask5, &q[(ci-1)*5]);
__m512d qk = _mm512_maskz_loadu_pd (bitmask5, &q[(ck-1)*5]);

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

#define __BITMASK5 0x1F //0b00011111


void intrinsic_grad(double *q, int ci, int nnghbrs, int* ckn, double *cf, double *gradq){

    // const __mmask8 bitmask5 = 0b00011111; 
    const __m512i vind_x = _mm512_set_epi64 ( 0, 0, 0, 12,  9, 6, 3, 0);
    const __m512i vind_y = _mm512_set_epi64 ( 0, 0, 0, 13, 10, 7, 4, 1);
    const __m512i vind_z = _mm512_set_epi64 ( 0, 0, 0, 14, 11, 8, 5, 2);

    __m512d qi = _mm512_maskz_loadu_pd (__BITMASK5, &q[(ci-1)*5]);
    __m512d grdQx = {};
    __m512d grdQy = {};
    __m512d grdQz = {};

    for (int kcell =0; kcell < nnghbrs; kcell++){
        int ck = ckn[kcell];
        __m512d qk = _mm512_maskz_loadu_pd (__BITMASK5, &q[(ck-1)*5]);
        __m512d dq = _mm512_sub_pd (qk, qi);   // dq = qk - qi

        grdQx = _mm512_fmadd_pd(dq, _mm512_set1_pd (cf[(kcell*3) + 0]) , grdQx);
        grdQy = _mm512_fmadd_pd(dq, _mm512_set1_pd (cf[(kcell*3) + 1]) , grdQy);
        grdQz = _mm512_fmadd_pd(dq, _mm512_set1_pd (cf[(kcell*3) + 2]) , grdQz);

    }

    // Transpose (using scatter) and write the gradient back to memory
    _mm512_mask_i64scatter_pd (&gradq[0], __BITMASK5, vind_x, grdQx, 8);
    _mm512_mask_i64scatter_pd (&gradq[0], __BITMASK5, vind_y, grdQy, 8);
    _mm512_mask_i64scatter_pd (&gradq[0], __BITMASK5, vind_z, grdQz, 8);
    
}

void transpose_3x5(double *cT, double *c){
    __m512d __cT1 = _mm512_loadu_pd(&cT[0]);
    __m512d __cT2 = _mm512_maskz_loadu_pd(0b01111111,&cT[8]);
    // Top four bits (4-7) do nothing.  if p1[3] = 1 we draw from c2, 0 we draw from c1
    // bottom three bits (0-2) are the destination in the destination register 
    // funny thing is if we store c1 and c2 as a single 15 item array (where c2[0] = c[8])
    // then we can just use the index as are value.
    const __m512i p1 = _mm512_setr_epi64( 0, 5, 10,  1, 6, 11,  2,  7);
    const __m512i p2 = _mm512_setr_epi64(12, 3,  8, 13, 4,  9, 14, 15);
    
    __m512d c1 = _mm512_permutex2var_pd(__cT1, p1, __cT2);
    __m512d c2 = _mm512_permutex2var_pd(__cT1, p2, __cT2);

    _mm512_storeu_pd (&c[0], c1);
    _mm512_mask_storeu_pd(&c[8], 0b01111111, c2);

}