; ~18+?c / 16 yuv pixels
;
;
%define vsize 64
%include "regs.asm"

global lab2yuv_f32_avx512
extern lab2yuv_f32_ivb_data

section .text

lab2yuv_f32_avx512:
push   rbp
push   rbx

%ifdef WINABI
push   rsi
push   rdi
mov    r4,[rsp+72]
mov    r5,[rsp+80]
mov    r6,[rsp+88]
sub    rsp,88
%else
mov    r6,[rsp+24]
%endif

sub    r4,r1
sub    r5,r1
sub    r6,r1

%ifdef WINABI
vmovups [rsp-80],xmm15
vmovups [rsp-64],xmm14
vmovups [rsp-48],xmm13
vmovups [rsp-32],xmm12
vmovups [rsp-16],xmm11
vmovups [rsp   ],xmm10
vmovups [rsp+16], xmm9
vmovups [rsp+32], xmm8
vmovups [rsp+48], xmm7
vmovups [rsp+64], xmm6
%endif

lea    r10,[r0+r1]
lea    r0,[rel lab2yuv_f32_ivb_data+80h]
sub    r2,r1
sub    r3,r1

vmovups        B,[r1+r4]
vmovups        G,[r1+r5]
vmovups        D,[r1+r6]

vbroadcastss   C,[r0-128]

align 16
.1:
add    r1,64

vpcmpgtd       k1,B,[r0-124]{1to16}
vmulps         P,B,[r0-60h]{1to16}
vfmadd132ps    B,C,[r0-120]{1to16}

vmovaps        H,G
vmovaps        R,D

vfmadd132ps    G,P,[r0]{1to16}
vfmadd132ps    D,P,[r0+20h]{1to16}
vfmadd132ps    H,B,[r0-40h]{1to16}
vfmadd132ps    R,B,[r0-20h]{1to16}

vmulps         A,B,B
vmulps         B,B,A
vpcmpgtd       k2,H,[r0-116]{1to16}
vmulps         A,H,H
vmulps         H,H,A
vpcmpgtd       k3,R,[r0-116]{1to16}
vmulps         A,R,R
vmulps         R,R,A

vmovdqa32      P{k1},B
vmovdqa32      G{k2},H
vmovdqa32      D{k3},R

vmulps         I,P,[r0+40h]{1to16}
vfmadd231ps    I,G,[r0-112]{1to16}
vfmadd231ps    I,D,[r0-104]{1to16}

vmulps         T,P,[r0+60h]{1to16}
vfmadd231ps    T,G,[r0-108]{1to16}
vfmadd231ps    T,D,[r0-100]{1to16}

vmovups        B,[r1+r4]
vmovups        G,[r1+r5]
vmovups        D,[r1+r6]
vmovups        [r1   -40h],P
vmovups        [r1+r2-40h],I
vmovups        [r1+r3-40h],T
cmp    r10,r1
ja     .1

%ifdef WINABI
vmovups xmm15,[rsp-80]
vmovups xmm14,[rsp-64]
vmovups xmm13,[rsp-48]
vmovups xmm12,[rsp-32]
vmovups xmm11,[rsp-16]
vmovups xmm10,[rsp   ]
vmovups xmm9, [rsp+16]
vmovups xmm8, [rsp+32]
vmovups xmm7, [rsp+48]
vmovups xmm6, [rsp+64]
add    rsp,88
pop    rdi
pop    rsi
%endif

pop    rbx
pop    rbp
vzeroupper
ret
align 16