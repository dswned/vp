; ~27c / 8 yuv pixels
;
;
%define vsize 32
%include "regs.asm"

global lab2yuv_f32_ivb
extern lab2yuv_f32_ivb_data

section .text

lab2yuv_f32_ivb:
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

vmovups        F,[r1+r4]
vmovups        H,[r1+r5]
vmovups        N,[r1+r6]

vbroadcastss   G,[r0-128]
vbroadcastss   C,[r0-124]
vbroadcastss   D,[r0-120]
vbroadcastss   E,[r0-116]
vbroadcastss   K,[r0-112]
vbroadcastss   R,[r0-108]
vbroadcastss   U,[r0-104]
vbroadcastss   V,[r0-100]

align 16
.1:
add    r1,32

vcmpgtps       B,F,C
vmulps         P,F,[r0-60h]

vmulps         F,F,D
vaddps         F,F,G

vmulps         I,H,[r0-40h]
vaddps         I,I,F

vmulps         T,N,[r0-20h]
vaddps         T,T,F

vmulps         A,F,F
vmulps         F,F,A

vmulps         H,H,[r0]
vaddps         H,P,H

vmulps         N,N,[r0+20h]
vaddps         N,N,P

vblendvps      P,P,F,B
vmovups        F,[r1+r4]

vcmpgtps       B,I,E
vmulps         A,I,I
vmulps         I,I,A
vblendvps      H,H,I,B

vcmpgtps       B,T,E
vmulps         A,T,T
vmulps         T,T,A
vblendvps      N,N,T,B

vmulps         I,P,[r0+40h]
vmulps         A,K,H
vaddps         I,I,A

vmulps         T,P,[r0+60h]
vmulps         A,R,H
vaddps         T,T,A

vmulps         A,N,U
vaddps         I,I,A

vmulps         A,N,V
vaddps         T,T,A

vmovups        H,[r1+r5]
vmovups        N,[r1+r6]
vmovups        [r1   -20h],P
vmovups        [r1+r2-20h],I
vmovups        [r1+r3-20h],T
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