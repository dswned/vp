%include "yregs.asm"
%ifidn __OUTPUT_FORMAT__, win64
%define WINABI
%define r0 rcx
%define r1 rdx
%define r2 r8
%define r3 r9
%define r4 rsi
%define r5 rdi
%define r6 rbp
%define r7 rbx
%else
%define r0 rdi
%define r1 rsi
%define r2 rdx
%define r3 rcx
%define r4 r8
%define r5 r9
%define r6 rbp
%define r7 rbx
%endif

%ifdef WINABI
export lab2yuv_709_32f_ivb
%endif

extern cvtcolor_dt
global lab2yuv_709_32f_ivb

section .text

lab2yuv_709_32f_ivb:
push           rbp
push           rbx

%ifdef WINABI
push           rsi
push           rdi
mov            r4,[rsp+72]
mov            r5,[rsp+80]
mov            r6,[rsp+88]
sub            rsp,88
%else
mov            r6,[rsp+24]
%endif

sub            r4,r1
sub            r5,r1
sub            r6,r1

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

lea            r10,[r0+r1]
lea            r0,[rel cvtcolor_dt+1520]
sub            r2,r1
sub            r3,r1

vmovaps        E,[r0+30h]
vmovaps        K,[r0+70h]
vbroadcastss   U,[r0-236]
vmovaps        c,[r0-10h]
vmovaps        H,[r0+50h]
vmovaps        R,[r0-30h]
vmovups        f,[r1+r4]
vmovups        n,[r1+r5]
vinsertf128    N,N,[r6+r1],1

align 16
.1:
add            r1,16

vmulps         p,f,c
vaddps         p,p,[r0-80h]
vmulps         b,f,[r0]
vcmpgtps       f,f,[r0-70h]

vinsertf128    A,P,p,1
vmulps         G,N,E
vaddps         G,G,A

vmulps         v,p,p
vmulps         v,v,p

vinsertf128    A,B,b,1
vmulps         N,N,K
vaddps         N,N,A

vcmpgtps       I,G,U
vmulps         T,G,G
vmulps         T,T,G

vandps         T,T,I
vpand          p,v,f
vandnps        I,I,N
vpandn         f,f,b
vorps          T,T,I
vpor           p,p,f

vmulps         n,p,[r0-60h]
vmulps         I,T,H
vextractf128   g,I,1
vaddps         i,i,n
vmulps         v,p,[r0+20h]
vmulps         T,T,R
vmovups        f,[r1+r4]
vextractf128   d,T,1
vaddps         t,t,v

vmovups        n,[r1+r5]
vaddps         i,i,g
vaddps         t,t,d
vinsertf128    N,N,[r6+r1],1

vmovups        [r1   -10h],p
vmovups        [r1+r2-10h],i
vmovups        [r1+r3-10h],t
cmp            r10,r1
ja             .1

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
add            rsp,88
pop            rdi
pop            rsi
%endif

pop            rbx
pop            rbp
vzeroupper
ret
align 16