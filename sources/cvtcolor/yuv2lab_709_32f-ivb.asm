; L, a, b scaled to fit rgb in 0..1,-.5+.5
; Ivy Bridge specific code, ~78c / 4 yuv pixels
;
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
export yuv2lab_709_32f_ivb
%endif

extern cvtcolor_dt
global yuv2lab_709_32f_ivb

section .text

yuv2lab_709_32f_ivb:
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
lea            r0,[rel cvtcolor_dt+1312]

vbroadcastf128 F,[r1+r4]
vbroadcastf128 H,[r1+r5]
vmulps         H,H,[r0-60h]
vaddps         H,H,F
vbroadcastf128 N,[r6+r1]
vmulps         N,N,[r0+20h]
vaddps         N,N,H
vextractf128   v,N,1

lea            r7,[rel cvtcolor_dt+160]
sub            r2,r1
sub            r3,r1

align 16
.1:
add            r1,16

vbroadcastss   a,[r7-32]
vpsrld         b,f,23
vpand          c,f,a
vpand          h,a,n

vbroadcastss   t,[r7-36]
vpsrld         g,n,23
vpand          r,v,a
vpand          p,t,f

vbroadcastss   a,[r7-28]
vpsrld         d,v,23
vpand          i,n,t
vpand          t,v,t

vpsrldq        e,f,2
vpor           p,p,a
vpor           i,i,a
vpor           t,t,a
vbroadcastss   a,[r7-24]
vpsrldq        k,n,2
vpor           c,c,a
vpor           h,h,a
vpor           r,r,a
vbroadcastss   a,[r7-40]
vpsrldq        u,v,2
vpand          e,e,a
vpand          k,k,a
vpand          u,u,a

vbroadcastss   a,[r7-20]
vsubps         c,c,p
vpmulld        p,b,a
vpsrld         p,p,12

vsubps         h,h,i
vpmulld        i,g,a
vpsrld         i,i,12

vsubps         r,r,t
vpmulld        t,d,a
vpsrld         t,t,12

vinsertf128    H,H,r,1

vmovd          eax,e
vmovd          a,[r7+rax]
vpextrd        eax,e,1
vpinsrd        a,a,[r7+rax],1
vpextrd        eax,e,2
vpinsrd        a,a,[r7+rax],2
vpextrd        eax,e,3
vpinsrd        r,a,[r7+rax],3

vmovd          eax,k
vmovd          a,[r7+rax]
vpextrd        eax,k,1
vpinsrd        a,a,[r7+rax],1
vpextrd        eax,k,2
vpinsrd        a,a,[r7+rax],2
vpextrd        eax,k,3
vpinsrd        v,a,[r7+rax],3

vmovd          eax,u
vmovd          a,[r7+rax]
vpextrd        eax,u,1
vpinsrd        a,a,[r7+rax],1
vpextrd        eax,u,2
vpinsrd        a,a,[r7+rax],2
vpextrd        eax,u,3
vpinsrd        a,a,[r7+rax],3
vinsertf128    V,V,a,1

vbroadcastss   A,[r7-48]
vmulps         c,r,c
vpsubd         b,b,p
vpsubd         g,g,i
vmulps         r,c,a
vpsubd         d,d,t
vpsubd         b,b,p
vmulps         H,H,V
vpsubd         g,g,i
vpsubd         d,d,t
vmulps         V,H,A
vpsubd         b,b,p
vpaddd         p,p,[r7-40h]
vpslld         b,b,7
vpsubd         g,g,i
vpaddd         i,i,[r7-40h]
vpslld         g,g,7
vpsubd         d,d,t
vpaddd         t,t,[r7-40h]
vpslld         d,d,7
vpaddd         b,b,e
vpaddd         g,k,g
vpaddd         d,u,d
vpslld         p,p,23
vpcmpgtd       a,b,[r7-10h]
vpand          b,b,a
vpslld         i,i,23
vpcmpgtd       a,g,[r7-10h]
vpand          g,g,a
vpslld         t,t,23
vpcmpgtd       a,d,[r7-10h]
vpand          d,d,a
vinsertf128    I,I,t,1

vmovaps        A,[r7-80h]
vaddps         r,r,a
vaddps         V,V,A

vmovaps        A,[r7-60h]
vmulps         r,r,c
vaddps         t,r,a
vmulps         V,V,H
vaddps         V,V,A

vmovd          eax,b
vmovq          e,[r7+rax*2]
vpextrd        eax,b,1
vpinsrq        e,e,[r7+rax*2],1
vpextrd        eax,b,2
vmovq          k,[r7+rax*2]
vpextrd        eax,b,3
vpinsrq        b,k,[r7+rax*2],1

vshufps        k,e,b,0x88
vshufps        e,e,b,0xdd

vmovd          eax,g
vmovq          a,[r7+rax*2]
vpextrd        eax,g,1
vpinsrq        a,a,[r7+rax*2],1
vmovd          eax,d
vmovq          b,[r7+rax*2]
vpextrd        eax,d,1
vpinsrq        b,b,[r7+rax*2],1
vinsertf128    B,A,b,1

vpextrd        eax,g,2
vmovq          a,[r7+rax*2]
vpextrd        eax,g,3
vpinsrq        a,a,[r7+rax*2],1
vpextrd        eax,d,2
vmovq          g,[r7+rax*2]
vpextrd        eax,d,3
vpinsrq        g,g,[r7+rax*2],1

vinsertf128    A,A,g,1
vmulps         k,k,p

vmulps         t,t,c

vmulps         V,V,H
vshufps        G,B,A,0x88

vmulps         G,G,I
vbroadcastf128 H,[r1+r5]

vmulps         e,p,e
vshufps        B,B,A,0xdd

vmulps         B,B,I
vbroadcastf128 R,[r6+r1]

vmulps         t,t,k
vaddps         t,t,e
vmulps         V,V,G
vaddps         V,V,B
vbroadcastss   B,[r0-32]

vmulps         p,f,b
vaddps         t,t,k

vbroadcastss   A,[r0+80]
vmulps         I,N,B
vaddps         V,V,G
vmulps         H,H,[r0-60h]

vaddps         p,a,p
vcmpgtps       b,f,[r0-40h]

vbroadcastf128 F,[r1+r4]

vaddps         I,I,A
vcmpgtps       N,N,[r0-40h]
vmulps         R,R,[r0+20h]

vpand          t,t,b
vandps         V,V,N

vaddps         H,H,F
vpandn         b,b,p
vandnps        I,N,I

vpor           t,t,b
vorps          I,I,V

vextractf128   c,I,1
vsubps         i,i,t
vmulps         p,t,[r0-70h]
vsubps         t,t,c
vmulps         i,i,[r0+10h]
vaddps         N,R,H
vextractf128   v,N,1

vmulps         t,t,[r0-10h]
vsubps         p,p,[r0+40h]

vmovups        [r1   -10h],p
vmovups        [r1+r2-10h],i
vmovups        [r1+r3-10h],t
cmp            r10,r1
ja             .1
; todo: last 4-7
; % 0 - last 4; 1 - dec ptr, rep & last 4; 2...

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