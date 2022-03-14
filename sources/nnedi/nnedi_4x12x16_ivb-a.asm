;
;
;
%define vsize 32
%include "regs.asm"

%macro PFMA 7
vpermilps B, %2, %3
vmulps A, B, [r7+64*%1+32*0+00A0h]
vaddps %4, %4, A
vmulps A, B, [r7+64*%1+32*1+00A0h]
vaddps %5, %5, A
vmulps A, B, [r7+64*%1+32*0+0CA0h]
vaddps %6, %6, A
vmulps A, B, [r7+64*%1+32*1+0CA0h]
vaddps %7, %7, A
%endmacro

%macro PFMA_0 3
PFMA %{1:3},G,I,C,E
%endmacro

%macro PFMA_1 3
PFMA %{1:3},F,K,D,H
%endmacro

section .rdata align=64

k_80000000: times 8 dd 0x80000000
k_7fffffff: times 8 dd 0x7fffffff
k_3f800000: times 8 dd 0x3f800000
k_20000000: dd 0x20000000

global nnedi_4x12x16_ivb

section .text

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 35-254c 200fma+12load+19p5+
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
nnedi_4x12x16_ivb:
push           rbp
push           rbx

%ifdef WINABI
push           rsi
push           rdi
sub            rsp,88
%else
%endif

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

lea r4, [r0+12*4]
lea r6, [r4+r4*2]
lea r7, [r3+160]
sub r1, r2
add r0, r2

align 16
.1:
; input = 00 01 02 03
;         04 05 06 07
;         08 09 10 11
;         12 13 14 15
; x  = [00:01:02:03|00:01:02:03] * [aaaa|bbbb]
; y  = [00:01:02:03|00:01:02:03] * [cccc|dddd]
; x += [04:05:06:07|04:05:06:07] * [aaaa|bbbb]
;
; x        = [a00a04a08a12:a01a05a09a13:a02a06a10a14:a03a07a11a15|b00b04b08b12:b01b05b09b13:b02b06b10b14:b03b07b11b15]
;
; unpcklps = [a00a04a08a12:c00c04c08c12:a01a05a09a13:c01c05c09c13|b00b04b08b12:d00d04d08d12:b01b05b09b13:d01d05d09d13]
; unpckhps = [a02a06a10a14:c02c06c10c14:a03a07a11a15:c03c07c11c15|b02b06b10b14:d02d06d10d14:b03b07b11b15:d03d07d11d15]
;
vbroadcastf128 P, [r2+16]
vmulps      F, P, [r7+32*0-80h]
vbroadcastf128 R, [r2+r4*1+16]
vmulps      H, P, [r7+32*1-80h]
vbroadcastf128 U, [r2+r4*2+16]
vbroadcastf128 V, [r2+r6*1+16]
vmulps      A, R, [r7+32*2-80h]
vaddps      F, F, A
vmulps      A, R, [r7+32*3-80h]
vaddps      H, H, A
vmulps      A, U, [r7+32*4-80h]
vaddps      F, F, A
vmulps      A, U, [r7+32*5-80h]
vaddps      H, H, A
vmulps      A, V, [r7+32*6-80h]
vaddps      F, F, A
vmulps      A, V, [r7+32*7-80h]
vaddps      H, H, A
vinsertps   t, f, f, 93
vxorps      F, F, T
vunpcklps   A, F, H
vunpckhps   F, F, H
vaddps      F, F, A
vxorps      A, F, [rel k_80000000]
vunpcklpd   H, F, A
vunpckhpd   F, F, A
vaddps      F, F, H
vbroadcastss   A, [r7-160]; bias (th1)
vaddps      F, F, A
vandps      A, F, [rel k_7fffffff]
vaddps      A, A, [rel k_3f800000]
vrcpps      A, A
vmulps      F, F, [r7+80h]; pk2
vmulps      F, F, A
vextractf128   a, F, 1
vaddps      f, f, a
vmovhlps    a, f, f
vaddps      f, f, a
vshufps     a, f, f, 1
vaddss      a, a, f

.2:
vcomiss     a, [r7-156]; bias (th2)
jbe .3

;vmovups t, [rel k_3f800000]
vextractps  [r1+r2], t, 1
add r2,4
cmp r0,r2
je  .4
jmp .1

align 16
.3:
vpermilps      E, P, 0
vmulps      G, E, [r7+64*4+32*0+00A0h]
vmulps      I, E, [r7+64*4+32*1+00A0h]
vmulps      C, E, [r7+64*4+32*0+0CA0h]
vmulps      E, E, [r7+64*4+32*1+0CA0h]
vpermilps      H, P, 055h
vmulps      F, H, [r7+64*5+32*0+00A0h]
vmulps      K, H, [r7+64*5+32*1+00A0h]
vmulps      D, H, [r7+64*5+32*0+0CA0h]
vmulps      H, H, [r7+64*5+32*1+0CA0h]
PFMA_0  6, P, 0AAh
PFMA_1  7, P, 0FFh

vbroadcastf128 P, [r2]
PFMA_0 16, R, 0
PFMA_1 17, R, 055h
PFMA_0 18, R, 0AAh
PFMA_1 19, R, 0FFh

vbroadcastf128 R, [r2+r4*1]
PFMA_0 28, U, 0
PFMA_1 29, U, 055h
PFMA_0 30, U, 0AAh
PFMA_1 31, U, 0FFh

vbroadcastf128 U, [r2+r4*2]
PFMA_0 40, V, 0
PFMA_1 41, V, 055h
PFMA_0 42, V, 0AAh
PFMA_1 43, V, 0FFh

vbroadcastf128 V, [r2+r6*1]
add r2, 32
PFMA_0  0, P, 0
PFMA_1  1, P, 055h
PFMA_0  2, P, 0AAh
PFMA_1  3, P, 0FFh

vbroadcastf128 P, [r2]
PFMA_0 12, R, 0
PFMA_1 13, R, 055h
PFMA_0 14, R, 0AAh
PFMA_1 15, R, 0FFh

vbroadcastf128 R, [r2+r4*1]
PFMA_0 24, U, 0
PFMA_1 25, U, 055h
PFMA_0 26, U, 0AAh
PFMA_1 27, U, 0FFh

vbroadcastf128 U, [r2+r4*2]
PFMA_0 36, V, 0
PFMA_1 37, V, 055h
PFMA_0 38, V, 0AAh
PFMA_1 39, V, 0FFh

vbroadcastf128 V, [r2+r6*1]
sub r2, 12
PFMA_0  8, P, 0
PFMA_1  9, P, 055h
PFMA_0 10, P, 0AAh
PFMA_1 11, P, 0FFh

	vbroadcastf128 P, [r2]
PFMA_0 20, R, 0
PFMA_1 21, R, 055h
PFMA_0 22, R, 0AAh
PFMA_1 23, R, 0FFh

	vbroadcastf128 R, [r2+r4*1]
PFMA_0 32, U, 0
PFMA_1 33, U, 055h
PFMA_0 34, U, 0AAh
PFMA_1 35, U, 0FFh

vaddps      G, G, F
	vmulps      F, P, [r7+32*0-80h]
vaddps      E, E, H
	vmulps      H, P, [r7+32*1-80h]
vpermilps      B, V, 0
vmulps      A, B, [r7+64*44+32*0+00A0h]
vaddps      G, G, A
vmulps      A, B, [r7+64*44+32*1+00A0h]
vaddps      I, I, A
vaddps      I, I, K
vpermilps      K, V, 0AAh
vmulps      A, B, [r7+64*44+32*0+0CA0h]
vaddps      C, C, A
vaddps      C, C, D
vpermilps      D, V, 0FFh
vmulps      A, B, [r7+64*44+32*1+0CA0h]
vaddps      E, E, A

	vbroadcastf128 U, [r2+r4*2]
	vmulps      A, R, [r7+32*2-80h]
	vaddps      F, F, A
	vmulps      A, R, [r7+32*3-80h]
	vaddps      H, H, A
vpermilps      B, V, 055h
vmulps      A, B, [r7+64*45+32*0+00A0h]
vaddps      G, G, A
vmulps      A, B, [r7+64*45+32*1+00A0h]
vaddps      I, I, A
vmulps      A, B, [r7+64*45+32*0+0CA0h]
vaddps      C, C, A
vmulps      A, B, [r7+64*45+32*1+0CA0h]
vaddps      E, E, A

	vbroadcastf128 V, [r2+r6*1]
	vmulps      A, U, [r7+32*4-80h]
	vaddps      F, F, A
	vmulps      A, U, [r7+32*5-80h]
	vaddps      H, H, A
vmulps      A, K, [r7+64*46+32*0+00A0h]
vaddps      G, G, A
vmulps      A, K, [r7+64*46+32*1+00A0h]
vaddps      I, I, A
vmulps      A, K, [r7+64*46+32*0+0CA0h]
vaddps      C, C, A
vmulps      A, K, [r7+64*46+32*1+0CA0h]
vaddps      E, E, A

	vmulps      A, V, [r7+32*6-80h]
	vaddps      F, F, A
	vmulps      A, V, [r7+32*7-80h]
	vaddps      H, H, A
vmulps      A, D, [r7+64*47+32*0+00A0h]
vaddps      G, G, A
vmulps      A, D, [r7+64*47+32*1+00A0h]
vaddps      I, I, A
vmulps      A, D, [r7+64*47+32*0+0CA0h]
vaddps      C, C, A
vmulps      A, D, [r7+64*47+32*1+0CA0h]
vaddps      E, E, A

	vinsertps   t, f, f, 93
	vxorps      F, F, T
vmulps      G, G, G
vmulps      I, I, I
	vunpcklps   A, F, H
	vunpckhps   F, F, H
	vaddps      F, F, A
vmulps      C, C, G
vaddps      G, G, I

	vxorps      A, F, [rel k_80000000]
	vunpcklpd   H, F, A
	vunpckhpd   F, F, A
	vaddps      F, F, H
vmulps      E, E, I
vaddps      E, E, C

	vbroadcastss   A, [r7-160]
	vaddps      F, F, A
	vandps      A, F, [rel k_7fffffff]
	vaddps      A, A, [rel k_3f800000]
	vrcpps      A, A
vhaddps     G, G, E

	vmulps      F, F, [r7+80h]
	vmulps      F, F, A
vextractf128   a, G, 1
vaddps      g, g, a
vhaddps     g, g, g

	vextractf128   a, F, 1
	vaddps      f, f, a
;	vhaddps     f, f, f
;	vhaddps     a, f, f
	vmovhlps    a, f, f
	vaddps      f, f, a
	vshufps     a, f, f, 1
	vaddss      a, a, f
vmovshdup   e, g
vmaxss      g, g, [rel k_20000000]

vdivss      e, e, g
vmovd       [r1+r2-20], e
sub r2, 16
cmp r0, r2
je  .4
jmp .2

.4:

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