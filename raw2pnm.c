/*
 * Next Video Tool for Linux* OS - Video4Linux2 API tool for developers.
 *
 * Copyright (c) 2017 Intel Corporation. All Rights Reserved.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License version
 * 2 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdarg.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#include <stdint.h>
#include "linux/videodev2.h"

#include "extradefs.h"

#define MIN(a,b)	((a) <= (b) ? (a) : (b))
#define MAX(a,b)	((a) >= (b) ? (a) : (b))
#define CLAMP(a,lo,hi)	((a) < (lo) ? (lo) : (a) > (hi) ? (hi) : (a))
#define CLAMPB(a)	CLAMP(a, 0, 255)

static char *name = "raw2pnm";

static int verbosity = 2;

struct symbol_list {
	int id;
	const char *symbol;
};
#define SYMBOL_END	{ -1, NULL }

#ifndef V4L2_PIX_FMT_YYUV420_V32
#define V4L2_PIX_FMT_YYUV420_V32 v4l2_fourcc('y', 'V', '0', 'A') /* Vectorized NV12 */
#endif

#define V4L2_PIX_FMT	"V4L2_PIX_FMT_"
#define PIXFMT(id)	{ V4L2_PIX_FMT_##id, (#id) }
static const struct symbol_list pixelformats[] = {
	PIXFMT(RGB332),
	PIXFMT(RGB444),
	PIXFMT(ARGB444),
	PIXFMT(XRGB444),
	PIXFMT(RGB555),
	PIXFMT(ARGB555),
	PIXFMT(XRGB555),
	PIXFMT(RGB565),
	PIXFMT(RGB555X),
	PIXFMT(ARGB555X),
	PIXFMT(XRGB555X),
	PIXFMT(RGB565X),
	PIXFMT(BGR666),
	PIXFMT(BGR24),
	PIXFMT(RGB24),
	PIXFMT(BGR32),
	PIXFMT(ABGR32),
	PIXFMT(XBGR32),
	PIXFMT(RGB32),
	PIXFMT(ARGB32),
	PIXFMT(XRGB32),
	PIXFMT(GREY),
	PIXFMT(Y4),
	PIXFMT(Y6),
	PIXFMT(Y10),
	PIXFMT(Y12),
	PIXFMT(Y16),
	PIXFMT(Y10BPACK),
	PIXFMT(PAL8),
	PIXFMT(UV8),
	PIXFMT(YVU410),
	PIXFMT(YVU420),
	PIXFMT(YUYV),
	PIXFMT(YYUV),
	PIXFMT(YVYU),
	PIXFMT(UYVY),
	PIXFMT(VYUY),
	PIXFMT(YUV422P),
	PIXFMT(YUV411P),
	PIXFMT(Y41P),
	PIXFMT(YUV444),
	PIXFMT(YUV555),
	PIXFMT(YUV565),
	PIXFMT(YUV32),
	PIXFMT(YUV410),
	PIXFMT(YUV420),
	PIXFMT(HI240),
	PIXFMT(HM12),
	PIXFMT(M420),
	PIXFMT(NV12),
	PIXFMT(NV21),
	PIXFMT(NV16),
	PIXFMT(NV61),
	PIXFMT(NV24),
	PIXFMT(NV42),
	PIXFMT(NV12M),
	PIXFMT(NV21M),
	PIXFMT(NV16M),
	PIXFMT(NV61M),
	PIXFMT(NV12MT),
	PIXFMT(NV12MT_16X16),
	PIXFMT(YUV420M),
	PIXFMT(YVU420M),
	PIXFMT(SBGGR8),
	PIXFMT(SGBRG8),
	PIXFMT(SGRBG8),
	PIXFMT(SRGGB8),
	PIXFMT(SBGGR10),
	PIXFMT(SGBRG10),
	PIXFMT(SGRBG10),
	PIXFMT(SRGGB10),
	PIXFMT(SBGGR10P),
	PIXFMT(SGBRG10P),
	PIXFMT(SGRBG10P),
	PIXFMT(SRGGB10P),
	PIXFMT(SBGGR10ALAW8),
	PIXFMT(SGBRG10ALAW8),
	PIXFMT(SGRBG10ALAW8),
	PIXFMT(SRGGB10ALAW8),
	PIXFMT(SBGGR10DPCM8),
	PIXFMT(SGBRG10DPCM8),
	PIXFMT(SGRBG10DPCM8),
	PIXFMT(SRGGB10DPCM8),
	PIXFMT(SBGGR12),
	PIXFMT(SGBRG12),
	PIXFMT(SGRBG12),
	PIXFMT(SRGGB12),
	PIXFMT(SBGGR12P),
	PIXFMT(SGBRG12P),
	PIXFMT(SGRBG12P),
	PIXFMT(SRGGB12P),
	PIXFMT(SBGGR14),
	PIXFMT(SGBRG14),
	PIXFMT(SGRBG14),
	PIXFMT(SRGGB14),
	PIXFMT(SBGGR14P),
	PIXFMT(SGBRG14P),
	PIXFMT(SGRBG14P),
	PIXFMT(SRGGB14P),
	PIXFMT(SBGGR16),
	PIXFMT(SBGGR8_16V32),
	PIXFMT(SGBRG8_16V32),
	PIXFMT(SGRBG8_16V32),
	PIXFMT(SRGGB8_16V32),
	PIXFMT(SBGGR10V32),
	PIXFMT(SGBRG10V32),
	PIXFMT(SGRBG10V32),
	PIXFMT(SRGGB10V32),
	PIXFMT(SBGGR12V32),
	PIXFMT(SGBRG12V32),
	PIXFMT(SGRBG12V32),
	PIXFMT(SRGGB12V32),
	PIXFMT(SBGGR8V32),
	PIXFMT(SGBRG8V32),
	PIXFMT(SGRBG8V32),
	PIXFMT(SRGGB8V32),
	PIXFMT(UYVY_V32),
	PIXFMT(YUYV420_V32),
	PIXFMT(MJPEG),
	PIXFMT(JPEG),
	PIXFMT(DV),
	PIXFMT(MPEG),
	PIXFMT(H264),
	PIXFMT(H264_NO_SC),
	PIXFMT(H264_MVC),
	PIXFMT(H263),
	PIXFMT(MPEG1),
	PIXFMT(MPEG2),
	PIXFMT(MPEG4),
	PIXFMT(XVID),
	PIXFMT(VC1_ANNEX_G),
	PIXFMT(VC1_ANNEX_L),
	PIXFMT(VP8),
	PIXFMT(CPIA1),
	PIXFMT(WNVA),
	PIXFMT(SN9C10X),
	PIXFMT(SN9C20X_I420),
	PIXFMT(PWC1),
	PIXFMT(PWC2),
	PIXFMT(ET61X251),
	PIXFMT(SPCA501),
	PIXFMT(SPCA505),
	PIXFMT(SPCA508),
	PIXFMT(SPCA561),
	PIXFMT(PAC207),
	PIXFMT(MR97310A),
	PIXFMT(JL2005BCD),
	PIXFMT(SN9C2028),
	PIXFMT(SQ905C),
	PIXFMT(PJPG),
	PIXFMT(OV511),
	PIXFMT(OV518),
	PIXFMT(STV0680),
	PIXFMT(TM6000),
	PIXFMT(CIT_YYVYUY),
	PIXFMT(KONICA420),
	PIXFMT(JPGL),
	PIXFMT(SE401),
	PIXFMT(S5C_UYVY_JPG),
	PIXFMT(SMIAPP_META_8),
	PIXFMT(YYUV420_V32),
	PIXFMT(PRIV_MAGIC),
	PIXFMT(FLAG_PREMUL_ALPHA),
	SYMBOL_END
};

static void print(int lvl, char *msg, ...)
{
	va_list ap;

	if (verbosity < lvl)
		return;

	va_start(ap, msg);
	vprintf(msg, ap);
	va_end(ap);
}

static void error(char *msg, ...)
{
	FILE *f = stdout;
	va_list ap;
	int e = errno;
 
	va_start(ap, msg);
	fprintf(f, "%s: ", name);
	vfprintf(f, msg, ap);
	if (e)
		fprintf(f, ": %s (%i)", strerror(e), e);
	fprintf(f, "\n");
	va_end(ap);
	exit(1);
}

static void usage(void)
{
	print(1,"Usage: %s [-x width] [-y height] [-s stride] [-f format] [-b skip bytes] [inputfile] [outputfile]\n", name);
}

static const char *symbol_str(int id, const struct symbol_list list[])
{
	static char buffer[200];
	int i;

	for (i=0; list[i].symbol; i++)
		if (list[i].id == id)
			break;

	if (list[i].symbol) {
		if (id < 1000)
			sprintf(buffer, "%s [%i]", list[i].symbol, id);
		else
			sprintf(buffer, "%s [0x%08X]", list[i].symbol, id);
	} else {
		sprintf(buffer, "%i", id);
	}
	return buffer;
}

static void symbol_dump(const char *prefix, const struct symbol_list *list)
{
	int i;
	for (i = 0; list[i].symbol != NULL; i++)
		print(0, "%s%s [0x%08X]\n", prefix, list[i].symbol, list[i].id);
}

static int symbol_get(const struct symbol_list *list, const char **symbol)
{
	const char *start;
	const char *end;
	int r, i;

	if (!symbol || !symbol[0])
		error("symbol missing");

	start = *symbol;

	if (isdigit(start[0])) {
		r = strtol(start, (char **)&end, 0);
		if (end == start)
			error("zero-length symbol value");
	} else {
		end = start;
		while (isalnum(*end) || *end == '_') end++;
		if (start == end)
			error("zero-length symbol");
		if (!list)
			error("only numeric value allowed");
		for (i=0; list[i].symbol!=NULL; i++) {
			if (strncasecmp(start, list[i].symbol, end-start) == 0)
				break;
		}
		if (list[i].symbol == NULL)
			error("symbol `%s' not found", start);
		r = list[i].id;
	}

	*symbol = end;
	return r;
}

static long fsize(FILE *f)
{
	int i;
	long p, s, r = 0;

	p = ftell(f);
	if (p == -1)
		return -1;
	i = fseek(f, 0, SEEK_END);
	if (i == -1)
		r = -1;
	s = ftell(f);
	if (s == -1)
		r = -1;
	i = fseek(f, p, SEEK_SET);
	if (i != p)
		r = -1;

	return r == 0 ? s : r;
}

static void *duplicate_buffer(void *buffer, int size, int new_size)
{
	void *b = calloc(1, new_size);
	if (!b)
		error("out of memory, can not allocate %i bytes", new_size);
	memcpy(b, buffer, MIN(size, new_size));
	if (new_size > size)
		print(1, "warning: input buffer too small by %i bytes, setting the rest to zero\n",
			new_size - size);
	return b;
}

static void inline yuv_to_rgb(unsigned char rgb[3], int y, int cb, int cr)
{
	static const int R = 0;
	static const int G = 1;
	static const int B = 2;
	int u = cb;
	int v = cr;
#if 0
	/* http://www.fourcc.org/fccyvrgb.php */
	rgb[B] = CLAMPB(1.164*(y - 16)                   + 2.018*(u - 128));
	rgb[G] = CLAMPB(1.164*(y - 16) - 0.813*(v - 128) - 0.391*(u - 128));
	rgb[R] = CLAMPB(1.164*(y - 16) + 1.596*(v - 128));
#else
	/* http://en.wikipedia.org/wiki/YUV
	 * conversion from Y'UV to RGB (NTSC version): */
	int c = y - 16;
	int d = u - 128;
	int e = v - 128;
	rgb[R] = CLAMPB((298*c + 409*e + 128) >> 8);
	rgb[G] = CLAMPB((298*c - 100*d - 208*e + 128) >> 8);
	rgb[B] = CLAMPB((298*c + 516*d + 128) >> 8);
#endif
}

static inline unsigned int get_word(unsigned char *src, int bpp)
{
	unsigned int lo = src[0];
	if (bpp == 1) {
		return lo;
	} else if (bpp == 2) {
		unsigned int hi = src[1];
		unsigned int w = (hi << 8) | lo;
		return w;
	} else error("bad bpp");
	return 0;
}

static inline unsigned int get_packed_word(unsigned char *src, int x)
{
	int val = 0;
	int shift = 2;
	int xoff = x % 4;
	int bitoff = xoff * 2;

	//printf("off %d, %d\n", xoff, bitoff);
	val = src[0] << shift | (src[4-xoff] && (0x3 << bitoff) >> bitoff);
	//printf("val 0x%x\n", val);
	return val;
}


/* Return 24 bits-per-pixel RGB image */
static int convert(void *in_buffer, int in_size, int width, int height, int stride, __u32 format, void *out_buffer)
{
	static const int dbpp = 3;
	int y, x, r, g, b, bpp, shift, g1, r1, b1;
	int oddrow, oddpix, initrow, initpix, lumaofs, chromaord, subsample;
	unsigned char *src = NULL, *dst = NULL;
	unsigned char *s, *u;
	unsigned char *d = out_buffer;
	unsigned int dstride = width * dbpp;

	switch (format) {
	case V4L2_PIX_FMT_YUYV:
	case V4L2_PIX_FMT_UYVY:
	case V4L2_PIX_FMT_YVYU:
	case V4L2_PIX_FMT_VYUY:
		if (stride <= 0) stride = width * 2;
		s = src = duplicate_buffer(in_buffer, in_size, stride * height);
		lumaofs = (format==V4L2_PIX_FMT_YUYV || format==V4L2_PIX_FMT_YVYU) ? 0 : 1;
		chromaord = (format==V4L2_PIX_FMT_YUYV || format==V4L2_PIX_FMT_UYVY) ? 0 : 1;
		for (y = 0; y < height; y++) {
			unsigned char *s1 = s;
			unsigned char *d1 = d;
			int cb = 0, cr = 0;
			for (x = 0; x < width; x++) {
				int b = s1[lumaofs];
				if ((x & 1) == chromaord)
					cb = s1[lumaofs^1];
				else
					cr = s1[lumaofs^1];
				yuv_to_rgb(d1, b, cb, cr);
				s1 += 2;
				d1 += dbpp;
			}
			s += stride;
			d += dstride;
		}
		break;

	case V4L2_PIX_FMT_NV12:
	case V4L2_PIX_FMT_NV21:
	case V4L2_PIX_FMT_NV24:
	case V4L2_PIX_FMT_NV42:
		subsample = (format == V4L2_PIX_FMT_NV12 ||
			     format == V4L2_PIX_FMT_NV21) ? 2 : 1;
		chromaord = (format == V4L2_PIX_FMT_NV12 ||
			     format == V4L2_PIX_FMT_NV24) ? 0 : 1;
		if (stride <= 0) stride = width;
		s = src = duplicate_buffer(in_buffer, in_size, stride * height*3/subsample);
		u = &s[height * stride];
		for (y = 0; y < height; y++) {
			unsigned char *s1 = s;
			unsigned char *u1 = u;
			unsigned char *d1 = d;
			for (x = 0; x < width; x++) {
				int b = *s1;
				int cb = u1[chromaord];
				int cr = u1[chromaord ^ 1];
				yuv_to_rgb(d1, b, cb, cr);
				s1 += 1;
				d1 += dbpp;
				if (subsample == 1 ||
				   (subsample == 2 && (x & 1)))
					u1 += 2;
			}
			s += stride;
			if (subsample == 1 ||
			   (subsample == 2 && (y & 1)))
				u += 2 * stride / subsample;
			d += dstride;
		}
		break;

	case V4L2_PIX_FMT_YYUV420_V32: {
		static const int VEC_SIZE = (64 + 2*32 + 64) * 2;	/* In bytes */
		static const int LUMA_SHIFT = 8;			/* In theory 8, 7 gives brighter image */
		static const int CHROMA_SHIFT = 6;			/* Should be verified */

		if (width & 63) print(0, "WARNING: width would be better to be multiple of 64\n");
		if (height & 1) error("height must be multiple of 2");

		stride = width;		/*	 Stride on the input buffer is meaningless, so overwrite it */
		s = src = duplicate_buffer(in_buffer, in_size, 2 * width * height * 3/2);
		d = dst = calloc(1, stride * height * 3/2);

		for (y = 0; y < height; y += 2) {
			for (x = 0; x < width; x += 64) {
				int y0, x0, p, c;
				/* Luma */
				for (y0 = 0; y0 < 2; y0++) for (x0 = 0; x0 < 64; x0++) {
					uint16_t *s0 = (uint16_t *)s;
					int x1 = (x0 & 31) | (y0 << 5);
					int y1 = (x0 & 32) >> 5;
					p = s0[y1 * 128 + x1] >> LUMA_SHIFT;
					d[(y+y0)*stride + x+x0] = CLAMPB(p);
				}
				/* Chroma */
				for (x0 = 0; x0 < 32; x0++) for (c = 0; c < 2; c++) {
					int16_t *s0 = (int16_t *)s;
					p = s0[64 + x0 + 32 * c];
					d[stride*height + stride*(y/2) + x + x0*2 + c] =
							CLAMPB((p >> CHROMA_SHIFT) + 128);
				}
				s += VEC_SIZE;
			}
		}
		r = convert(dst, stride * height * 3/2, width, height, stride, V4L2_PIX_FMT_NV12, out_buffer);
		if (r) error("conversion failed");
		break;
	}

	case V4L2_PIX_FMT_GREY:
	case V4L2_PIX_FMT_Y16:
		bpp = (format == V4L2_PIX_FMT_Y16) ? 2 : 1;
		if (stride <= 0) stride = width * bpp;
		s = src = duplicate_buffer(in_buffer, in_size, stride * height);
		for (y = 0; y < height; y++) {
			unsigned char *s1 = s;
			unsigned char *d1 = d;
			for (x = 0; x < width; x++) {
				int b = s1[0];
				if (bpp == 2) {
					b |= s1[1] << 8;
					if (b > 1023) error("Y16 image not in range 0..1023");
					b >>= 2;
				}
				d1[0] = b;
				d1[1] = b;
				d1[2] = b;
				s1 += bpp;
				d1 += dbpp;
			}
			s += stride;
			d += dstride;
		}
		break;

	case V4L2_PIX_FMT_BGR24:
	case V4L2_PIX_FMT_RGB24:
		bpp = 3;
		if (stride <= 0) stride = width * bpp;
		s = src = duplicate_buffer(in_buffer, in_size, stride * height);
		for (y = 0; y < height; y++) {
			unsigned char *s1 = s;
			unsigned char *d1 = d;
			for (x = 0; x < width; x++) {
				if (format == V4L2_PIX_FMT_RGB24) {
					d1[0] = s1[0];
					d1[2] = s1[2];
				} else {
					d1[0] = s1[2];
					d1[2] = s1[0];
				}
				d1[1] = s1[1];
				s1 += bpp;
				d1 += dbpp;
			}
			s += stride;
			d += dstride;
		}
		break;

	case V4L2_PIX_FMT_SBGGR10V32:
	case V4L2_PIX_FMT_SGBRG10V32:
	case V4L2_PIX_FMT_SGRBG10V32:
	case V4L2_PIX_FMT_SRGGB10V32:
	case V4L2_PIX_FMT_SBGGR12V32:
	case V4L2_PIX_FMT_SGBRG12V32:
	case V4L2_PIX_FMT_SGRBG12V32:
	case V4L2_PIX_FMT_SRGGB12V32:
		if (format == V4L2_PIX_FMT_SBGGR10V32) format = V4L2_PIX_FMT_SBGGR10;
		if (format == V4L2_PIX_FMT_SGBRG10V32) format = V4L2_PIX_FMT_SGBRG10;
		if (format == V4L2_PIX_FMT_SGRBG10V32) format = V4L2_PIX_FMT_SGRBG10;
		if (format == V4L2_PIX_FMT_SRGGB10V32) format = V4L2_PIX_FMT_SRGGB10;
		if (format == V4L2_PIX_FMT_SBGGR12V32) format = V4L2_PIX_FMT_SBGGR12;
		if (format == V4L2_PIX_FMT_SGBRG12V32) format = V4L2_PIX_FMT_SGBRG12;
		if (format == V4L2_PIX_FMT_SGRBG12V32) format = V4L2_PIX_FMT_SGRBG12;
		if (format == V4L2_PIX_FMT_SRGGB12V32) format = V4L2_PIX_FMT_SRGGB12;

		if (width & 1) error("must be even image width");
		if (height & 1) error("must be even image height");

		bpp = 2;
		if (stride <= 0) stride = width * bpp;
		s = src = duplicate_buffer(in_buffer, in_size, stride * height);
		d = dst = calloc(1, stride * height);

		for (y = 0; y < height; y += 2) {
			unsigned char *d0 = d;
			unsigned char *s0 = s;
			for (x = 0; x < width; x += 64) {
				for (b = 0; b < 32; b++) {
					if (x + b*2 >= width) break;
					d0[b*4+0] = s0[b*2+0];			// Gr
					d0[b*4+1] = s0[b*2+1];
					d0[b*4+2] = s0[b*2+64];			// R
					d0[b*4+3] = s0[b*2+65];
					d0[stride+b*4+0] = s0[b*2+128];		// B
					d0[stride+b*4+1] = s0[b*2+129];
					d0[stride+b*4+2] = s0[b*2+192];		// Gb
					d0[stride+b*4+3] = s0[b*2+193];
				}
				s0 += 4 * 32 * 2;
				d0 += 64 * bpp;
			}
			d += 2 * stride;
			s += 2 * stride;
		}

		r = convert(dst, stride * height, width, height, stride, format, out_buffer);
		if (r) error("conversion failed");
		break;

	case V4L2_PIX_FMT_SBGGR8:		/* 1 byte per pixel, little endian */
	case V4L2_PIX_FMT_SGBRG8:
	case V4L2_PIX_FMT_SRGGB8:
	case V4L2_PIX_FMT_SGRBG8:
	case V4L2_PIX_FMT_SBGGR10:		/* 2 bytes per pixel, little endian */
	case V4L2_PIX_FMT_SGBRG10:
	case V4L2_PIX_FMT_SRGGB10:
	case V4L2_PIX_FMT_SGRBG10:
	case V4L2_PIX_FMT_SBGGR12:		/* 2 bytes per pixel, little endian */
	case V4L2_PIX_FMT_SGBRG12:
	case V4L2_PIX_FMT_SRGGB12:
	case V4L2_PIX_FMT_SGRBG12:
	case V4L2_PIX_FMT_SBGGR14:		/* 2 bytes per pixel, little endian */
	case V4L2_PIX_FMT_SGBRG14:
	case V4L2_PIX_FMT_SRGGB14:
	case V4L2_PIX_FMT_SGRBG14:
	case V4L2_PIX_FMT_SBGGR16:
		if (format == V4L2_PIX_FMT_SBGGR8 ||
		    format == V4L2_PIX_FMT_SBGGR10 ||
		    format == V4L2_PIX_FMT_SBGGR12 ||
		    format == V4L2_PIX_FMT_SBGGR14 ||
		    format == V4L2_PIX_FMT_SBGGR16) { initrow = 1; initpix = 0; } else
		if (format == V4L2_PIX_FMT_SGBRG8 ||
		    format == V4L2_PIX_FMT_SGBRG10 ||
		    format == V4L2_PIX_FMT_SGBRG12 ||
		    format == V4L2_PIX_FMT_SGBRG14) { initrow = 1; initpix = 1; } else
		if (format == V4L2_PIX_FMT_SRGGB8 ||
		    format == V4L2_PIX_FMT_SRGGB10 ||
		    format == V4L2_PIX_FMT_SRGGB12 ||
		    format == V4L2_PIX_FMT_SRGGB14) { initrow = 0; initpix = 1; } else
		if (format == V4L2_PIX_FMT_SGRBG8 ||
		    format == V4L2_PIX_FMT_SGRBG10 ||
		    format == V4L2_PIX_FMT_SGRBG12 ||
		    format == V4L2_PIX_FMT_SGRBG14) { initrow = 0; initpix = 0; }
		if (format == V4L2_PIX_FMT_SBGGR8 ||
	 	    format == V4L2_PIX_FMT_SGBRG8 ||
		    format == V4L2_PIX_FMT_SRGGB8 ||
		    format == V4L2_PIX_FMT_SGRBG8) {
			bpp = 1;	/* 8 bits per pixel */
			shift = 0;
		} else if (format == V4L2_PIX_FMT_SBGGR10 ||
			   format == V4L2_PIX_FMT_SGBRG10 ||
		 	   format == V4L2_PIX_FMT_SRGGB10 ||
			   format == V4L2_PIX_FMT_SGRBG10) {
			bpp = 2;	/* 10 bits per pixel */
			shift = 2;
		} else if (format == V4L2_PIX_FMT_SBGGR12 ||
			   format == V4L2_PIX_FMT_SGBRG12 ||
		 	   format == V4L2_PIX_FMT_SRGGB12 ||
			   format == V4L2_PIX_FMT_SGRBG12) {
			bpp = 2;
			shift = 4;	/* 12 bits per pixel */
		} else if (format == V4L2_PIX_FMT_SBGGR14 ||
			   format == V4L2_PIX_FMT_SGBRG14 ||
		 	   format == V4L2_PIX_FMT_SRGGB14 ||
			   format == V4L2_PIX_FMT_SGRBG14) {
			bpp = 2;
			shift = 6;	/* 14 bits per pixel */
		} else {
			bpp = 2;
			shift = 8;	/* 16 bits per pixel */
		}

		if (stride <= 0) stride = width * bpp;
		s = src = duplicate_buffer(in_buffer, in_size, stride * height);
		r = g = b = 0;
		oddrow = initrow;
		for (y = 0; y < height; y++) {
			unsigned char *s1 = s;
			unsigned char *d1 = d;
			oddpix = initpix;
			for (x = 0; x < width; x++) {
				if (!oddrow) {
					if (!oddpix) g = get_word(s1, bpp);
					if ( oddpix) r = get_word(s1, bpp);
					if (!oddpix && y > 0) b = get_word(s1 - stride, bpp);
				} else {
					if (!oddpix) b = get_word(s1, bpp);
					if ( oddpix) g = get_word(s1, bpp);
					if ( oddpix && y > 0) r = get_word(s1 - stride, bpp);
				}
				d1[0] = r >> shift;
				d1[1] = g >> shift;
				d1[2] = b >> shift;
				s1 += bpp;
				d1 += dbpp;
				oddpix ^= 1;
			}
			s += stride;
			d += dstride;
			oddrow ^= 1;
		}
		break;
	case V4L2_PIX_FMT_SBGGR10P:		/* 2 bytes per pixel, little endian */
	case V4L2_PIX_FMT_SGBRG10P:
	case V4L2_PIX_FMT_SRGGB10P:
	case V4L2_PIX_FMT_SGRBG10P:
		if (format == V4L2_PIX_FMT_SBGGR10P ||
		    format == V4L2_PIX_FMT_SBGGR12P ||
		    format == V4L2_PIX_FMT_SBGGR14P) { initrow = 1; initpix = 0; } else
		if (format == V4L2_PIX_FMT_SGBRG10P ||
		    format == V4L2_PIX_FMT_SGBRG12P ||
		    format == V4L2_PIX_FMT_SGBRG14P) { initrow = 1; initpix = 1; } else
		if (format == V4L2_PIX_FMT_SRGGB10P ||
		    format == V4L2_PIX_FMT_SRGGB12P ||
		    format == V4L2_PIX_FMT_SRGGB14P) { initrow = 0; initpix = 1; } else
		if (format == V4L2_PIX_FMT_SGRBG10P ||
		    format == V4L2_PIX_FMT_SGRBG12P ||
		    format == V4L2_PIX_FMT_SGRBG14P) { initrow = 0; initpix = 0; }
		bpp = 2;	/* 10 bits per pixel */
		shift = 2;
		if (stride <= 0) stride = width / 4 * 5;
		printf("stride %d\n", stride);
		s = src = duplicate_buffer(in_buffer, in_size, stride * height);
		r = g = b = r1 = g1 = b1 = 0;
		oddrow = initrow;
		for (y = 0; y < height; y++) {
			unsigned char *s1 = s;
			unsigned char *d1 = d;
			oddpix = initpix;
			for (x = 0; x < width; x++) {
				/* packed raw sRGB / Bayer formats with 10 bits
				per sample. Every four consecutive samples are
				packed into 5 bytes. Each of the first 4 bytes
				contain the 8 high order bits of the pixels, and
				the 5th byte contains the 2 least significants
				bits of each pixel, in the same order.
				*/
				if (!oddrow) {
					if (!oddpix) g = get_packed_word(s1, x);
					if ( oddpix) r = get_packed_word(s1, x);
					if (!oddpix && y > 0) b = get_packed_word(s1 - stride, x);
				} else {
					if (!oddpix) b = get_packed_word(s1, x);
					if ( oddpix) g = get_packed_word(s1, x);
					if ( oddpix && y > 0) r = get_packed_word(s1 - stride, x);
				}
				//printf("%d,%d 0x%.4x 0x%.4x 0x%.4x %d\n", y,x,r,g,b, oddpix);

				d1[0] = r >> shift;
				d1[1] = g >> shift;
				d1[2] = b >> shift;
				s1 += 1;
				/* jump over every 5th byte to skip the low part data */
				if (x % 4 == 3)
					s1 += 1;
				d1 += dbpp;
				oddpix ^= 1;
				//if (x > 10)
				//	return -1;
			}
			//printf("s1 %d vs. %d\n", width, s1-s);
			s += stride;
			d += dstride;
			oddrow ^= 1;
		}
		break;

	default:
		errno = EINVAL;
		return -1;
	}

	free(src);
	free(dst);
	return 0;
}

int main(int argc, char *argv[])
{
	char *in_name = NULL;
	char *out_name = NULL;
	void *in_buffer;
	void *out_buffer;
	int in_size, out_size;
	FILE *f;
	int i;

	int opt;
	int width = -1;
	int height = -1;
	int stride = -1;
	int skip = 0;
	__u32 format = 0;

	while ((opt = getopt(argc, argv, "hf:x:y:s:b:")) != -1) {
		switch (opt) {
		case 'f': {
			const char *t = optarg;
			format = symbol_get(pixelformats, &t);
			break;
		}
		case 'x':
			width = atoi(optarg);
			break;
		case 'y':
			height = atoi(optarg);
			break;
		case 's':
			stride = atoi(optarg);
			break;
		case 'b':
			skip = atoi(optarg);
			break;
		default:
			usage();
			print(1, "Available formats:\n");
			symbol_dump(V4L2_PIX_FMT, pixelformats);
			return -1;
	        }
	}

	if (optind+2 > argc)
		error("input or output filename missing");

	in_name = argv[optind++];
	out_name = argv[optind++];

	print(1, "Reading file `%s', %ix%i stride %i format %s, skip %i\n",
		in_name, width, height, stride, symbol_str(format, pixelformats), skip);
	f = fopen(in_name, "rb");
	if (!f) error("failed opening file");
	in_size = fsize(f);
	if (in_size == -1) error("error checking file size");
	print(2, "File size %i bytes, data size %i bytes\n", in_size, in_size - skip);
	if (skip >= in_size) error("no data left to read");
	if (fseek(f, skip, SEEK_SET) != 0) error("seeking failed");
	in_size -= skip;
	in_buffer = malloc(in_size);
	if (!in_buffer) error("can not allocate input buffer");
	i = fread(in_buffer, in_size, 1, f);
	if (i != 1) error("failed reading file");
	fclose(f);

	out_size = width * height * 3;
	out_buffer = malloc(out_size);
	if (!out_buffer) error("can not allocate output buffer");
	i = convert(in_buffer, in_size, width, height, stride, format, out_buffer);
	if (i < 0) error("failed to convert image");
	free(in_buffer);

	print(1, "Writing file `%s', %i bytes\n", out_name, out_size);
	f = fopen(out_name, "wb");
	if (!f) error("failed opening file");
	i = fprintf(f, "P6\n%i %i 255\n", width, height);
	if (i < 0) error("can not write file header");
	i = fwrite(out_buffer, out_size, 1, f);
	if (i != 1) error("failed writing file");
	fclose(f);

	free(out_buffer);
	return 0;
}
